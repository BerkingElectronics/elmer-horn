"""
Elmer-Horn MCP-Server
=====================

Lokaler MCP-Server, der Claude erlaubt, rotationssymmetrische Hörner in
Elmer FEM zu simulieren:

    Geometrie (Gmsh .geo)  ->  Netz (Gmsh -> ElmerGrid)  ->
    Frequenz-Sweep (ElmerSolver, Helmholtz, achsensymmetrisch)  ->
    Auswertung (SPL unter Winkeln, Halsimpedanz)  ->  Plots

Modell (Version 1):
  * Horn in unendlicher Schallwand (Halbraum), achsensymmetrisch
  * Hals: Kolben mit konstanter Schnelle v0 = 1 m/s (Pegel sind daher
    "dB re 1 m/s Halsschnelle", für Vergleiche zwischen Varianten gedacht)
  * Hornwand + Schallwand: schallhart
  * Außenrand (Viertelkreis, Radius R): absorbierend (Plane Wave BC)

Alle Projekte liegen in HORN_WORKDIR (Standard: ~/ElmerHorn).
"""

from __future__ import annotations

import cmath
import json
import math
import os
import re
import shutil
import subprocess
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

try:  # MCP-SDK 1.x
    from mcp.server.fastmcp import FastMCP, Image
except ImportError:  # MCP-SDK 2.x (FastMCP heißt dort MCPServer)
    from mcp.server.mcpserver import MCPServer as FastMCP, Image

# --------------------------------------------------------------------------
# Konstanten / Konfiguration
# --------------------------------------------------------------------------

RHO = 1.205        # Luftdichte kg/m^3 (20 °C)
C0 = 343.0         # Schallgeschwindigkeit m/s
V0 = 1.0           # Halsschnelle m/s
P_REF = 20e-6      # Bezugsschalldruck Pa
POLAR_ANGLES = [0, 15, 30, 45, 60, 75, 90]

# Tags der Physical Groups in Gmsh = Boundary-Nummern in Elmer
BND_THROAT, BND_WALL, BND_BAFFLE, BND_OUTER, BND_AXIS = 1, 2, 3, 4, 5

WORKDIR = Path(os.environ.get("HORN_WORKDIR", str(Path.home() / "ElmerHorn"))).expanduser().resolve()
WORKDIR.mkdir(parents=True, exist_ok=True)

mcp = FastMCP("elmer-horn")


# --------------------------------------------------------------------------
# Hilfsfunktionen: Programme finden
# --------------------------------------------------------------------------

def _exe_name(name: str) -> str:
    return name + ".exe" if os.name == "nt" else name


def _version_key(p: Path) -> tuple:
    """Sortiert 'Elmer 26.2-Release' vor 'Elmer 9.0-Release' (numerisch statt alphabetisch)."""
    return tuple(int(n) for n in re.findall(r"\d+", p.name)) or (0,)


def _find_exe(name: str, env_var: str) -> str | None:
    """Sucht ein Programm über Umgebungsvariable, PATH, ELMER_HOME und
    typische Installationsordner unter Windows."""
    cands: list[Path] = []
    if os.environ.get(env_var):
        cands.append(Path(os.environ[env_var]))
    w = shutil.which(name)
    if w:
        cands.append(Path(w))
    exe = _exe_name(name)
    if os.environ.get("ELMER_HOME"):
        cands.append(Path(os.environ["ELMER_HOME"]) / "bin" / exe)
    if os.name == "nt":
        bases = [os.environ.get("ProgramFiles"), os.environ.get("ProgramFiles(x86)"),
                 os.environ.get("LOCALAPPDATA"), "C:\\"]
        for base in [b for b in bases if b]:
            bp = Path(base)
            if not bp.exists():
                continue
            for d in sorted(bp.glob("Elmer*"), key=_version_key, reverse=True):
                cands.append(d / "bin" / exe)
            for d in sorted(bp.glob("gmsh*"), key=_version_key, reverse=True):
                cands.append(d / exe)
                cands.append(d / "bin" / exe)
    for c in cands:
        try:
            if c.is_file():
                return str(c)
        except OSError:
            pass
    return None


def _tools() -> dict[str, str | None]:
    return {
        "ElmerSolver": _find_exe("ElmerSolver", "ELMER_SOLVER"),
        "ElmerGrid": _find_exe("ElmerGrid", "ELMER_GRID"),
        "gmsh": _find_exe("gmsh", "GMSH_EXE"),
    }


def _require(tool: str) -> str:
    p = _tools()[tool]
    if not p:
        raise RuntimeError(
            f"{tool} wurde nicht gefunden. Bitte zum PATH hinzufügen oder die "
            f"Umgebungsvariable {'GMSH_EXE' if tool == 'gmsh' else 'ELMER_SOLVER/ELMER_GRID'} "
            f"in der Connector-Konfiguration setzen.")
    return p


def _run(cmd: list[str], cwd: Path, timeout: int = 3600) -> subprocess.CompletedProcess:
    env = dict(os.environ)
    env.setdefault("OMP_NUM_THREADS", "1")
    kwargs = {}
    if os.name == "nt":
        kwargs["creationflags"] = 0x08000000  # CREATE_NO_WINDOW
    return subprocess.run(cmd, cwd=str(cwd), capture_output=True, text=True,
                          timeout=timeout, env=env, errors="replace", **kwargs)


def _tail(text: str, n: int = 40) -> str:
    lines = (text or "").strip().splitlines()
    return "\n".join(lines[-n:])


# --------------------------------------------------------------------------
# Projekte
# --------------------------------------------------------------------------

def _proj(name: str) -> Path:
    if not re.fullmatch(r"[A-Za-z0-9_\-]{1,64}", name or ""):
        raise ValueError("Projektname nur aus Buchstaben, Ziffern, _ und - (max. 64 Zeichen).")
    return WORKDIR / name


def _safe_file(proj: Path, filename: str) -> Path:
    p = (proj / filename).resolve()
    if proj.resolve() not in p.parents and p != proj.resolve():
        raise ValueError("Datei muss innerhalb des Projektordners liegen.")
    return p


def _load_cfg(proj: Path) -> dict:
    f = proj / "project.json"
    if not f.exists():
        raise FileNotFoundError(f"Projekt '{proj.name}' existiert nicht (erst create_horn ausführen).")
    return json.loads(f.read_text(encoding="utf-8"))


# --------------------------------------------------------------------------
# Hornprofile
# --------------------------------------------------------------------------

def horn_profile(profile: str, throat_r: float, length: float,
                 mouth_r: float | None = None, coverage_deg: float | None = None,
                 points: list[list[float]] | None = None, n: int = 40) -> list[tuple[float, float]]:
    """Liefert Liste (r, z) in Metern. Hals bei z = -length, Mündung bei z = 0."""
    profile = profile.lower()
    pts: list[tuple[float, float]] = []
    if profile == "custom":
        if not points or len(points) < 2:
            raise ValueError("Für profile='custom' mindestens 2 Punkte [z_mm, r_mm] angeben.")
        raw = sorted((float(z) / 1000, float(r) / 1000) for z, r in points)
        z0, zL = raw[0][0], raw[-1][0]
        if zL <= z0:
            raise ValueError("Profil braucht steigende z-Werte.")
        return [(r, z - zL) for z, r in raw]
    m_hyp = None
    if profile == "hyperbolic":   # Hypex nach Salmon, T = 0.7; m per Bisektion
        T = 0.7
        lo, hi = 1e-9, 200.0 / length
        for _ in range(200):
            mid = 0.5 * (lo + hi)
            if throat_r * (math.cosh(mid * length) + T * math.sinh(mid * length)) > mouth_r:
                hi = mid
            else:
                lo = mid
        m_hyp = 0.5 * (lo + hi)
    for i in range(n + 1):
        t = i / n
        z = t * length
        if profile == "conical":
            r = throat_r + (mouth_r - throat_r) * t
        elif profile == "exponential":
            r = throat_r * (mouth_r / throat_r) ** t
        elif profile in ("os", "oblate_spheroidal"):
            a = math.radians(coverage_deg / 2.0)
            r = math.sqrt(throat_r ** 2 + (z * math.tan(a)) ** 2)
        elif profile == "hyperbolic":
            r = throat_r * (math.cosh(m_hyp * z) + 0.7 * math.sinh(m_hyp * z))
        else:
            raise ValueError("Unbekanntes Profil. Erlaubt: conical, exponential, os, hyperbolic, custom")
        pts.append((r, z - length))
    return pts


def write_geo(path: Path, prof: list[tuple[float, float]], R: float, h: float, h_throat: float) -> None:
    """Schreibt das Gmsh-Skript. Koordinaten: x = Radius r, y = Achse z."""
    L = -prof[0][1]
    lines = ["// automatisch erzeugt vom Elmer-Horn-Connector",
             "Mesh.MshFileVersion = 2.2;",
             "Mesh.Algorithm = 6;",
             f"Mesh.CharacteristicLengthMax = {h:.6g};"]
    pid = 1
    lines.append(f"Point({pid}) = {{0, {-L:.9g}, 0, {h_throat:.6g}}};")  # Achse am Hals
    axis_throat = pid
    prof_ids = []
    nprof = len(prof)
    for i, (r, z) in enumerate(prof):
        pid += 1
        # Netzgröße entlang des Horns vom Hals bis zur Mündung verlaufen lassen
        hs = h_throat + (h - h_throat) * (i / max(nprof - 1, 1))
        lines.append(f"Point({pid}) = {{{r:.9g}, {z:.9g}, 0, {hs:.6g}}};")
        prof_ids.append(pid)
    pid += 1; p_baffle = pid
    lines.append(f"Point({pid}) = {{{R:.9g}, 0, 0, {h:.6g}}};")
    pid += 1; p_center = pid
    lines.append(f"Point({pid}) = {{0, 0, 0, {h:.6g}}};")
    pid += 1; p_top = pid
    lines.append(f"Point({pid}) = {{0, {R:.9g}, 0, {h:.6g}}};")

    lines += [
        f"Line(1) = {{{axis_throat}, {prof_ids[0]}}};",                        # Hals
        f"Spline(2) = {{{', '.join(map(str, prof_ids))}}};",                   # Hornwand
        f"Line(3) = {{{prof_ids[-1]}, {p_baffle}}};",                          # Schallwand
        f"Circle(4) = {{{p_baffle}, {p_center}, {p_top}}};",                   # Außenrand
        f"Line(5) = {{{p_top}, {axis_throat}}};",                              # Achse
        "Curve Loop(1) = {1, 2, 3, 4, 5};",
        "Plane Surface(1) = {1};",
        f"Physical Curve({BND_THROAT}) = {{1}};",
        f"Physical Curve({BND_WALL}) = {{2}};",
        f"Physical Curve({BND_BAFFLE}) = {{3}};",
        f"Physical Curve({BND_OUTER}) = {{4}};",
        f"Physical Curve({BND_AXIS}) = {{5}};",
        "Physical Surface(1) = {1};",
    ]
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


# --------------------------------------------------------------------------
# Elmer-Solver-Input
# --------------------------------------------------------------------------

def write_sif(path: Path, f: float, tag: str, direct_method: str = "umfpack",
              absorber: str = "planewave") -> None:
    w = 2 * math.pi * f
    flux = w * RHO * V0   # |dp/dn| = omega * rho * v  (Euler-Gleichung)
    if absorber == "planewave":
        outer = "  Plane Wave BC = Logical True"
    elif absorber == "impedance":
        # spezifische Impedanz rho*c -> Elmer erwartet hier die Impedanz Z
        outer = f"  Wave Impedance 1 = Real {RHO * C0:.6g}\n  Wave Impedance 2 = Real 0.0"
    else:
        raise ValueError("absorber muss 'planewave' oder 'impedance' sein")
    sif = f"""! automatisch erzeugt vom Elmer-Horn-Connector, f = {f:.6g} Hz
Header
  CHECK KEYWORDS Warn
  Mesh DB "." "mesh"
End

Simulation
  Max Output Level = 3
  Coordinate System = "Axi Symmetric"
  Simulation Type = Steady state
  Steady State Max Iterations = 1
  Frequency = Real {f:.9g}
  Angular Frequency = Real {w:.9g}
End

Constants
End

Body 1
  Target Bodies(1) = 1
  Equation = 1
  Material = 1
End

Material 1
  Sound Speed = Real {C0}
  Density = Real {RHO}
End

Equation 1
  Active Solvers(1) = 1
End

Solver 1
  Equation = Helmholtz
  Procedure = "HelmholtzSolve" "HelmholtzSolver"
  Variable = -dofs 2 Pressure Wave
  Frequency = Real {f:.9g}
  Angular Frequency = Real {w:.9g}
  Linear System Solver = Direct
  Linear System Direct Method = {direct_method}
End

Solver 2
  Exec Solver = After Simulation
  Equation = SaveLine
  Procedure = "SaveData" "SaveLine"
  Filename = "f_{tag}.dat"
End

! Hals: Kolbenquelle mit Schnelle {V0} m/s
Boundary Condition 1
  Target Boundaries(1) = {BND_THROAT}
  Wave Flux 1 = Real 0.0
  Wave Flux 2 = Real {flux:.9g}
  Save Line = Logical True
End

! Hornwand + Schallwand: schallhart (natürliche RB, nichts zu setzen)
Boundary Condition 2
  Target Boundaries(2) = {BND_WALL} {BND_BAFFLE}
End

! Außenrand: absorbierend
Boundary Condition 3
  Target Boundaries(1) = {BND_OUTER}
{outer}
  Save Line = Logical True
End
"""
    path.write_text(sif, encoding="utf-8")


# --------------------------------------------------------------------------
# Ergebnisse einlesen / auswerten
# --------------------------------------------------------------------------

def read_saveline(dat: Path) -> list[dict]:
    """Liest eine SaveLine-Datei samt .names-Datei. Gibt Zeilen mit x, y, p (komplex) zurück."""
    names_file = Path(str(dat) + ".names")
    cols: dict[int, str] = {}
    if names_file.exists():
        for line in names_file.read_text(errors="replace").splitlines():
            m = re.match(r"^\s*(\d+)\s*:\s*(.+?)\s*$", line)
            if m:
                cols[int(m.group(1)) - 1] = m.group(2).lower()

    def find(pred, default=None):
        for i, n in cols.items():
            if pred(n):
                return i
        return default

    ix = find(lambda n: n.startswith("coordinate 1"))
    iy = find(lambda n: n.startswith("coordinate 2"))
    ip1 = find(lambda n: "pressure" in n and n.rstrip().endswith("1") and "coordinate" not in n)
    ip2 = find(lambda n: "pressure" in n and n.rstrip().endswith("2") and "coordinate" not in n)
    rows = []
    for line in dat.read_text(errors="replace").splitlines():
        parts = line.split()
        if not parts:
            continue
        try:
            vals = [float(v) for v in parts]
        except ValueError:
            continue
        if None in (ix, iy, ip1, ip2):
            # Fallback: letzte 5 Spalten = x y z p1 p2
            x, y, p1, p2 = vals[-5], vals[-4], vals[-2], vals[-1]
        else:
            x, y, p1, p2 = vals[ix], vals[iy], vals[ip1], vals[ip2]
        rows.append({"x": x, "y": y, "p": complex(p1, p2)})
    return rows


def _interp(xs: list[float], ys: list[complex], x: float) -> complex:
    if x <= xs[0]:
        return ys[0]
    if x >= xs[-1]:
        return ys[-1]
    for i in range(1, len(xs)):
        if xs[i] >= x:
            t = (x - xs[i - 1]) / (xs[i] - xs[i - 1]) if xs[i] != xs[i - 1] else 0
            return ys[i - 1] + (ys[i] - ys[i - 1]) * t
    return ys[-1]


def analyse(rows: list[dict], L: float, R: float) -> dict:
    throat = [r for r in rows if abs(r["y"] + L) < 1e-6 + 1e-5 * L]
    outer = [r for r in rows if abs(math.hypot(r["x"], r["y"]) - R) < 2e-3 * R and r["y"] >= -1e-9]
    res: dict = {}
    # Halsimpedanz (flächengewichteter Mittelwert, Achsensymmetrie -> Gewicht r)
    if len(throat) >= 2:
        th = sorted({round(r["x"], 12): r for r in throat}.values(), key=lambda r: r["x"])
        num = 0j
        den = 0.0
        for a, b in zip(th, th[1:]):
            dr = b["x"] - a["x"]
            num += 0.5 * (a["p"] * a["x"] + b["p"] * b["x"]) * dr
            den += 0.5 * (a["x"] + b["x"]) * dr
        p_avg = num / den if den else th[0]["p"]
        z = p_avg / (RHO * C0 * V0)
        note = ""
        if z.real < 0:   # Strahlungswiderstand muss positiv sein -> Vorzeichenkonvention
            z = -z
            note = "sign_flipped"
        res["z_throat"] = [z.real, z.imag]
        if note:
            res["note"] = note
    # SPL unter Winkeln (0° = Achse, 90° = entlang der Schallwand)
    if len(outer) >= 2:
        ou = sorted(outer, key=lambda r: math.degrees(math.atan2(r["x"], r["y"])))
        th = [math.degrees(math.atan2(r["x"], r["y"])) for r in ou]
        ps = [r["p"] for r in ou]
        spl = {}
        for a in POLAR_ANGLES:
            p = _interp(th, ps, a)
            prms = abs(p) / math.sqrt(2)
            spl[str(a)] = 20 * math.log10(max(prms, 1e-30) / P_REF) + 20 * math.log10(R / 1.0)
        res["spl_1m"] = spl
    res["n_throat_nodes"] = len(throat)
    res["n_outer_nodes"] = len(outer)
    return res


def _find_output(proj: Path, fname: str) -> Path | None:
    hits = sorted(proj.rglob(fname), key=lambda p: p.stat().st_mtime, reverse=True)
    return hits[0] if hits else None


def _freq_list(f_start: float, f_stop: float, ppo: int) -> list[float]:
    n = int(math.floor(math.log2(f_stop / f_start) * ppo + 1e-9))
    fs = [f_start * 2 ** (i / ppo) for i in range(n + 1)]
    if fs[-1] < f_stop * 0.999:
        fs.append(f_stop)
    return [round(f, 2) for f in fs]


# --------------------------------------------------------------------------
# MCP-Werkzeuge
# --------------------------------------------------------------------------

@mcp.tool()
def check_setup() -> str:
    """Prüft, ob ElmerSolver, ElmerGrid und gmsh gefunden werden, und zeigt
    Arbeitsordner und CPU-Kerne. Immer zuerst ausführen."""
    t = _tools()
    out = [f"Arbeitsordner: {WORKDIR}", f"CPU-Kerne (logisch): {os.cpu_count()}"]
    for k, v in t.items():
        out.append(f"{k}: {v or 'NICHT GEFUNDEN'}")
    if t["gmsh"]:
        try:
            r = _run([t["gmsh"], "--version"], WORKDIR, timeout=30)
            out.append(f"gmsh-Version: {(r.stdout or r.stderr).strip()}")
        except Exception as e:  # noqa: BLE001
            out.append(f"gmsh-Version: Fehler {e}")
    try:
        import matplotlib  # noqa: F401
        out.append("matplotlib: ok (Plots möglich)")
    except ImportError:
        out.append("matplotlib: fehlt (pip install matplotlib für Plots)")
    return "\n".join(out)


@mcp.tool()
def create_horn(name: str, profile: str, throat_diameter_mm: float, length_mm: float,
                mouth_diameter_mm: float | None = None, coverage_angle_deg: float | None = None,
                profile_points_mm: list[list[float]] | None = None,
                f_max_hz: float = 5000.0, domain_radius_m: float | None = None,
                elements_per_wavelength: float = 10.0) -> str:
    """Erzeugt ein rotationssymmetrisches Horn in unendlicher Schallwand und vernetzt es.

    profile: 'conical' | 'exponential' | 'hyperbolic' (brauchen mouth_diameter_mm),
             'os' (oblate spheroidal, braucht coverage_angle_deg = voller Abstrahlwinkel),
             'custom' (profile_points_mm = [[z_mm, r_mm], ...], z vom Hals = 0 zur Mündung).
    f_max_hz: höchste Frequenz, für die das Netz fein genug sein soll.
    domain_radius_m: Radius des Luftraums vor der Schallwand (Standard: automatisch).
    """
    proj = _proj(name)
    proj.mkdir(parents=True, exist_ok=True)
    rt = throat_diameter_mm / 2000.0
    L = length_mm / 1000.0
    rm = mouth_diameter_mm / 2000.0 if mouth_diameter_mm else None
    p = profile.lower()
    if p in ("conical", "exponential", "hyperbolic") and not rm:
        raise ValueError("mouth_diameter_mm fehlt.")
    if p in ("os", "oblate_spheroidal") and not coverage_angle_deg:
        raise ValueError("coverage_angle_deg fehlt (voller Winkel, z. B. 90).")
    prof = horn_profile(p, rt, L, rm, coverage_angle_deg, profile_points_mm)
    rt, L = prof[0][0], -prof[0][1]
    rm = prof[-1][0]

    lam_min = C0 / f_max_hz
    h = lam_min / elements_per_wavelength
    h_throat = min(h, rt / 4)
    R = domain_radius_m or max(4 * rm, 0.6, 2 * L)
    if R <= rm * 1.2:
        raise ValueError(f"Luftraum-Radius {R} m zu klein für Mündungsradius {rm} m.")

    write_geo(proj / "mesh.geo", prof, R, h, h_throat)
    gmsh = _require("gmsh")
    grid = _require("ElmerGrid")
    r1 = _run([gmsh, "mesh.geo", "-2", "-format", "msh22", "-o", "mesh.msh"], proj, timeout=900)
    if r1.returncode != 0 or not (proj / "mesh.msh").exists():
        return "Gmsh-Fehler:\n" + _tail(r1.stdout + "\n" + r1.stderr)
    if (proj / "mesh").exists():
        shutil.rmtree(proj / "mesh")
    r2 = _run([grid, "14", "2", "mesh.msh", "-autoclean"], proj, timeout=900)
    header = proj / "mesh" / "mesh.header"
    if not header.exists():
        return "ElmerGrid-Fehler:\n" + _tail(r2.stdout + "\n" + r2.stderr)
    nn, ne, nb = (header.read_text().split() + ["?", "?", "?"])[:3]

    cfg = {"name": name, "profile": p, "throat_r": rt, "mouth_r": rm, "length": L,
           "coverage_deg": coverage_angle_deg, "R": R, "h": h, "f_max_mesh": f_max_hz,
           "profile_points": prof, "created": time.strftime("%Y-%m-%d %H:%M:%S")}
    (proj / "project.json").write_text(json.dumps(cfg, indent=2), encoding="utf-8")
    fc = C0 / (2 * math.pi * rm) if rm else float("nan")
    return (f"Projekt '{name}' angelegt in {proj}\n"
            f"Profil {p}: Hals Ø {2*rt*1000:.1f} mm, Mündung Ø {2*rm*1000:.1f} mm, Länge {L*1000:.1f} mm\n"
            f"Luftraum R = {R:.3f} m, Elementgröße {h*1000:.1f} mm (Hals {h_throat*1000:.1f} mm)\n"
            f"Netz: {nn} Knoten, {ne} Elemente, {nb} Randelemente\n"
            f"Faustwert Mündung (ka = 1): ca. {fc:.0f} Hz")


@mcp.tool()
def run_sweep(name: str, f_start_hz: float = 200.0, f_stop_hz: float = 5000.0,
              points_per_octave: int = 6, workers: int | None = None,
              direct_method: str = "umfpack", absorber: str = "planewave") -> str:
    """Rechnet einen Frequenz-Sweep. Jede Frequenz läuft als eigener ElmerSolver-Prozess,
    `workers` davon parallel (Standard: Hälfte der logischen Kerne).
    direct_method: umfpack | mumps | banded.  absorber: planewave | impedance."""
    proj = _proj(name)
    cfg = _load_cfg(proj)
    solver = _require("ElmerSolver")
    freqs = _freq_list(f_start_hz, f_stop_hz, points_per_octave)
    warn = ""
    if f_stop_hz > cfg["f_max_mesh"] * 1.05:
        warn = (f"WARNUNG: Netz nur bis {cfg['f_max_mesh']:.0f} Hz ausgelegt, "
                f"Ergebnisse darüber sind unzuverlässig.\n")
    (proj / "runs").mkdir(exist_ok=True)
    workers = workers or max(1, (os.cpu_count() or 2) // 2)

    def one(f: float) -> tuple[float, dict]:
        tag = f"{f:08.2f}".replace(".", "_")
        sif = proj / "runs" / f"f_{tag}.sif"
        write_sif(sif, f, tag, direct_method, absorber)
        rel = sif.relative_to(proj).as_posix()
        t0 = time.time()
        r = _run([solver, rel], proj, timeout=7200)
        (proj / "runs" / f"f_{tag}.log").write_text(r.stdout + "\n" + r.stderr, encoding="utf-8")
        out = _find_output(proj, f"f_{tag}.dat")
        if r.returncode != 0 or out is None:
            return f, {"error": _tail(r.stdout + "\n" + r.stderr, 25), "returncode": r.returncode}
        a = analyse(read_saveline(out), cfg["length"], cfg["R"])
        a["seconds"] = round(time.time() - t0, 1)
        return f, a

    results: dict[str, dict] = {}
    errors = []
    t_start = time.time()
    with ThreadPoolExecutor(max_workers=workers) as ex:
        futs = [ex.submit(one, f) for f in freqs]
        for fu in as_completed(futs):
            f, a = fu.result()
            results[f"{f:.2f}"] = a
            if "error" in a:
                errors.append((f, a["error"]))
    res_file = proj / "results.json"
    res_file.write_text(json.dumps({"freqs": freqs, "data": results}, indent=2), encoding="utf-8")
    ok = len(freqs) - len(errors)
    msg = (f"{warn}Sweep fertig: {ok}/{len(freqs)} Frequenzen erfolgreich, "
           f"{time.time() - t_start:.0f} s mit {workers} parallelen Prozessen.\n"
           f"Ergebnisse: {res_file}")
    if errors:
        f, e = sorted(errors)[0]
        msg += f"\n\nErster Fehler bei {f} Hz (Log in runs/):\n{e}"
    return msg


@mcp.tool()
def get_results(name: str) -> str:
    """Tabelle der Ergebnisse: SPL (dB re 1 m/s Halsschnelle, auf 1 m umgerechnet) unter
    0–90° sowie normierte Halsimpedanz (Re/Im, bezogen auf rho*c). Re ≈ 1 = gute Hornbeladung."""
    proj = _proj(name)
    data = json.loads((proj / "results.json").read_text(encoding="utf-8"))
    head = "f [Hz] | " + " | ".join(f"{a}°" for a in POLAR_ANGLES) + " | Re(Z) | Im(Z)"
    lines = [head, "-" * len(head)]
    for f in data["freqs"]:
        d = data["data"].get(f"{f:.2f}", {})
        if "error" in d or "spl_1m" not in d:
            lines.append(f"{f:7.0f} | FEHLER")
            continue
        spl = " | ".join(f"{d['spl_1m'][str(a)]:5.1f}" for a in POLAR_ANGLES)
        z = d.get("z_throat", [float("nan")] * 2)
        lines.append(f"{f:7.0f} | {spl} | {z[0]:5.2f} | {z[1]:5.2f}")
    if any(d.get("note") == "sign_flipped" for d in data["data"].values()):
        lines.append("\nHinweis: Vorzeichen der Impedanz wurde korrigiert (Konvention von Elmer).")
    return "\n".join(lines)


@mcp.tool()
def plot_results(name: str) -> Image:
    """Erzeugt ein PNG mit Frequenzgängen unter Winkeln, normierter Richtwirkung und
    Halsimpedanz und gibt es als Bild zurück (Datei: <projekt>/results.png)."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    proj = _proj(name)
    data = json.loads((proj / "results.json").read_text(encoding="utf-8"))
    fs, spl, zr, zi = [], {a: [] for a in POLAR_ANGLES}, [], []
    for f in data["freqs"]:
        d = data["data"].get(f"{f:.2f}", {})
        if "spl_1m" not in d:
            continue
        fs.append(f)
        for a in POLAR_ANGLES:
            spl[a].append(d["spl_1m"][str(a)])
        z = d.get("z_throat", [float("nan")] * 2)
        zr.append(z[0]); zi.append(z[1])
    if not fs:
        raise RuntimeError("Keine gültigen Ergebnisse zum Plotten.")

    cmap = plt.get_cmap("viridis")
    fig, ax = plt.subplots(3, 1, figsize=(9, 11), sharex=True)
    for i, a in enumerate(POLAR_ANGLES):
        col = cmap(i / (len(POLAR_ANGLES) - 1))
        ax[0].semilogx(fs, spl[a], color=col, label=f"{a}°")
        norm = [s - s0 for s, s0 in zip(spl[a], spl[0])]
        ax[1].semilogx(fs, norm, color=col, label=f"{a}°")
    ax[0].set_ylabel("SPL @1 m [dB]\n(re 1 m/s Halsschnelle)")
    ax[0].legend(ncol=4, fontsize=8)
    ax[1].set_ylabel("relativ zu 0° [dB]")
    ax[1].set_ylim(-30, 5)
    ax[2].semilogx(fs, zr, color="#2b6cb0", label="Re(Z)/ρc")
    ax[2].semilogx(fs, zi, color="#718096", ls="--", label="Im(Z)/ρc")
    ax[2].axhline(1.0, color="#a0aec0", lw=0.8)
    ax[2].set_ylabel("Halsimpedanz (normiert)")
    ax[2].set_xlabel("Frequenz [Hz]")
    ax[2].legend(fontsize=8)
    for a_ in ax:
        a_.grid(True, which="both", alpha=0.3)
    fig.suptitle(f"Horn '{name}'")
    fig.tight_layout()
    out = proj / "results.png"
    fig.savefig(out, dpi=110)
    plt.close(fig)
    return Image(path=str(out))


@mcp.tool()
def plot_field(name: str, freqs_hz: list[float], direct_method: str = "umfpack") -> Image:
    """Rechnet das komplette Schallfeld für 1–4 Frequenzen und rendert pro Frequenz
    Pegelbild (SPL) und Wellenbild (Momentandruck), an der Achse gespiegelt.
    Die VTU-Dateien bleiben im Ordner <projekt>/mesh und lassen sich in ParaView öffnen.
    Bild: <projekt>/field.png"""
    from fieldplot import plot_fields

    proj = _proj(name)
    cfg = _load_cfg(proj)
    solver = _require("ElmerSolver")
    freqs = [float(f) for f in freqs_hz][:4]
    (proj / "runs").mkdir(exist_ok=True)
    files = []
    for f in freqs:
        tag = f"{int(round(f)):05d}"
        sif = proj / "runs" / f"field_{tag}.sif"
        write_sif(sif, f, "field_" + tag, direct_method)
        txt = sif.read_text(encoding="utf-8")
        # SaveLine durch VTU-Ausgabe ersetzen
        txt = re.sub(r"Solver 2\n.*?\nEnd\n",
                     "Solver 2\n  Exec Solver = After Simulation\n  Equation = ResultOutput\n"
                     "  Procedure = \"ResultOutputSolve\" \"ResultOutputSolver\"\n"
                     f"  Output File Name = \"field_{tag}\"\n  Output Format = String \"vtu\"\n"
                     "  Ascii Output = Logical True\nEnd\n", txt, flags=re.S)
        sif.write_text(txt, encoding="utf-8")
        r = _run([solver, sif.relative_to(proj).as_posix()], proj, timeout=7200)
        out = _find_output(proj, f"field_{tag}_t0001.vtu")
        if r.returncode != 0 or out is None:
            raise RuntimeError(f"Elmer-Fehler bei {f} Hz:\n" + _tail(r.stdout + "\n" + r.stderr, 25))
        files.append((f, out))
    png = plot_fields(files, proj / "project.json", proj / "field.png", f"Horn '{name}'")
    return Image(path=str(png))


@mcp.tool()
def list_projects() -> str:
    """Listet alle Hornprojekte im Arbeitsordner."""
    out = []
    for p in sorted(WORKDIR.iterdir()):
        if (p / "project.json").exists():
            c = json.loads((p / "project.json").read_text(encoding="utf-8"))
            has = "mit Ergebnissen" if (p / "results.json").exists() else "ohne Ergebnisse"
            out.append(f"{p.name}: {c['profile']}, Hals Ø {2000*c['throat_r']:.0f} mm, "
                       f"Mündung Ø {2000*c['mouth_r']:.0f} mm, L {1000*c['length']:.0f} mm ({has})")
    return "\n".join(out) or "Noch keine Projekte."


@mcp.tool()
def read_project_file(name: str, filename: str, max_chars: int = 20000) -> str:
    """Liest eine Textdatei aus einem Projekt (z. B. mesh.geo, runs/f_01000_00.log, project.json)."""
    proj = _proj(name)
    p = _safe_file(proj, filename)
    txt = p.read_text(encoding="utf-8", errors="replace")
    return txt if len(txt) <= max_chars else txt[:max_chars] + f"\n... ({len(txt)} Zeichen gesamt)"


@mcp.tool()
def write_project_file(name: str, filename: str, content: str) -> str:
    """Schreibt eine Textdatei in ein Projekt (z. B. eigene .sif oder geänderte mesh.geo)."""
    proj = _proj(name)
    proj.mkdir(parents=True, exist_ok=True)
    p = _safe_file(proj, filename)
    p.parent.mkdir(parents=True, exist_ok=True)
    p.write_text(content, encoding="utf-8")
    return f"Geschrieben: {p} ({len(content)} Zeichen)"


@mcp.tool()
def run_command(name: str, program: str, args: list[str], timeout_s: int = 3600) -> str:
    """Führt ElmerSolver, ElmerGrid oder gmsh mit eigenen Argumenten im Projektordner aus
    (für Sonderfälle und Fehlersuche). program: 'ElmerSolver' | 'ElmerGrid' | 'gmsh'."""
    if program not in ("ElmerSolver", "ElmerGrid", "gmsh"):
        raise ValueError("Nur ElmerSolver, ElmerGrid oder gmsh erlaubt.")
    proj = _proj(name)
    proj.mkdir(parents=True, exist_ok=True)
    r = _run([_require(program)] + list(args), proj, timeout=timeout_s)
    return f"Exit-Code {r.returncode}\n--- Ausgabe (Ende) ---\n" + _tail(r.stdout + "\n" + r.stderr, 60)


if __name__ == "__main__":
    mcp.run()
