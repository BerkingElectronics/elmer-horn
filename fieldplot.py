"""Feldbilder aus Elmer-VTU-Dateien (ASCII) rendern – ohne VTK, nur mit matplotlib."""
from __future__ import annotations

import json
import math
import re
from pathlib import Path

P_REF = 20e-6


def read_vtu(path: Path):
    txt = Path(path).read_text(errors="replace")
    arrays = re.findall(r"<DataArray([^>]*)>(.*?)</DataArray>", txt, re.S)
    data = {}
    pts = None
    for attrs, body in arrays:
        m = re.search(r'Name="([^"]*)"', attrs)
        name = m.group(1).lower() if m else ""
        vals = body.split()
        if name == "connectivity":
            data["conn"] = [int(v) for v in vals]
        elif name == "offsets":
            data["offsets"] = [int(v) for v in vals]
        elif name == "types":
            data["types"] = [int(v) for v in vals]
        elif "pressure" in name and "p" not in data:
            data["p"] = [float(v) for v in vals]
        elif not name and pts is None:
            pts = [float(v) for v in vals]
    xs, ys = pts[0::3], pts[1::3]
    p = data["p"]
    ncomp = len(p) // len(xs)
    pc = [complex(p[i * ncomp], p[i * ncomp + 1]) for i in range(len(xs))]
    tris = []
    start = 0
    for off, typ in zip(data["offsets"], data["types"]):
        c = data["conn"][start:off]
        if typ == 5:                       # lineares Dreieck
            tris.append(c)
        elif typ == 22:                    # quadratisches Dreieck -> Eckknoten
            tris.append(c[:3])
        start = off
    return xs, ys, pc, tris


def plot_fields(vtu_files: list[tuple[float, Path]], project_json: Path, out: Path,
                title: str = "") -> Path:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import matplotlib.tri as mtri

    cfg = json.loads(Path(project_json).read_text())
    prof = cfg["profile_points"]            # [[r, z], ...]
    R = cfg["R"]
    n = len(vtu_files)
    fig, axes = plt.subplots(n, 2, figsize=(13, 4.6 * n), squeeze=False)
    for row, (f, path) in enumerate(vtu_files):
        xs, ys, pc, tris = read_vtu(path)
        # Spiegeln an der Achse: horizontal = Achse z, vertikal = ±r
        X = ys + ys
        Y = xs + [-x for x in xs]
        N = len(xs)
        T = tris + [[a + N, b + N, c + N] for a, b, c in tris]
        P = pc + pc
        tri = mtri.Triangulation(X, Y, T)
        spl = [20 * math.log10(max(abs(p) / math.sqrt(2), 1e-12) / P_REF) for p in P]
        smax = max(spl)
        re_ = [p.real for p in P]
        # Farbskala am Außenraum ausrichten (sonst dominiert der Hals)
        outside = sorted(abs(v) for v, x in zip(re_, X) if x > 0.05 * R)
        lim = (outside[int(0.98 * (len(outside) - 1))] if outside else max(abs(v) for v in re_)) or 1.0

        ax = axes[row][0]
        tc = ax.tripcolor(tri, spl, shading="gouraud", cmap="magma", vmin=smax - 40, vmax=smax)
        fig.colorbar(tc, ax=ax, label="SPL [dB]")
        ax.set_title(f"{f:.0f} Hz – Schalldruckpegel")
        ax2 = axes[row][1]
        tc2 = ax2.tripcolor(tri, re_, shading="gouraud", cmap="RdBu_r", vmin=-lim, vmax=lim)
        fig.colorbar(tc2, ax=ax2, label="Druck (Momentanwert) [Pa]")
        ax2.set_title(f"{f:.0f} Hz – Wellenbild")
        for a in (ax, ax2):
            for sgn in (1, -1):
                a.plot([z for r, z in prof], [sgn * r for r, z in prof], color="#e2e8f0", lw=1.6)
                a.plot([0, 0], [sgn * prof[-1][0], sgn * R], color="#e2e8f0", lw=1.6)
            a.set_aspect("equal")
            a.set_xlim(prof[0][1] - 0.02, R)
            a.set_ylim(-R, R)
            a.set_facecolor("#1a202c")
            a.set_xlabel("Achse z [m]")
            a.set_ylabel("Radius r [m]")
    if title:
        fig.suptitle(title)
    fig.tight_layout()
    fig.savefig(out, dpi=100)
    plt.close(fig)
    return out
