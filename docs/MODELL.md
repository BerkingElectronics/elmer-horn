# Das Simulationsmodell im Detail

Diese Seite erklärt, was der Connector physikalisch rechnet, welche Annahmen dahinterstecken
und wie man die Ergebnisse einordnet. Sie ist bewusst so geschrieben, dass man sie auch
ohne FEM-Vorwissen nachvollziehen kann.

## 1. Die Grundidee: Helmholtz-Gleichung

Schall in Luft ist eine kleine Druckschwankung, die sich als Welle ausbreitet. Wenn man
nur **eine Frequenz zur Zeit** betrachtet (eingeschwungener Zustand, Sinuston), wird aus
der Wellengleichung die **Helmholtz-Gleichung**:

```
∇²p + k²·p = 0        mit  k = ω / c = 2πf / c
```

`p` ist hier der komplexe Schalldruck: Betrag = Amplitude, Winkel = Phase.
Die Finite-Elemente-Methode (FEM) zerlegt den Luftraum in viele kleine Dreiecke und löst
diese Gleichung näherungsweise auf jedem davon. Das übernimmt der `HelmholtzSolver` von Elmer.

Für einen Frequenzgang wird diese Rechnung einfach für viele Frequenzen wiederholt
(`run_sweep`). Da die Frequenzen unabhängig voneinander sind, laufen sie parallel.

Verwendete Konstanten (Luft bei 20 °C):

| Größe | Wert |
|---|---|
| Dichte ρ | 1,205 kg/m³ |
| Schallgeschwindigkeit c | 343 m/s |
| Bezugsschalldruck | 20 µPa |

## 2. Achsensymmetrie: 3D-Horn, 2D-Rechnung

Ein rundes Horn sieht in jedem Längsschnitt gleich aus. Elmer kann das mit dem
Koordinatensystem `"Axi Symmetric"` ausnutzen: Gerechnet wird nur ein **2D-Halbschnitt**
(Radius r nach außen, Achse z nach vorne), das Ergebnis gilt aber für das volle
rotationssymmetrische 3D-Horn. Das macht die Rechnung um Größenordnungen schneller als echtes 3D.

Konsequenz: Alles muss rotationssymmetrisch sein – Horn, Quelle und Schallfeld.

## 3. Die Geometrie

```
   r
   ↑
 R ┼  ·  ·  ·
   │          ·        Außenrand: Viertelkreis um den
   │ Schall-     ·     Mittelpunkt der Mündung, Radius R
   │ wand          ·
   │                 ·
 rₘ┼ ← Mündung        ·
  ╱│                   ·
 ╱ │   Luftraum         ·
╱  │                    ·
│Hals                   ·
┴──┴────────────────────┴──→ z   (Symmetrieachse)
−L  0                    R
```

Links von z = 0 liegt das Horn (Hornwand von r₀ am Hals bis rₘ an der Mündung),
rechts davon der halbkugelförmige Luftraum vor der Schallwand.

Die Randabschnitte bekommen in Gmsh feste Nummern („Physical Groups“), die in Elmer
direkt als Randbedingungs-Nummern weiterverwendet werden:

| Nr. | Rand | Randbedingung |
|---|---|---|
| 1 | Hals | Kolbenquelle (Schnelle vorgegeben) |
| 2 | Hornwand | schallhart |
| 3 | Schallwand | schallhart |
| 4 | Außenrand | absorbierend |
| 5 | Symmetrieachse | Symmetrie (ergibt sich automatisch) |

Der Hals liegt bei z = −L, die Mündung bei z = 0, bündig in der Schallwand.

### Profile

Alle Profile liefern den Radius r in Abhängigkeit vom Abstand z vom Hals
(r₀ = Halsradius, r_m = Mündungsradius, L = Länge, t = z/L):

| Profil | Formel |
|---|---|
| Konisch | r = r₀ + (r_m − r₀) · t |
| Exponentiell | r = r₀ · (r_m / r₀)^t |
| Hyperbolisch (Hypex, Salmon) | r = r₀ · (cosh(m·z) + T·sinh(m·z)), T = 0,7; m wird per Bisektion so bestimmt, dass bei z = L genau r_m herauskommt |
| Oblate Spheroidal (OS) | r = √(r₀² + (z · tan α)²), α = halber Abstrahlwinkel |
| Eigenes Profil | Stützpunkte aus `profile_points_mm`, in Gmsh als Spline verbunden |

### Größe des Luftraums

Standard: R = max(4 · r_m, 0,6 m, 2 · L). R muss deutlich größer als die Mündung sein,
damit der Außenrand im Fernfeld liegt (siehe Abschnitt 6).

## 4. Das Netz

FEM braucht genügend Elemente pro Wellenlänge, sonst wird die Welle „verschmiert“.
Faustregel: mindestens 6–10 Elemente pro Wellenlänge.

```
λ_min = c / f_max
h     = λ_min / elements_per_wavelength      (Standard: 10)
h_Hals = min(h, r₀ / 4)
```

Am Hals ist das Netz feiner (mindestens vier Elemente über den Halsradius), weil dort
die Impedanz ausgewertet wird. Entlang der Hornwand wird die Elementgröße linear vom
Hals- zum Normalwert übergeblendet.

Beispiel: f_max = 5 kHz → λ = 68,6 mm → h ≈ 6,9 mm.

Gmsh erzeugt das Netz (Algorithmus 6, „Frontal-Delaunay“), ElmerGrid wandelt es in
Elmers Format um.

## 5. Randbedingungen

### Hals: Kolbenquelle

Am Hals bewegt sich die Luft wie ein starrer Kolben mit der Schnelle **v₀ = 1 m/s**.
Über die Euler-Gleichung hängt die Schnelle mit dem Druckgradienten zusammen:

```
∂p/∂n = −iωρ · v_n      →      |∂p/∂n| = ω · ρ · v₀
```

In Elmer wird das als `Wave Flux` (Imaginärteil) gesetzt.

**Wichtig:** Weil die Schnelle fest vorgegeben ist und nicht die Kraft eines Treibers,
sind die Pegel „dB re 1 m/s Halsschnelle“. Ein echter Treiber liefert keine konstante
Schnelle: Seine Schnelle hängt von der Last ab, also von genau der Halsimpedanz, die
das Horn erzeugt. Die Pegelkurven zeigen deshalb das **Verhalten des Horns**, nicht den
Frequenzgang des fertigen Lautsprechers.

### Hornwand und Schallwand: schallhart

Keine Bewegung durch die Wand → ∂p/∂n = 0. Das ist in der FEM die „natürliche“
Randbedingung, die automatisch gilt, wenn man nichts setzt.

### Außenrand: absorbierend

Der echte Raum ist unendlich groß, das Netz aber nicht. Damit am Außenrand keine
Reflexionen zurücklaufen, wird dort eine **Plane-Wave-Randbedingung** gesetzt: Sie lässt
eine Welle, die senkrecht auftrifft, ungehindert hinaus. Weil der Außenrand ein Kreisbogen
um die Mündung ist, treffen die Wellen im Fernfeld fast senkrecht auf – dann funktioniert
das sehr gut. Alternativ: `absorber="impedance"` (Wellenimpedanz ρc).

## 6. Auswertung

Elmer schreibt mit `SaveLine` den Druck an allen Knoten von Hals und Außenrand in Textdateien.
Daraus berechnet der Connector:

### Halsimpedanz

Der Druck wird über die Halsfläche gemittelt (bei Achsensymmetrie mit dem Radius
gewichtet, weil äußere Ringe mehr Fläche haben) und durch ρ·c·v₀ geteilt:

```
z = p_mittel / (ρ · c · v₀)
```

Das ist die **normierte spezifische Strahlungsimpedanz** am Hals:

- **Re(z)** = Strahlungswiderstand → wie viel Leistung tatsächlich abgestrahlt wird.
  Re(z) ≈ 1 bedeutet: Das Horn verhält sich am Hals wie eine unendlich lange Röhre,
  der Treiber wird optimal belastet.
- **Im(z)** = Blindanteil → mitschwingende Luftmasse, die keine Leistung abstrahlt.
- Unter der **Grenzfrequenz** des Horns (und unterhalb von ka ≈ 1 an der Mündung) bricht
  Re(z) ein, und Welligkeit zeigt Reflexionen an der Mündung.

Falls Re(z) negativ herauskommt, wird das Vorzeichen umgedreht (Vorzeichenkonvention der
Zeitabhängigkeit in Elmer). `get_results` weist dann darauf hin.

### Pegel unter Winkeln

Auf dem Außenrand (Radius R) wird der Druck unter 0°, 15°, …, 90° interpoliert
(0° = Achse, 90° = entlang der Schallwand) und in einen Pegel umgerechnet:

```
SPL(1 m) = 20 · log10( |p| / √2 / 20 µPa ) + 20 · log10( R / 1 m )
```

Der zweite Term rechnet vom Radius R auf 1 m um und nimmt dafür an, dass der Druck mit
1/r abfällt. Das stimmt nur im **Fernfeld**. Bei großen Mündungen und hohen Frequenzen
liegt das Fernfeld weiter draußen (grob ab D²/λ, D = Mündungsdurchmesser) – dann sind
besonders die Werte außerhalb der Achse nur Näherungen. Abhilfe: `domain_radius_m`
größer wählen (kostet Rechenzeit).

### Frequenzraster

Die Frequenzen sind logarithmisch verteilt, mit `points_per_octave` Punkten pro Oktave
(Standard 6): f_i = f_start · 2^(i / ppo).

## 7. Was das Modell nicht kann

| Nicht abgebildet | Auswirkung |
|---|---|
| Endliche Schallwand / Gehäuse | Keine Kantenbeugung; im echten Gehäuse sieht die Richtwirkung unten herum anders aus |
| Nicht-runde Hörner | Rechteckige Hörner, unterschiedliche Winkel horizontal/vertikal: braucht 3D |
| Treiber | Kein Kompressionstreiber, keine Kompressionskammer, kein Phase Plug, keine TSP |
| Verluste | Luft ist verlustfrei; bei langen, engen Hörnern (Bass) werden Pegel überschätzt |
| Nichtlinearität | Bei sehr hohen Pegeln am Hals entstehen Verzerrungen, die die lineare Helmholtz-Gleichung nicht kennt |

Für den Vergleich von Hornprofilen, die Wahl von Länge und Mündung und ein Gefühl für
Richtwirkung und Ladung ist das Modell trotzdem sehr aussagekräftig.

## 8. Weiterlesen

- Elmer Models Manual, Kapitel „Helmholtz Equation“: <https://www.elmerfem.org/>
- Gmsh-Dokumentation: <https://gmsh.info/doc/texinfo/gmsh.html>
- Earl Geddes: *Audio Transducers* (zu OS-Waveguides und Hornladung)
- Vincent Salmon (1946): Generalized Plane Wave Horn Theory (Hypex-Hörner)
