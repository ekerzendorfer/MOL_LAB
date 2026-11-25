# ⚛️ MOL_LAB – Moleküllabor (für Unterricht & Fortbildung)

Dieses Notebook zeigt, wie aus einer **2D-Struktur (SMILES)** ein **3D-Molekül** wird – inkl. Optimierung, Surface und Vergleich.

## ▶ Start (nur 1 Klick)
[![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/ekerzendorfer/MOL_LAB/blob/main/MOL_LAB.ipynb)

## So funktioniert’s (für Schüler:innen)
1. Oben im Startbereich optional **SMILES** ändern (oder Standard lassen).
2. **▶ START MOL_LAB** einmal ausführen (Run-Button der Zelle).
3. Im Tool:
   - **1 2D** → zeigt die Strukturformel
   - **2 3D** → erzeugt ein räumliches Modell
   - **3 MMFF/UFF** → optimiert die Geometrie
   - **4 Surface** → zeigt Ladungsverteilung (rot/blau)
   - Optional: **Adv xTB**, **Vergleich**, **Torsionen**, **Achse/Zoom**

## Tipps
- **SMILES** bekommst du z.B. aus Ketcher/MolView oder aus PubChem.
- Wenn etwas nicht gleich läuft: **Runtime → Restart runtime** und nochmal START drücken.

## Dateien im Repo
- `MOL_LAB.ipynb` – Start-Notebook (ein Klick)
- `molekuellabor.py` – die App (UI + Logik)
- `README.md` – diese Anleitung

## Hinweis zu „Screenshots“
PNG-Screenshots aus WebGL können in Colab manchmal leer sein.
→ Nutze **3D-HTML Export**, öffne die HTML lokal und mache dann einen Screenshot.


### Konformersuche: geometrische Duplikate vermeiden (Entwicklungsbranch, 25.11.2026)

Die bisherige Einstellung `ETKDG.pruneRmsThresh = 0.75` verwarf möglicherweise unterschiedliche Startkonformationen vor der Energieminimierung. Anschließend wurden die energieärmsten zehn Kandidaten auch dann ausgegeben, wenn sie nach MMFF-Optimierung dieselbe Geometrie besaßen. Im Entwicklungsbranch `fix-konformere-postopt-rms` wird zuerst eine vollständige Stichprobe eingebettet (`pruneRmsThresh = -1`) und minimiert; erst danach werden maximal zehn **geometrisch verschiedene lokale Minima** mittels symmetrieberücksichtigendem Schweratom-RMSD (Grenze 0,30 Å) ausgewählt.

Bei Cyclohexan findet die normale Minimasuche insbesondere Sessel und ggf. Twist-Boat; **Wanne und Halbsessel** sind keine stabilen Minima und werden damit nicht automatisch gefunden. Für einen **gezielten Vergleich aller vier Ringformen** sollte DIDACTIC V0.8 verwendet werden. Zwei Ergebnisse mit gleicher Energie können räumlich unterschiedlich sein; identische Geometrien sollen aber nicht als verschiedene Treffer ausgegeben werden. Beim Overlay einer zweimal ausgewählten Konformer-Nummer wird kein irreführender Vergleich mehr angezeigt.

**Entwicklungs-Notebook (lädt bewusst den Code dieses Branches):** https://colab.research.google.com/github/ekerzendorfer/MOL_LAB/blob/fix-konformere-postopt-rms/MOL_LAB.ipynb . Das stabile Notebook auf `main` bleibt bis zur Abnahme unverändert.

## Entwicklungsstand V0.2 – Startformulare, Minima-Vielfalt und Geometriemessung (25.11.2026)

**Start:** `MOL_LAB.ipynb` zeigt den Start standardmäßig als Colab-Formular (`cellView=form`). `MOL_LAB_MW024.ipynb` enthält bereits `START_SMILES="C1CCCCC1"`, `START_NAME="Cyclohexan"`; das MW-024-Entwicklungsnotebook verlinkt diese vorbereitete Colab-Datei. Das Notebook lädt weiterhin **denselben** Code aus `molekuellabor.py`.

**Konformersuche:** steuerbare 60, 120 oder 240 ETKDG-Startgeometrien; nach MMFF/UFF-Minimierung RMSD-Deduplikation, optional bis +40 kcal/mol Energiefenster. Anzeige von maximal zehn unterschiedlichen, lokalen Minima: entweder energieaufsteigend oder über den Energiebereich verteilt. Das Energiefenster ist **keine** Suche nach Übergangszuständen. Für größere (>18 Schweratome) Moleküle maximal 120 Starts; ein zeitliches Laufzeitlimit ist weiterhin aktiv.

**Geometrie messen:** einklappbares Panel unabhängig von der Konformersuche, Atomwahl mittels nummerierter Dropdowns, Entfernung (2 Atome, Å), Bindungswinkel (3 gebundene Atome, °), Torsionswinkel (4 aufeinanderfolgende gebundene Atome, °). Ausgewählte Atome werden in eigenem 3D-Fenster beschriftet und farbig markiert; H kann bei Bedarf angezeigt werden. Wahl der Quelle: aktuelles 3D-Modell oder vorhandener Konformer A/B. Die Funktion ist nicht an die 30-Schweratom-Grenze der Konformersuche gebunden (etwa Atorvastatin ist geometrisch messbar, wenn ein 3D-Modell vorliegt).

**Bewusst verschoben:** separate manuelle Mausrotation der beiden Overlay-Modelle und direkte Atomwahl durch Anklicken der 3D-Grafik. Die Dropdown-Auswahl ist zunächst die robuste Colab-Version. Noch nicht auf `main` veröffentlicht; Browser-/Colab-Tests folgen nach automatischer Regression.

### Werkstatt-Moleküle und Konformerfamilien (25.11.2026)

Das MOL_LAB-Beispiel-Dropdown enthält jetzt alle **12 unterschiedlichen Moleküle** der MOLEKÜL-WERKSTATT: Ethanol, Dimethylether, Wasser, Schwefelwasserstoff, Propanal, Aceton, Kohlenstoffdioxid, Ammoniak, Methan, Cyclohexan, Benzol und Koffein. Grundlage: `learning_paths[].molecules` und `free_examples` aus `ekerzendorfer/MOLEKUELWERKSTATT/data/learning_paths.json` (Stand des Entwicklungszweigs `mw024-cyclohexan-v08`, 25.11.2026). Duplikate nach SMILES werden entfernt; das bereits eingebaute Koffein wird nur einmal angeboten. Die Werkstatt-Beispieldaten sind als offline-verwendbarer Snapshot in `molekuellabor.py` enthalten; **bei späteren neuen Forschungswegen ist dieser Snapshot bewusst zu aktualisieren**. Neue Beispiele in der Werkstatt erscheinen daher nicht automatisch im MOL_LAB-Dropdown.

Für umfangreichere Konformeruntersuchungen wurden **Butan, n-Pentan und n-Hexan** ergänzt. Bei Butan werden die drei lokalen Minima als anti, gauche (+) und gauche (−) bezeichnet; Cyclohexan zeigt Sessel und die beiden Twist-Boat-Varianten A/B. Damit sind trotz drei räumlicher Minima nur zwei unterschiedliche energetische Familien vertreten.

Der Auswahlmodus `Vielfalt über das Energiefenster` verteilt die **maximal zehn ausgewählten Minima nach ihrer Energie über den gefundenen Bereich**, falls **mehr als zehn** verschiedene Minima gefunden wurden. Bei höchstens zehn wird jedes gefundene Minimum bereits angezeigt; ein Moduswechsel oder eine große Energieobergrenze kann dann **keine zusätzlichen Minima erfinden**. Die Oberfläche gibt in diesem Fall einen eigenen Hinweis aus. Wanne und Halbsessel werden nur im geführten DIDACTIC-Vergleich als nicht-minimale Referenzgeometrien behandelt.

Beim Start mit einem vorbereiteten Colab-Notebook (insbesondere MW-024) wird nicht nur SMILES/Name, sondern auch das passende Beispiel im Dropdown ausgewählt.
