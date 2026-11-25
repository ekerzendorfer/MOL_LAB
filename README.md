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
