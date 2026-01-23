# ⚛️ MOL_LAB – Moleküllabor (RDKit + py3Dmol, Colab)

Interaktives Moleküllabor für Unterricht & Fortbildung:
- SMILES → 2D → 3D (Ball-and-Stick)
- Optimierung: MMFF/UFF + Advanced xTB (GFN2)
- Partial-Charge Surface (rwb)
- Vergleich: MMFF vs xTB (Overlay)
- Torsionsvergleich + Zoom/Markierung der größten Torsionsachse
- Export: SDF / PDB / 3D-HTML

## ▶️ Start (Open in Colab)
[![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/ekerzendorfer/MOL_LAB/blob/main/MOL_LAB.ipynb)

## Kurz-Anleitung (für Kolleg:innen)
1. Auf „Open in Colab“ klicken
2. In Colab: **Runtime → Run all** (oder die Startzelle ausführen)
3. Im Tool SMILES einfügen → **„1 2D“** → **„2 3D“** → optimieren / Surface / Vergleich

## Dateien im Repository
- **molekuellabor.py** – kompletter App-Code (UI + Logik)
- **MOL_LAB.ipynb** – Launcher-Notebook (lädt `molekuellabor.py` und startet die App)
- **README.md** – diese Anleitung

## Hinweis zu „Screenshots“
Direktes PNG-Screenshotting aus WebGL kann in Colab/Browsern manchmal leer sein.
→ Verwende **„3D-HTML“ Export**, öffne die HTML lokal und mache dann einen Screenshot (Win+Shift+S).
