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

