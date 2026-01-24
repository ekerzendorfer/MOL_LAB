# ==========================================
# ⚛️ MOL_LAB — Moleküllabor (v2 + Enhancements)
# ==========================================

import subprocess, sys, io, os, base64, math, time, shutil, stat

APP_VERSION = "MOL_LAB v2.1 (2026-01-24)"
APP_CREDIT  = "Erstellt von Erich Kerzendorfer mit ChatGPT"

# ---- Minimaler Import/Install-Helper ----
def _ensure_import(pkg, pip_name=None):
    try:
        __import__(pkg)
        return True
    except ImportError:
        name = pip_name or pkg
        subprocess.check_call([sys.executable, "-m", "pip", "install", "-q", name])
        return True

_ensure_import("rdkit", "rdkit")
_ensure_import("py3Dmol", "py3Dmol")
_ensure_import("ipywidgets", "ipywidgets")

from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors, Crippen
from rdkit.Chem import rdDepictor
from rdkit.Chem.Draw import rdMolDraw2D
from rdkit.Chem import rdMolTransforms

import py3Dmol
import ipywidgets as widgets
from ipywidgets import Layout
from IPython.display import display, SVG
from google.colab import files


class MoleculeApp:
    def __init__(self):
        self.mol = None
        self.mol_2d = None
        self.last_2d_svg = None
        self._viewer = None

        # Snapshots
        self.mol_mmff = None
        self.mol_xtb  = None

        # Cache für Overlay & Torsion
        self._cmp_cache = None     # (m1, m2, max_d) (ohne H, m2 aligned)
        self._max_torsion = None   # (i,j,k,l)

        # Beispiele
        self.examples = {
            "Koffein (Default)": {
                "name": "Caffeine",
                "smiles": "Cn1cnc2n(C)c(=O)n(C)c(=O)c12"
            },
            "Orforglipron": {
                "name": "Orforglipron",
                "smiles": "Fc1c(C)cc(N2N=C3CCN(C(=O)C=4N(c5ccc(N6CCOC(C)(C)C6)cc5C4)[C@@H](CCC)C7=NOC(O)=N7)[C@@H](C)C3=C2N8C=CN(C8=O)c9ccc%10N(N=Cc%10c9F)C)cc1C"
            },
            "Chlorophyll a": {
                "name": "Chlorophyll_a",
                "smiles": "CCC1=C(C2=NC1=CC3=C(C4=C([N-]3)C(=C5[C@H]([C@@H](C(=N5)C=C6C(=C(C(=C2)[N-]6)C=C)C)C)CCC(=O)OC/C=C(\\C)/CCCC(C)CCCC(C)CCCC(C)C)[C@H](C4=O)C(=O)OC)C)C.[Mg+2]"
            }
        }

        # Darstellung
        self.SPHERE_SCALE = 0.25
        self.STICK_RADIUS = 0.14
        self.OVERLAY_STICK_RADIUS = 0.18
        self.SURFACE_OPACITY = 0.55

        # Zoom: moderat, auf Achse zentriert
        self.ZOOM_FACTOR_AXIS = 1.85

        # Konformer / Optimierung
        self.NUM_CONFS = 12
        self.MINIMIZE_QUICK = 250
        self.MINIMIZE_FINAL = 1800

        # ---------------- UI ----------------
        style = {'description_width': 'initial'}

        self.info_box = widgets.HTML()

        self.dd_example = widgets.Dropdown(
            options=list(self.examples.keys()),
            value="Koffein (Default)",
            description="Beispiel:",
            layout=Layout(width="98%"),
            style=style
        )

        ex = self.examples[self.dd_example.value]
        self.input_smiles = widgets.Textarea(
            value=ex["smiles"],
            description="SMILES:",
            layout=Layout(width="98%", height="90px"),
            style=style
        )
        self.input_name = widgets.Text(
            value=ex["name"],
            description="Name:",
            layout=Layout(width="98%"),
            style=style
        )

        self.chk_show_h = widgets.Checkbox(
            value=True,
            description="H anzeigen (3D)",
            indent=False,
            layout=Layout(width="98%")
        )

        # 3D Darstellungsmodus (nicht für Overlay)
        self.view_mode = widgets.ToggleButtons(
            options=[("Ball & Stick", "ballstick"), ("Spacefill", "spacefill")],
            value="ballstick",
            description="3D:",
            style={'description_width': 'initial'},
            layout=Layout(width="98%")
        )

        def _on_view_mode_change(change):
            if self.mol is not None:
                # nur normale Ansicht neu rendern
                self.show_mol()
        self.view_mode.observe(_on_view_mode_change, names="value")

        # Buttons (schmäler, für Grid)
        b = dict(layout=Layout(width="100%"))
        self.btn_load = widgets.Button(description="1 2D", icon="search", **b)
        self.btn_3d   = widgets.Button(description="2 3D", icon="cube", **b)
        self.btn_opt  = widgets.Button(description="3 MMFF/UFF", icon="magic", **b)
        self.btn_ana  = widgets.Button(description="4 Surface", icon="atom", **b)

        self.btn_xtb  = widgets.Button(description="Adv xTB", icon="flask", **b)
        self.btn_cmp  = widgets.Button(description="Vergleich", icon="balance-scale", **b)
        self.btn_tor  = widgets.Button(description="Torsionen", icon="repeat", **b)
        self.btn_zoom = widgets.Button(description="Zoom Achse", icon="search-plus", **b)
        self.btn_axis = widgets.Button(description="Achse", icon="compress-arrows-alt", **b)

        self.btn_sdf  = widgets.Button(description="SDF", icon="download", **b)
        self.btn_pdb  = widgets.Button(description="PDB", icon="download", **b)
        self.btn_html = widgets.Button(description="3D-HTML", icon="file-code", **b)

        # Outputs
        self.out_2d = widgets.Output(layout={'border': '1px solid #999', 'height': '250px', 'padding': '4px'})
        self.out_3d = widgets.Output(layout={'border': '1px solid #999', 'height': '950px', 'padding': '4px'})
        self.out_log = widgets.Output(layout={'border': '1px solid #ddd', 'height': '170px', 'overflow_y': 'scroll', 'padding': '4px'})

        # Events
        self.dd_example.observe(self.on_example_change, names="value")
        self.btn_load.on_click(self.load_2d)
        self.btn_3d.on_click(self.make_3d)
        self.btn_opt.on_click(self.optimize_mmff)
        self.btn_ana.on_click(self.analyze_surface)

        self.btn_xtb.on_click(self.optimize_xtb_micromamba)
        self.btn_cmp.on_click(self.compare_mmff_xtb)
        self.btn_tor.on_click(self.compare_torsions)
        self.btn_zoom.on_click(self.zoom_to_max_torsion_axis)
        self.btn_axis.on_click(self.mark_torsion_axis)

        self.btn_sdf.on_click(self.save_sdf)
        self.btn_pdb.on_click(self.save_pdb)
        self.btn_html.on_click(self.save_3d_html)

        # Info-Box initial befüllen
        self._refresh_info_box()
        self.log(f"ℹ️ {APP_VERSION}")

        # --- START-Parameter aus Launcher übernehmen (falls gesetzt) ---
        try:
            start_smiles = globals().get("MOLLAB_START_SMILES", None)
            start_name   = globals().get("MOLLAB_START_NAME", None)
            auto_load_2d = globals().get("MOLLAB_AUTO_LOAD_2D", False)

            if start_smiles:
                self.input_smiles.value = str(start_smiles).strip()
            if start_name:
                self.input_name.value = str(start_name).strip()

            if auto_load_2d and start_smiles:
                self.load_2d(None)
        except Exception:
            pass

    # ---------- Logging ----------
    def log(self, txt):
        with self.out_log:
            print(txt)

    def clear_3d(self):
        self.out_3d.clear_output()
        self._viewer = None

    # ---------- Status / Info Box ----------
    def _rdkit_status(self):
        try:
            import rdkit
            path = getattr(rdkit, "__file__", "") or ""
            src = "pip"
            if "/content/rdkit_prefix" in path or "/conda" in path or "/mamba" in path:
                src = "conda-forge (micromamba)"
            ver = getattr(rdkit, "__version__", "unknown")
            return f"RDKit: {ver} • {src}"
        except Exception:
            return "RDKit: ❌ nicht geladen"

    def _xtb_status(self):
        if os.path.exists("/content/xtb_prefix/bin/xtb"):
            return "xTB (GFN2): ✅ verfügbar"
        return "xTB (GFN2): ⏳ nicht installiert (erst nach Adv xTB)"

    def _refresh_info_box(self):
        rd = self._rdkit_status()
        xt = self._xtb_status()
        self.info_box.value = f"""
        <div style="
            border:1px solid #ddd; border-radius:10px; padding:10px;
            background:#fafafa; font-size:12px; line-height:1.35;">
          <div style="font-weight:700;">{APP_VERSION}</div>
          <div style="color:#555;">{APP_CREDIT}</div>
          <div style="margin-top:6px;">{rd}</div>
          <div>{xt}</div>
        </div>
        """

    # ---------- Example ----------
    def on_example_change(self, change):
        ex = self.examples[change["new"]]
        self.input_smiles.value = ex["smiles"]
        self.input_name.value   = ex["name"]

    # ---------- 2D ----------
    def draw_2d_svg(self, mol):
        rdDepictor.SetPreferCoordGen(True)
        rdDepictor.Compute2DCoords(mol)
        drawer = rdMolDraw2D.MolDraw2DSVG(520, 240)
        opts = drawer.drawOptions()
        opts.bondLineWidth = 2
        opts.minFontSize = 14
        rdMolDraw2D.PrepareAndDrawMolecule(drawer, mol)
        drawer.FinishDrawing()
        return drawer.GetDrawingText()

    def load_2d(self, b):
        # Log leeren ist für Unterricht oft angenehm
        self.out_log.clear_output()

        # 3D-Fenster leeren, damit klar ist: wir starten neu
        self.clear_3d()

        smi = self.input_smiles.value.strip()
        mol = Chem.MolFromSmiles(smi)
        if not mol:
            return self.log("❌ Ungültiger SMILES.")

        self.mol_2d = Chem.Mol(mol)
        self.mol = None
        self.mol_mmff = None
        self.mol_xtb  = None
        self._cmp_cache = None
        self._max_torsion = None

        self.out_2d.clear_output()
        self.last_2d_svg = self.draw_2d_svg(self.mol_2d)
        with self.out_2d:
            display(SVG(self.last_2d_svg))

        self.log("✅ 2D-Struktur geladen.")
        self._refresh_info_box()

    # ---------- 3D ----------
    def make_3d(self, b):
        if self.mol_2d is None:
            return self.log("⚠️ Erst 2D laden.")

        self.clear_3d()
        base = Chem.AddHs(self.mol_2d)

        params = AllChem.ETKDGv3()
        params.randomSeed = 0xF00D

        ids = list(AllChem.EmbedMultipleConfs(base, numConfs=self.NUM_CONFS, params=params))
        if not ids:
            return self.log("❌ 3D-Embedding fehlgeschlagen.")

        bestE, bestId = None, ids[0]
        for cid in ids:
            try:
                ff = AllChem.UFFGetMoleculeForceField(base, confId=cid)
                ff.Minimize(maxIts=self.MINIMIZE_QUICK)
                E = ff.CalcEnergy()
                if bestE is None or E < bestE:
                    bestE, bestId = E, cid
            except:
                pass

        mol = Chem.Mol(base)
        mol.RemoveAllConformers()
        mol.AddConformer(base.GetConformer(bestId), assignId=True)
        self.mol = mol

        self.log("✅ 3D-Struktur erzeugt.")
        self.show_mol()

    # ---------- MMFF/UFF ----------
    def optimize_mmff(self, b):
        if self.mol is None:
            return self.log("⚠️ Erst 3D erzeugen.")

        mol = Chem.Mol(self.mol)
        used = "UFF"

        try:
            props = AllChem.MMFFGetMoleculeProperties(mol)
            if props:
                ff = AllChem.MMFFGetMoleculeForceField(mol, props)
                ff.Minimize(maxIts=self.MINIMIZE_FINAL)
                used = "MMFF94"
            else:
                ff = AllChem.UFFGetMoleculeForceField(mol)
                ff.Minimize(maxIts=self.MINIMIZE_FINAL)
        except:
            ff = AllChem.UFFGetMoleculeForceField(mol)
            ff.Minimize(maxIts=self.MINIMIZE_FINAL)

        self.mol = mol
        self.mol_mmff = Chem.Mol(self.mol)

        self.log(f"✅ Optimiert ({used}).")
        self.show_mol()

    # ---------- Viewer ----------
    def _mol_to_sdf(self, mol):
        buf = io.StringIO()
        w = Chem.SDWriter(buf)
        w.write(mol)
        w.close()
        return buf.getvalue()

    def _apply_style_normal(self, v):
        """Style für normale 3D-Ansichten (nicht Overlay)"""
        mode = getattr(self, "view_mode", None)
        mode = mode.value if mode else "ballstick"

        if mode == "spacefill":
            # Spacefill: Kugeln im VDW-Stil (ohne Sticks)
            v.setStyle({'sphere': {'scale': 1.0}})
        else:
            # Ball & Stick
            v.setStyle({
                'stick':  {'radius': self.STICK_RADIUS},
                'sphere': {'scale': self.SPHERE_SCALE}
            })

    def show_mol(self):
        if self.mol is None:
            return

        self.clear_3d()
        v = py3Dmol.view(width="100%", height=950)

        m = self.mol if self.chk_show_h.value else Chem.RemoveHs(self.mol)
        v.addModel(self._mol_to_sdf(m), "sdf")

        self._apply_style_normal(v)

        v.zoomTo()
        v.setBackgroundColor("white")

        self._viewer = v
        with self.out_3d:
            v.show()

    # ---------- Surface ----------
    def analyze_surface(self, b):
        if self.mol is None:
            return self.log("⚠️ Kein Molekül.")

        # Wichtig: an derselben Stelle rendern
        self.clear_3d()

        try:
            AllChem.ComputeGasteigerCharges(self.mol)
        except:
            pass

        v = py3Dmol.view(width="100%", height=950)
        v.addModel(self._mol_to_sdf(self.mol), "sdf")

        # Darstellung je nach Modus
        self._apply_style_normal(v)

        v.addSurface(py3Dmol.VDW, {
            'opacity': self.SURFACE_OPACITY,
            'colorscheme': {'prop': 'partialCharge', 'gradient': 'rwb'}
        })

        v.zoomTo()
        v.setBackgroundColor("white")

        self._viewer = v
        with self.out_3d:
            v.show()

        self.log("✅ Surface berechnet.")

    # ---------- Micromamba ----------
    def _ensure_micromamba(self):
        target = "/content/bin/micromamba"
        os.makedirs("/content/bin", exist_ok=True)

        if os.path.isfile(target):
            os.chmod(target, os.stat(target).st_mode | stat.S_IXUSR | stat.S_IXGRP | stat.S_IXOTH)
            return target

        cmd = "curl -Ls https://micro.mamba.pm/api/micromamba/linux-64/latest | tar -xvj -C /content"
        r = subprocess.run(cmd, shell=True, capture_output=True, text=True)

        candidates = [
            "/content/bin/micromamba",
            "/content/bin/bin/micromamba",
            "/content/micromamba",
        ]

        found = None
        for c in candidates:
            if os.path.isfile(c):
                found = c
                break

        if not found:
            raise RuntimeError("Micromamba Installation fehlgeschlagen:\n" + (r.stderr[:1200] or r.stdout[:1200]))

        if found != target:
            shutil.move(found, target)

        os.chmod(target, os.stat(target).st_mode | stat.S_IXUSR | stat.S_IXGRP | stat.S_IXOTH)
        return target

    def _ensure_xtb_prefix(self):
        mm = self._ensure_micromamba()
        prefix = "/content/xtb_prefix"
        if os.path.exists(os.path.join(prefix, "bin", "xtb")):
            return prefix

        self.log("🧪 Installiere xTB (conda-forge) ...")
        subprocess.check_call([mm, "create", "-y", "-p", prefix, "-c", "conda-forge", "xtb"],
                              stdout=subprocess.DEVNULL)
        return prefix

    # ---------- xTB ----------
    def _write_xyz(self, mol, path):
        conf = mol.GetConformer(0)
        lines = [str(mol.GetNumAtoms()), "mol"]
        for i, a in enumerate(mol.GetAtoms()):
            p = conf.GetAtomPosition(i)
            lines.append(f"{a.GetSymbol()} {p.x:.8f} {p.y:.8f} {p.z:.8f}")
        with open(path, "w") as f:
            f.write("\n".join(lines) + "\n")

    def _read_xyz(self, path):
        with open(path) as f:
            lines = f.readlines()[2:]
        out = []
        for ln in lines:
            p = ln.split()
            out.append((p[0], float(p[1]), float(p[2]), float(p[3])))
        return out

    def optimize_xtb_micromamba(self, b):
        if self.mol is None:
            return self.log("⚠️ Erst 3D erzeugen.")

        self.log("🧪 xTB Optimierung läuft ...")
        mol_copy = Chem.Mol(self.mol)

        try:
            prefix = self._ensure_xtb_prefix()
        except Exception as e:
            self._refresh_info_box()
            return self.log(f"❌ xTB Setup fehlgeschlagen: {e}")

        xtb = os.path.join(prefix, "bin", "xtb")
        workdir = "/content/xtb_work"
        os.makedirs(workdir, exist_ok=True)

        xyz = os.path.join(workdir, "input.xyz")
        self._write_xyz(mol_copy, xyz)

        cmd = [xtb, "input.xyz", "--gfn", "2", "--opt"]
        res = subprocess.run(cmd, cwd=workdir, capture_output=True, text=True, timeout=600)

        opt_xyz = os.path.join(workdir, "xtbopt.xyz")
        if res.returncode != 0 or not os.path.exists(opt_xyz):
            self._refresh_info_box()
            return self.log("❌ xTB Lauf fehlgeschlagen (keine xtbopt.xyz).")

        coords = self._read_xyz(opt_xyz)
        conf = mol_copy.GetConformer(0)
        for i, (_, x, y, z) in enumerate(coords):
            conf.SetAtomPosition(i, Chem.rdGeometry.Point3D(x, y, z))

        self.mol = mol_copy
        self.mol_xtb = Chem.Mol(self.mol)

        self.log("✅ xTB Optimierung übernommen.")
        self._refresh_info_box()
        self.show_mol()

    # ---------- Vergleich (Overlay bewusst als Stick) ----------
    def compare_mmff_xtb(self, b):
        if self.mol_mmff is None or self.mol_xtb is None:
            return self.log("⚠️ Erst MMFF/UFF und danach xTB rechnen.")

        m1 = Chem.RemoveHs(Chem.Mol(self.mol_mmff))
        m2 = Chem.RemoveHs(Chem.Mol(self.mol_xtb))

        try:
            rmsd = AllChem.AlignMol(m2, m1)
        except Exception as e:
            return self.log(f"❌ AlignMol fehlgeschlagen: {e}")

        c1 = m1.GetConformer(0)
        c2 = m2.GetConformer(0)

        dists = []
        for i in range(m1.GetNumAtoms()):
            p1 = c1.GetAtomPosition(i)
            p2 = c2.GetAtomPosition(i)
            d = math.dist((p1.x, p1.y, p1.z), (p2.x, p2.y, p2.z))
            dists.append(d)

        max_d = max(dists) if dists else 0.0
        self._cmp_cache = (Chem.Mol(m1), Chem.Mol(m2), float(max_d))

        self.log(f"📏 Overlay (ohne H, aligned): RMSD = {rmsd:.3f} Å | maxΔ = {max_d:.3f} Å")
        self._render_overlay(zoom_atoms=None, draw_axis=False, zoom_factor=None)

    # ---------- Torsionen ----------
    def compare_torsions(self, b, top_n=15):
        if self.mol_mmff is None or self.mol_xtb is None:
            return self.log("⚠️ Erst MMFF/UFF und danach xTB rechnen.")

        m1 = Chem.RemoveHs(Chem.Mol(self.mol_mmff))
        m2 = Chem.RemoveHs(Chem.Mol(self.mol_xtb))
        try:
            AllChem.AlignMol(m2, m1)
        except:
            pass

        torsions = []
        for bond in m1.GetBonds():
            if bond.GetBondType() != Chem.BondType.SINGLE or bond.IsInRing():
                continue

            a = bond.GetBeginAtomIdx()
            bidx = bond.GetEndAtomIdx()

            na = [n.GetIdx() for n in m1.GetAtomWithIdx(a).GetNeighbors() if n.GetIdx() != bidx]
            nb = [n.GetIdx() for n in m1.GetAtomWithIdx(bidx).GetNeighbors() if n.GetIdx() != a]
            if not na or not nb:
                continue

            i, j, k, l = na[0], a, bidx, nb[0]
            try:
                a1 = rdMolTransforms.GetDihedralDeg(m1.GetConformer(), i, j, k, l)
                a2 = rdMolTransforms.GetDihedralDeg(m2.GetConformer(), i, j, k, l)
                dd = abs(((a2 - a1 + 180) % 360) - 180)  # 0..180
                torsions.append((dd, (i, j, k, l), a1, a2))
            except:
                pass

        torsions.sort(reverse=True, key=lambda x: x[0])
        if not torsions:
            self._max_torsion = None
            return self.log("⚠️ Keine geeigneten Torsionen gefunden.")

        self._max_torsion = torsions[0][1]
        dd0, idxs0, a1_0, a2_0 = torsions[0]
        self.log(f"🧷 Größte Torsion: {idxs0} | Δ≈{dd0:.1f}° (MMFF {a1_0:.1f}° → xTB {a2_0:.1f}°)")
        self.log(f"Top {min(top_n, len(torsions))} Torsionen (|Δ| in °):")
        for dd, idxs, a1, a2 in torsions[:top_n]:
            self.log(f"  {dd:6.1f}°  {idxs}   ({a1:7.1f}° → {a2:7.1f}°)")

    # ---------- Overlay render helper ----------
    def _render_overlay(self, zoom_atoms=None, draw_axis=False, zoom_factor=1.85):
        if self._cmp_cache is None:
            return self.log("⚠️ Erst Vergleich ausführen.")

        m1, m2, _ = self._cmp_cache

        self.clear_3d()
        v = py3Dmol.view(width="100%", height=950)

        # Overlay bleibt Stick (didaktisch)
        v.addModel(self._mol_to_sdf(m1), "sdf")
        v.setStyle({'model': 0}, {'stick': {'radius': self.OVERLAY_STICK_RADIUS, 'colorscheme': 'cyanCarbon'}})

        v.addModel(self._mol_to_sdf(m2), "sdf")
        v.setStyle({'model': 1}, {'stick': {'radius': self.OVERLAY_STICK_RADIUS, 'colorscheme': 'magentaCarbon'}})

        if draw_axis and self._max_torsion is not None:
            _, j, k, _ = self._max_torsion
            c = m2.GetConformer(0)
            pj = c.GetAtomPosition(int(j))
            pk = c.GetAtomPosition(int(k))

            v.addCylinder({
                'start': {'x': pj.x, 'y': pj.y, 'z': pj.z},
                'end':   {'x': pk.x, 'y': pk.y, 'z': pk.z},
                'radius': 0.22,
                'fromCap': 1,
                'toCap': 1,
                'color': 'orange'
            })
            v.addSphere({'center': {'x': pj.x, 'y': pj.y, 'z': pj.z}, 'radius': 0.38, 'color': 'orange'})
            v.addSphere({'center': {'x': pk.x, 'y': pk.y, 'z': pk.z}, 'radius': 0.38, 'color': 'orange'})

        if zoom_atoms:
            try:
                v.zoomTo({'atom': [int(a) for a in zoom_atoms]})
                if zoom_factor:
                    v.zoom(float(zoom_factor))
            except:
                v.zoomTo()
                if zoom_factor:
                    v.zoom(float(zoom_factor))
        else:
            v.zoomTo()

        v.setBackgroundColor("white")
        self._viewer = v
        with self.out_3d:
            v.show()

    # ---------- Zoom & Achse (zentriert auf j-k) ----------
    def zoom_to_max_torsion_axis(self, b):
        if self._cmp_cache is None:
            return self.log("⚠️ Zoom: Bitte zuerst Vergleich ausführen.")
        if self._max_torsion is None:
            return self.log("⚠️ Zoom: Bitte zuerst Torsionsvergleich ausführen.")

        _, j, k, _ = self._max_torsion
        self._render_overlay(zoom_atoms=[j, k], draw_axis=False, zoom_factor=self.ZOOM_FACTOR_AXIS)
        self.log(f"🔎 Zoom auf Torsionsachse (j-k)=({j},{k})")

    def mark_torsion_axis(self, b):
        if self._cmp_cache is None:
            return self.log("⚠️ Achse: Bitte zuerst Vergleich ausführen.")
        if self._max_torsion is None:
            return self.log("⚠️ Achse: Bitte zuerst Torsionsvergleich ausführen.")

        _, j, k, _ = self._max_torsion
        self._render_overlay(zoom_atoms=[j, k], draw_axis=True, zoom_factor=self.ZOOM_FACTOR_AXIS)
        self.log(f"🟧 Torsionsachse markiert (j-k)=({j},{k})")

    # ---------- Export ----------
    def save_sdf(self, b):
        if self.mol is None:
            return self.log("⚠️ Kein Molekül zum Export.")
        fn = f"{self.input_name.value}.sdf"
        w = Chem.SDWriter(fn)
        w.write(self.mol)
        w.close()
        self.log(f"⬇️ SDF: {fn}")
        files.download(fn)

    def save_pdb(self, b):
        if self.mol is None:
            return self.log("⚠️ Kein Molekül zum Export.")
        fn = f"{self.input_name.value}.pdb"
        w = Chem.PDBWriter(fn)
        w.write(self.mol)
        w.close()
        self.log(f"⬇️ PDB: {fn}")
        files.download(fn)

    def save_3d_html(self, b):
        if self._viewer is None:
            return self.log("⚠️ Kein 3D-Viewer aktiv.")
        fn = f"{self.input_name.value}_3D.html"
        html = self._viewer._make_html()
        with open(fn, "w", encoding="utf-8") as f:
            f.write(html)
        self.log(f"💾 3D-HTML: {fn}")
        files.download(fn)

    # ---------- UI Render ----------
    def render(self):
        btn_grid = widgets.GridBox(
            children=[
                self.btn_load, self.btn_3d, self.btn_opt,
                self.btn_ana,  self.btn_xtb, self.btn_cmp,
                self.btn_tor,  self.btn_zoom, self.btn_axis,
                self.btn_sdf,  self.btn_pdb,  self.btn_html
            ],
            layout=Layout(
                width="98%",
                grid_template_columns="1fr 1fr 1fr",
                grid_gap="4px 6px"
            )
        )

        left = widgets.VBox([
            widgets.HTML("<h3>⚛️ Moleküllabor</h3>"),
            self.info_box,
            self.dd_example,
            self.input_smiles,
            self.input_name,
            self.chk_show_h,
            self.view_mode,
            widgets.HTML("<b>Aktionen</b>"),
            btn_grid,
            widgets.HTML("<b>2D-Struktur</b>"),
            self.out_2d,
            widgets.HTML("<b>Log</b>"),
            self.out_log
        ], layout=Layout(width="30%", padding="10px", border="1px solid #ccc"))

        right = widgets.VBox([
            widgets.HTML("<b>3D / Surface / Vergleich</b>"),
            self.out_3d
        ], layout=Layout(width="70%", padding="10px"))

        display(widgets.HBox([left, right]))


# START
app = MoleculeApp()
app.render()
