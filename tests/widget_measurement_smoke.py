"""Exercise real ipywidgets controls and measurement callbacks without Colab/WebGL."""
import ast
import contextlib
import io
import time
from pathlib import Path
from types import SimpleNamespace

import ipywidgets as widgets
from IPython.display import SVG
from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors, Crippen, rdMolAlign, rdMolTransforms, rdDepictor
from rdkit.Chem.Draw import rdMolDraw2D

root=Path(__file__).resolve().parents[1]
code=(root/"molekuellabor.py").read_text(encoding="utf-8")
tree=ast.parse(code)
members=[x for x in tree.body
         if isinstance(x,ast.ClassDef) and x.name=="MoleculeApp"
         or isinstance(x,ast.FunctionDef) and x.name=="calculate_geometry_measurement"
         or (isinstance(x,ast.Assign) and any(
             isinstance(target,ast.Name) and target.id=="WORKSHOP_EXAMPLES"
             for target in x.targets))]
assert len(members)==3

class Viewer:
    def __init__(self):
        self.spheres=[]
        self.labels=[]
        self.lines=[]
        self.models=[]
        self.shown=False
    def addModel(self,text,*args):
        self.models.append(text)
    def setStyle(self,*args):pass
    def addSphere(self,item):self.spheres.append(item)
    def addLabel(self,title,item):self.labels.append((title,item))
    def addLine(self,item):self.lines.append(item)
    def zoomTo(self):pass
    def setBackgroundColor(self,*args):pass
    def show(self):self.shown=True

viewers=[]
def fake_view(**kwargs):
    obj=Viewer()
    viewers.append(obj)
    return obj
shown=[]
ctx={
    "Chem":Chem,"AllChem":AllChem,"Descriptors":Descriptors,
    "Crippen":Crippen,"rdMolAlign":rdMolAlign,
    "rdMolTransforms":rdMolTransforms,"rdDepictor":rdDepictor,
    "rdMolDraw2D":rdMolDraw2D,
    "py3Dmol":SimpleNamespace(view=fake_view),
    "widgets":widgets,"Layout":widgets.Layout,
    "SVG":SVG,"display":lambda *obj: shown.extend(obj),
    "io":io,"time":time
}
exec(compile(ast.fix_missing_locations(ast.Module(body=members,type_ignores=[])),
             "mol_lab AST","exec"),ctx)
with contextlib.redirect_stdout(io.StringIO()):
    app=ctx["MoleculeApp"]()
    app.render()
workshop=ctx["WORKSHOP_EXAMPLES"]
assert len(workshop)==12
assert len({x["smiles"] for x in workshop})==12
assert all(Chem.MolFromSmiles(x["smiles"]) for x in workshop)
assert len(app.examples)==17,app.examples.keys()
example_smiles=[x["smiles"] for x in app.examples.values()]
assert len(example_smiles)==len(set(example_smiles))
assert {x["smiles"] for x in workshop}.issubset(set(example_smiles))
assert all(name in [x["name"] for x in app.examples.values()]
           for name in ("Ethanol","Dimethylether","Schwefelwasserstoff",
                        "Propanal","Aceton","Cyclohexan","Benzol","Koffein"))
# A preselected molecule from MW-024 stays an editable dropdown choice.
cyclo_key=next(k for k,v in app.examples.items()
               if v["smiles"]=="C1CCCCC1")
app.dd_example.value=cyclo_key
assert app.input_name.value=="Cyclohexan"
assert app.input_smiles.value=="C1CCCCC1"
assert app.examples["Konformer-Training · n-Hexan"]["smiles"]=="CCCCCC"
assert app.measure_accordion.selected_index is None
assert app.measure_panel in app.measure_accordion.children
assert app.measure_accordion in shown[-1].children[1].children, "must be in wide right column"
assert app.measure_mode.value=="distance"
assert app.measure_atom_boxes[0].layout.display=="block"
assert app.measure_atom_boxes[1].layout.display=="block"
assert app.measure_atom_boxes[2].layout.display=="none"
assert not app.measure_atom_boxes[0].options

# Loading the SMILES alone should never measure a 2D geometry.
app.input_name.value="Butan"
app.input_smiles.value="CCCC"
with contextlib.redirect_stdout(io.StringIO()):
    app.load_2d(None)
assert app.mol_2d is not None and app.mol is None
assert not app.measure_atom_boxes[0].options
app.perform_geometry_measurement(None)
assert "Messung nicht möglich" in app.measure_result.value

with contextlib.redirect_stdout(io.StringIO()):
    app.make_3d_robust(None)
assert app.mol and app.mol.GetNumConformers()==1
assert len(app.measure_atom_boxes[0].options)==app.mol.GetNumAtoms()
app.measure_mode.value="torsion"
assert all(d.layout.display=="block" for d in app.measure_atom_boxes)
for index,box in enumerate(app.measure_atom_boxes):
    box.value=index
with contextlib.redirect_stdout(io.StringIO()):
    app.perform_geometry_measurement(None)
assert "Torsionswinkel" in app.measure_result.value
assert "°" in app.measure_result.value
assert viewers[-1].shown and len(viewers[-1].spheres)==4
assert len(viewers[-1].labels)==4
assert viewers[-1].models
app.measure_mode.value="distance"
assert app.measure_atom_boxes[2].layout.display=="none"
with contextlib.redirect_stdout(io.StringIO()):
    app.perform_geometry_measurement(None)
assert "Abstand" in app.measure_result.value
assert "Å" in app.measure_result.value
assert len(viewers[-1].lines)==1

# Preserve current measurement source and expose conformer A/B as options.
app.conf_sampling.value=60
with contextlib.redirect_stdout(io.StringIO()):
    app.compute_conformers(None)
assert app.confs
assert "conformer_a" in [v for _,v in app.measure_source.options]
app.measure_source.value="conformer_a"
assert len(app.measure_atom_boxes[0].options)==app.confs[0]["mol"].GetNumAtoms()
with contextlib.redirect_stdout(io.StringIO()):
    app.perform_geometry_measurement(None)
assert "Abstand" in app.measure_result.value
assert len(viewers[-1].labels)==2

# Loading a new molecule must invalidate stale 3D and previous conformers.
app.input_smiles.value="O"
with contextlib.redirect_stdout(io.StringIO()):
    app.load_2d(None)
assert not app.confs
assert app.mol is None
assert [x[1] for x in app.measure_source.options]==["current"]
assert not app.measure_atom_boxes[0].options
app.clear_geometry_measurement()
assert "zurückgesetzt" in app.measure_result.value
print("ipywidgets: Colab-style geometry panel, modes, labels, sources and resets – OK")
