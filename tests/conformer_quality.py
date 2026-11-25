"""RDKit-level regression of MOL_LAB minimum search; Colab UI is intentionally not imported."""
import ast
import json
import math
import time
from pathlib import Path
from types import SimpleNamespace

from rdkit import Chem
from rdkit.Chem import AllChem, rdMolAlign, rdMolTransforms

root=Path(__file__).resolve().parents[1]
source=(root/"molekuellabor.py").read_text(encoding="utf-8")
tree=ast.parse(source)
klass=next(x for x in tree.body if isinstance(x,ast.ClassDef) and x.name=="MoleculeApp")
method=next(x for x in klass.body if isinstance(x,(ast.FunctionDef,ast.AsyncFunctionDef))
            and x.name=="compute_conformers")
module=ast.fix_missing_locations(ast.Module(body=[method],type_ignores=[]))
ctx={"Chem":Chem,"AllChem":AllChem,"rdMolAlign":rdMolAlign,"rdMolTransforms":rdMolTransforms,"time":time}
exec(compile(module,str(root/"molekuellabor.py"),"exec"),ctx)
compute=ctx["compute_conformers"]

class Dummy:
    def __init__(self,smiles):
        self.mol_2d=Chem.MolFromSmiles(smiles)
        self.dd_conf_a=SimpleNamespace(options=[],value=None)
        self.dd_conf_b=SimpleNamespace(options=[],value=None)
        self.conf_view_mode=SimpleNamespace(value=None)
        self.msg=[]
        self.confs=[]
        self.view_count=0
        self.conf_sampling=SimpleNamespace(value=120)
        self.conf_energy_window=SimpleNamespace(value=15.0)
        self.conf_selection=SimpleNamespace(value="spread")
        self.conf_summary=SimpleNamespace(value="")
    def _conformer_eligibility(self,m):
        return True,[],{"heavy_atoms":m.GetNumHeavyAtoms()}
    def log(self,msg):
        self.msg.append(msg)
    def show_conformer_view(self,_):
        self.view_count+=1
    def _refresh_measure_sources(self):
        pass  # UI-only callback, unrelated to the pure RDKit test

def run(smiles):
    dummy=Dummy(smiles)
    compute(dummy,None)
    assert 2 <= len(dummy.confs) <= 10,(smiles,dummy.msg)
    assert dummy.view_count==1
    assert dummy.dd_conf_a.value==0 and dummy.dd_conf_b.value==1
    assert any("verschiedene Minima erkannt" in m for m in dummy.msg)
    assert all(d["dE"]>=0 for d in dummy.confs)
    assert all(math.isfinite(d["E"]) for d in dummy.confs)
    assert all(d["E"]==d["E"] for d in dummy.confs)
    assert all(dummy.confs[i]["E"]<=dummy.confs[i+1]["E"]+1e-7
               for i in range(len(dummy.confs)-1))
    mols=[Chem.RemoveHs(c["mol"]) for c in dummy.confs]
    for i in range(len(mols)):
        for j in range(i+1,len(mols)):
            rms=rdMolAlign.GetBestRMS(Chem.Mol(mols[i]),mols[j],maxMatches=128)
            assert rms >= .28, (smiles,i,j,rms)
    print(smiles,len(dummy.confs),"unique minima, relative E:",
          [round(x["dE"],2) for x in dummy.confs])
    return dummy

c=run("C1CCCCC1")
# Chair and twist-boat minima occur; the boat/half-chair remain absent.
assert c.confs[0]["dE"]==0.0
assert c.confs[0]["ring_label"]=="Sessel"
assert any(x["ring_label"].startswith("Twist-Boat") for x in c.confs)
assert sum("Twist-Boat · Variante" in x["ring_label"] for x in c.confs)==2
assert any("Twist-Boat" in str(option[0]) for option in c.dd_conf_a.options)

assert max(x["dE"] for x in c.confs)>3.0
assert "Warum nicht mehr?" in c.conf_summary.value
assert "zwei gleichenergetische" in c.conf_summary.value
b=run("CCCC")
assert sorted(x["ring_label"] for x in b.confs)==[
    "anti","gauche (+)","gauche (−)"
]
assert "zwei Energiegruppen" in b.conf_summary.value
run("ClCCCl")


# Higher-energy *minima* may be retained with a wider window and diverse mode.
test=Dummy("CCCCCC")
test.conf_energy_window.value=35.0
test.conf_sampling.value=240
test.conf_selection.value="spread"
compute(test,None)
assert test.confs and len(test.confs)<=10
assert "Höhere Energien" in test.conf_summary.value
assert all(c["dE"] <= 35.01 for c in test.confs)
assert all(test.confs[i]["E"]<=test.confs[i+1]["E"] for i in range(len(test.confs)-1))
print("Hexan broad min search:",len(test.confs),
      [round(c["dE"],2) for c in test.confs])

# The selection modes MUST differ when enough distinct minima exist:
# energy mode retains the ten lowest; spread mode includes a higher one.
lowest=Dummy("CCCCCC")
lowest.conf_energy_window.value=35.0
lowest.conf_sampling.value=240
lowest.conf_selection.value="energy"
compute(lowest,None)
assert len(lowest.confs)==10
assert max(c["dE"] for c in test.confs) > max(
    c["dE"] for c in lowest.confs
) + 0.5,(test.confs[-1]["dE"],lowest.confs[-1]["dE"])
assert test.confs[0]["dE"]==lowest.confs[0]["dE"]==0.0
print("Hexan: energy vs spread selects different ranges:",
      round(lowest.confs[-1]["dE"],2),
      round(test.confs[-1]["dE"],2))

notebook=json.loads((root/"MOL_LAB.ipynb").read_text(encoding="utf-8"))
code="".join(notebook["cells"][1]["source"])
assert 'try_download("main")' in code
assert 'try_download("fix-konformere-postopt-rms")' not in code
assert "params.pruneRmsThresh = -1.0" in source
assert "self.conf_sampling = widgets.Dropdown" in source
assert "self.conf_energy_window = widgets.FloatSlider" in source
assert "self.conf_selection = widgets.Dropdown" in source
assert "GetBestRMS" in source
assert "if ia == ib:" in source
assert notebook["cells"][1]["metadata"]["cellView"]=="form"

assert "RMSD nach Ausrichtung" in source
print("MOL_LAB: post-minimization RMSD, non-duplicate energy candidates and branch loading OK")

# MW-024 launcher: form only, preset SMILES/name, shared code loader on test branch.
mw=json.loads((root/"MOL_LAB_MW024.ipynb").read_text(encoding="utf-8"))
assert mw["cells"][1]["metadata"]["cellView"]=="form"
launch="".join(mw["cells"][1]["source"])
assert 'START_SMILES = "C1CCCCC1"' in launch
assert 'START_NAME = "Cyclohexan"' in launch
assert 'try_download("main")' in launch
assert 'try_download("fix-konformere-postopt-rms")' not in launch
assert 'display-mode: "form"' in launch
