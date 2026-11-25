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
    def _conformer_eligibility(self,m):
        return True,[],{"heavy_atoms":m.GetNumHeavyAtoms()}
    def log(self,msg):
        self.msg.append(msg)
    def show_conformer_view(self,_):
        self.view_count+=1

def run(smiles):
    dummy=Dummy(smiles)
    compute(dummy,None)
    assert 2 <= len(dummy.confs) <= 10,(smiles,dummy.msg)
    assert dummy.view_count==1
    assert dummy.dd_conf_a.value==0 and dummy.dd_conf_b.value==1
    assert any("geometrisch verschiedene Minima" in m for m in dummy.msg)
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
assert any(x["ring_label"]=="Twist-Boat" for x in c.confs)
assert any("Twist-Boat" in str(option[0]) for option in c.dd_conf_a.options)

assert max(x["dE"] for x in c.confs)>3.0
run("CCCC")
run("ClCCCl")

notebook=json.loads((root/"MOL_LAB.ipynb").read_text(encoding="utf-8"))
code="".join(notebook["cells"][1]["source"])
assert 'try_download("fix-konformere-postopt-rms")' in code
assert "params.pruneRmsThresh = -1.0" in source
assert "GetBestRMS" in source
assert "if ia == ib:" in source
assert "RMSD nach Ausrichtung" in source
print("MOL_LAB: post-minimization RMSD, non-duplicate energy candidates and branch loading OK")
