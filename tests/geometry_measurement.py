"""Independent RDKit geometries, including four-atom torsions, with UI-source checks."""
import ast
import math
from pathlib import Path

from rdkit import Chem
from rdkit.Chem import AllChem, rdMolTransforms

root=Path(__file__).resolve().parents[1]
s=(root/"molekuellabor.py").read_text(encoding="utf-8")
tree=ast.parse(s)
fn=next(x for x in tree.body if isinstance(x,ast.FunctionDef)
        and x.name=="calculate_geometry_measurement")
ctx={"rdMolTransforms":rdMolTransforms}
exec(compile(ast.fix_missing_locations(ast.Module(body=[fn],type_ignores=[])),
             "measure method","exec"),ctx)
measure=ctx["calculate_geometry_measurement"]

m=Chem.AddHs(Chem.MolFromSmiles("CCCC"))
assert AllChem.EmbedMolecule(m,randomSeed=173)==0
assert AllChem.MMFFOptimizeMolecule(m)==0
distance=measure(m,[0,1],"distance")
angle=measure(m,[0,1,2],"angle")
torsion=measure(m,[0,1,2,3],"torsion")
assert distance["unit"]=="Å" and 1.35<distance["value"]<1.7
assert angle["unit"]=="°" and 98<angle["value"]<122
assert torsion["unit"]=="°" and -180<=torsion["value"]<=180
assert abs(torsion["value"]-rdMolTransforms.GetDihedralDeg(
    m.GetConformer(),0,1,2,3))<1e-9
assert torsion["indices"]==(0,1,2,3)
for mode,ids in [("torsion",[0,1,2,2]),("angle",[0,0,1]),
                 ("torsion",[0,1,3,2]),("distance",[0,100])]:
    try:
        measure(m,ids,mode)
        raise AssertionError("Invalid "+mode+" ids accepted: "+str(ids))
    except ValueError: pass
no3d=Chem.MolFromSmiles("CCCC")
try:
    measure(no3d,[0,1],"distance")
    raise AssertionError("Unembedded geometry accepted")
except ValueError:pass

# Complex drug structures can be measured without conformer eligibility.
big=Chem.AddHs(Chem.MolFromSmiles(
 "CC(C)C1=C(C(=C(N1CC[C@H](C[C@H](CC(=O)O)O)O)C2=CC=C(C=C2)F)C3=CC=CC=C3)C(=O)NC4=CC=CC=C4"
))
assert big.GetNumHeavyAtoms()==41
assert AllChem.EmbedMolecule(big,randomSeed=173,maxAttempts=50)==0
assert measure(big,[0,1],"distance")["value"]>1.0

cls=next(x for x in tree.body if isinstance(x,ast.ClassDef)
         and x.name=="MoleculeApp")
methods={x.name for x in cls.body if isinstance(x,ast.FunctionDef)}
for name in ("_measurement_molecule","_refresh_measure_sources",
             "_refresh_measure_atoms","perform_geometry_measurement",
             "clear_geometry_measurement","_update_measure_visibility"):
    assert name in methods
assert "self.measure_accordion" in s
assert "self.measure_accordion," in s
assert "self.measure_source" in s
assert "v.addSphere" in s and "v.addLabel" in s
assert "self.confs = []" in s
print("Geometrie: Abstand/Winkel/Torsion, Fehlerfälle und 41-Schweratom-Wirkstoff – OK")
