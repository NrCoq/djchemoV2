from rdkit import Chem
from rdkit.Chem.Draw import SimilarityMaps
from rdkit.DataStructs.cDataStructs import ConvertToNumpyArray
import hashlib
import numpy as np
from matplotlib import pyplot as plot

m = hashlib.md5()

m.update(b"CH3OH")

print(m.hexdigest())

m = hashlib.md5()

refsmi = "O=C(O)[C@@H]3[C@@H](O)C[C@@]2(O)C[C@@H](O)C[C@@H](O)[C@H](O)CC[C@@H](O)C[C@@H](O)CC(=O)O[C@@H]" \
      "(C)[C@H](C)[C@H](O)[C@@H](C)C=CC=CC=CC=CC=CC=CC=C[C@H](O[C@@H]1O[C@H](C)[C@@H](O)[C@H](N)[C@@H]1O)C[C@@H]3O2"

smi = "CC1C=CC=CCCC=CC=CC=CC=CC(CC2C(C(CC(O2)(CC(C(CCC(CC(CC(CC(=O)OC(C(C1O)C)C)O)O)O)O)O)O)O)C(=O)O)OC3C(C(C(C(O3)C)O)N)O"

refmol = Chem.MolFromSmiles(refsmi)
mol = Chem.MolFromSmiles(smi)

fp = SimilarityMaps.GetMorganFingerprint(mol, fpType='bv')
fig, maxweight = SimilarityMaps.GetSimilarityMapForFingerprint(refmol, mol, SimilarityMaps.GetMorganFingerprint)


fp = SimilarityMaps.GetMorganFingerprint(mol, fpType='bv')
arr = np.zeros((0,), dtype=np.int8)
ConvertToNumpyArray(fp, arr)
print(arr)