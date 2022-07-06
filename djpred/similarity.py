from rdkit import Chem
from rdkit.Chem.Draw import SimilarityMaps
from rdkit.DataStructs.cDataStructs import ConvertToNumpyArray
import numpy as np


def calc_similarity(smiles):

	amphotericin_b = "O=C(O)[C@@H]3[C@@H](O)C[C@@]2(O)C[C@@H](O)C[C@@H](O)[C@H](O)CC[C@@H](O)C[C@@H](O)CC" \
		                 "(=O)O[C@@H](C)[C@H](C)[C@H](O)[C@@H](C)C=CC=CC=CC=CC=CC=CC=C[C@H](O[C@@H]1O[C@H](C)" \
		                 "[C@@H](O)[C@H](N)[C@@H]1O)C[C@@H]3O2"

	nystatin = "CC1C=CC=CCCC=CC=CC=CC=CC(CC2C(C(CC(O2)(CC(C(CCC(CC(CC(CC(=O)OC(C" \
		           "(C1O)C)C)O)O)O)O)O)O)O)C(=O)O)OC3C(C(C(C(O3)C)O)N)O"

	natamycin = "OC(=O)[C@@H]3[C@@H](O)C[C@@]2(O)C[C@@H](O)C[C@H]4O[C@@H]4/C=C/C(=O)O[C@H](C)C\C=C\C=C" \
		            "\C=C\C=C\[C@H](O[C@@H]1O[C@H](C)[C@@H](O)[C@H](N)[C@@H]1O)C[C@@H]3O2"

	filipin = "CCCCC[C@H](O)[C@H]1C(=O)O[C@H](C)[C@@H](O)\C=C\C=C\C=C\C=C\C=C(/C)[C@@H](O)C[C@H](O)C[C@H](O)" \
	          "C[C@H](O)C[C@H](O)C[C@H](O)C[C@@H]1O"

	hamycin = "CC1C=CC=CC=CC=CC=CC=CC=CC(CC2C(C(CC(O2)(CC(CC(CC(CC(CC(CC(CC(=O)OC1C(C)CCC(CC(=O)C3=CC=C(C=C3)N)" \
	          "O)O)O)O)O)O)O)O)O)C(=O)O)OC4C(C(C(C(O4)C)O)N)O"

	ampo_mol = Chem.MolFromSmiles(amphotericin_b, sanitize=True)
	nysta_mol = Chem.MolFromSmiles(nystatin, sanitize=True)
	natam_mol = Chem.MolFromSmiles(natamycin, sanitize=True)
	#filip_mol = Chem.MolFromSmiles(filipin, sanitize=True)
	#hamyc_mol = Chem.MolFromSmiles(hamycin, sanitize=True)
	query_mol = Chem.MolFromSmiles(smiles)

	ampo_fp = SimilarityMaps.GetRDKFingerprint(ampo_mol, fpType='bv', nBits=2048)
	nysta_fp = SimilarityMaps.GetRDKFingerprint(nysta_mol, fpType='bv', nBits=2048)
	natam_fp = SimilarityMaps.GetRDKFingerprint(natam_mol, fpType='bv', nBits=2048)
	#filip_fp = SimilarityMaps.GetRDKFingerprint(filip_mol, fpType='bv', nBits=2048)
	#hamyc_fp = SimilarityMaps.GetRDKFingerprint(hamyc_mol,fpType='bv', nBits=2048)
	query_fp = SimilarityMaps.GetRDKFingerprint(query_mol, fpType='bv', nBits=2048)

	bench_np = generate_benchmark_fp(ampo_fp, nysta_fp, natam_fp)

	return calculate_similarity(bench_np, query_fp)

def generate_benchmark_fp(ampo_fp, nysta_fp, natam_fp):

	ampo_np = np.zeros((0,), dtype=np.int8)
	nysta_np = np.zeros((0,), dtype=np.int8)
	natam_np = np.zeros((0,), dtype=np.int8)
	#filip_np = np.zeros((0,), dtype=np.int8)
	#hamyc_np = np.zeros((0,), dtype=np.int8)

	ConvertToNumpyArray(ampo_fp, ampo_np)
	ConvertToNumpyArray(nysta_fp, nysta_np)
	ConvertToNumpyArray(natam_fp, natam_np)
	#ConvertToNumpyArray(filip_fp, filip_np)
	#ConvertToNumpyArray(hamyc_fp, hamyc_np)


	for i in range(len(ampo_np)):
		if ampo_np[i]==nysta_np[i] & ampo_np[i]==natam_np[i]:
			continue
		else:
			# 3 serves as a 'this bit is not relevant to polyene function' marker
			ampo_np[i] = 0
	return ampo_np

def calculate_similarity(bench_np, query_fp):
	important_bits = 0
	matching_bits = 0

	query_np = np.zeros((0,), dtype=np.int8)
	ConvertToNumpyArray(query_fp, query_np)

	for i in range(len(bench_np)):
		if bench_np[i] == 1:
			important_bits += 1
			if query_np[i] == bench_np[i]:
				matching_bits += 1
	print(important_bits)
	print(matching_bits)
	return float(matching_bits)/float(important_bits)


