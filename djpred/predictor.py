from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors, Descriptors3D, rdMolDescriptors
from rdkit.Chem.EState import EState_VSA
from rdkit import DataStructs
import numpy as np
import os
from catboost import CatBoostClassifier, CatBoostRegressor


def predict(smiles):
	# set random seed for reproducibility
	from numpy.random import seed
	seed(8)

	BASE_FOLDER = os.getcwd() + '/djpred'

	mol = Chem.MolFromSmiles(smiles)
	smi = Chem.MolToSmiles(mol, isomericSmiles=False)

	results = {

	}

	properties = ["NR-AhR", "NR-AR", "NR-AR-LBD", "NR-Aromatase", "NR-ER",
				  "NR-ER-LBD", "NR-PPAR-gamma", "SR-ARE", "SR-ATAD5", "SR-HSE", "SR-MMP",
				  "SR-p53"]

	fname = BASE_FOLDER+"/nn_models/catboost_toxld50oralrats"

	m = CatBoostRegressor()
	m.load_model(fname=fname)
	descs = calc_descs(smiles)
	fp = calc_fp(smiles)
	results["LD50"] = pow(1.7, m.predict(descs))

	for prop in properties:
		m = CatBoostRegressor()
		m.load_model(BASE_FOLDER+"/nn_models/"+"catboost_"+prop)
		res = m.predict(fp)
		results[prop] = res > 0.5
	return results


def calc_fp(smiles):
	m = Chem.MolFromSmiles(smiles)
	fp = Chem.RDKFingerprint(Chem.MolFromSmiles(smiles))
	return np.array(list(fp))


def calc_descs(smiles):
	# There are 263 separate values here
	# + 2048 long Morgan fingerprint
	m = Chem.MolFromSmiles(smiles)
	descs = []


	tpsa = rdMolDescriptors.CalcTPSA(m)
	descs.append(tpsa)

	peoe = rdMolDescriptors.PEOE_VSA_(m)
	for x in peoe:
		descs.append(x)


	smr = rdMolDescriptors.SMR_VSA_(m)
	for x in smr:
		descs.append(x)

	slogp = rdMolDescriptors.SlogP_VSA_(m)
	for x in slogp:
		descs.append(x)

	evsa = EState_VSA.EState_VSA_(m)
	for x in evsa:
		descs.append(x)

	autocorr = rdMolDescriptors.CalcAUTOCORR2D(m)
	for x in autocorr:
		descs.append(x)

	molwt = Descriptors.ExactMolWt(m)
	descs.append(molwt)
	numAlRings = rdMolDescriptors.CalcNumAliphaticRings(m)
	descs.append(numAlRings)
	numArRings = rdMolDescriptors.CalcNumAromaticRings(m)
	descs.append(numArRings)
	numSatRings = rdMolDescriptors.CalcNumSaturatedRings(m)
	descs.append(numSatRings)
	armoCarbcycles = rdMolDescriptors.CalcNumAromaticCarbocycles(m)
	descs.append(armoCarbcycles)
	armoHetcycles = rdMolDescriptors.CalcNumAromaticHeterocycles(m)
	descs.append(armoHetcycles)
	minabpc = Descriptors.MinAbsPartialCharge(m)
	descs.append(minabpc)
	maxabpc = Descriptors.MaxAbsPartialCharge(m)
	descs.append(maxabpc)
	maxpc = Descriptors.MaxPartialCharge(m)
	descs.append(maxpc)
	minpc = Descriptors.MinPartialCharge(m)
	descs.append(minpc)

	chi0n = rdMolDescriptors.CalcChi0n(m)
	descs.append(chi0n)
	chi1n = rdMolDescriptors.CalcChi1n(m)
	descs.append(chi1n)
	chi2n = rdMolDescriptors.CalcChi2n(m)
	descs.append(chi2n)
	chi3n = rdMolDescriptors.CalcChi3n(m)
	descs.append(chi3n)
	chi4n = rdMolDescriptors.CalcChi4n(m)
	descs.append(chi4n)

	chi0v = rdMolDescriptors.CalcChi0v(m)
	descs.append(chi0v)
	chi1v = rdMolDescriptors.CalcChi1v(m)
	descs.append(chi1v)
	chi2v = rdMolDescriptors.CalcChi2v(m)
	descs.append(chi2v)
	chi3v = rdMolDescriptors.CalcChi3v(m)
	descs.append(chi3v)
	chi4v = rdMolDescriptors.CalcChi4v(m)
	descs.append(chi4v)

	fp = AllChem.GetMorganFingerprintAsBitVect(m, 3, 2048)
	for i in fp:
		descs.append(int(i))

	return np.array(descs)



amph_b_smi = "CC1C=CC=CC=CC=CC=CC=CC=CC(CC2C(C(CC(O2)(CC(CC(C(CCC(CC(CC(=O)OC(C(C1O)C)C)O)O)O)O)O)O)O)C(=O)O)OC3C(C(C(C(O3)C)O)N)O"

res = predict(amph_b_smi)

print(res)