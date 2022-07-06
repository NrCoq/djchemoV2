from rdkit import Chem
import hashlib


class MolValidate:
	def valid_smi(smi):
		return Chem.MolFromSmiles(smi, sanitize=True) is not None

	def get_mol(smi):
		return Chem.MolFromSmiles(smi, sanitize=True)

	def get_std_smiles(mol):
		return Chem.MolToSmiles(mol, isomericSmiles=False)

	def get_unique_id(mol):
		smi = Chem.MolToSmiles(mol)
		byte_smi = smi.encode('utf-8')
		return hashlib.md5(byte_smi).hexdigest()


