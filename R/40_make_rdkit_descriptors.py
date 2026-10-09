# ============================================================
# 40_make_rdkit_descriptors.py
# ============================================================

import csv
from rdkit import Chem
from rdkit.Chem import Descriptors, Crippen

def process(input_file, output_file):
	with open(input_file, newline="", encoding="utf-8") as f:
		reader = csv.DictReader(f)
		smiles_list = [row["smiles"] for row in reader]

	with open(output_file, "w", newline="", encoding="utf-8") as f:
		writer = csv.writer(f)
		writer.writerow(["smiles", "mw_rdkit", "logp_rdkit"])

		failed = 0

		for smiles in smiles_list:
			try:
				mol = Chem.MolFromSmiles(smiles)
			except Exception:
				mol = None

			if mol is None:
				writer.writerow([smiles, "", ""])
				failed += 1
			else:
				writer.writerow([
					smiles,
					Descriptors.MolWt(mol),
					Crippen.MolLogP(mol)
				])

	return len(smiles_list), failed

n_skinpix, fail_skinpix = process(
	"data/processed/40_skinpix_smiles.csv",
	"data/processed/40_skinpix_rdkit_descriptors.csv"
)

n_huskin, fail_huskin = process(
	"data/processed/40_huskindb_smiles.csv",
	"data/processed/40_huskindb_rdkit_descriptors.csv"
)

print("SkinPiX structures:", n_skinpix)
print("SkinPiX failed:", fail_skinpix)
print("HuSkin structures:", n_huskin)
print("HuSkin failed:", fail_huskin)