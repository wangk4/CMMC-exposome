import os,sys
import itertools
import multiprocessing
from rdkit import Chem
from rdkit import DataStructs
from rdkit import RDLogger
from rdkit.Chem import AllChem
from rdkit.Chem.AtomPairs import Torsions

def read_and_process_smilefile():
	out = {}
	file = "homo.chemical.smile.txt"
	with open(file) as f:
		for line in f:
			line = line.strip()
			tp = line.split("\t")
			if tp[1] == "NA":
				next
			else:
				if tp[0] in out.keys():
					next
				else:
					out[tp[0]]=tp[1]
	return(out)

def calculate_MorganFingerprints_similarity(c1,c2):
	m1    = Chem.MolFromSmiles(smile[c1])
	m2    = Chem.MolFromSmiles(smile[c2])
	RDLogger.DisableLog('rdApp.*')
	f1    = AllChem.GetMorganFingerprintAsBitVect(m1,2)
	f2    = AllChem.GetMorganFingerprintAsBitVect(m2,2)
	mgsim = DataStructs.FingerprintSimilarity(f1,f2)
	return(c1+"\t"+c2+"\t"+str(mgsim))

def calculate_TorsionFingerprints_similarity(c1,c2):
	m1    = Chem.MolFromSmiles(smile[c1])
	m2    = Chem.MolFromSmiles(smile[c2])
	RDLogger.DisableLog('rdApp.*')
	t1    = Torsions.GetTopologicalTorsionFingerprintAsIntVect(m1)
	t2    = Torsions.GetTopologicalTorsionFingerprintAsIntVect(m2)
	ttsim = DataStructs.DiceSimilarity(t1,t2)
	
	return(c1+"\t"+c2+"\t"+str(ttsim))

smile = read_and_process_smilefile()

if __name__ == '__main__':

	o1 = open("mt.txt",'w')
	chems = (list(smile.keys()))
	with multiprocessing.Pool(processes=10) as pool:
		out = pool.starmap(calculate_MorganFingerprints_similarity, itertools.combinations(chems, 2))
		pool.close()
		pool.join()
		for i in out:
			o1.write(i+"\n")
	o1.close()

	o2 = open("tt.txt",'w')
	with multiprocessing.Pool(processes=10) as pool:
		out = pool.starmap(calculate_TorsionFingerprints_similarity, itertools.combinations(chems, 2))
		pool.close()
		pool.join()
		for i in out:
			o2.write(i+"\n")
	o2.close()
