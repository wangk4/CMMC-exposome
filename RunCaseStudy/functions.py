import os,sys
import numpy as np
import itertools
import cirpy
import subprocess
import pubchempy as pcp
from rdkit import Chem
from rdkit import DataStructs
from rdkit import RDLogger
from rdkit.Chem import AllChem
from rdkit.Chem.AtomPairs import Torsions


def generate_chemical_similarities(path,chem):
	print("generating chemical similarities ...")
	op1 = open(path + "chem.mg.txt",'w')
	op2 = open(path + "chem.tt.txt",'w')
	sme = cirpy.resolve(chem,'smiles')
	totalfile = path + "homo.chemical.smile.txt"
	m1 = Chem.MolFromSmiles(sme)
	with open(totalfile) as f:
		for line in f:
			line = line.strip()
			tchem = line.split("\t")[0]
			tsme  = line.split("\t")[1]
			if tsme == "NA":
				op1.write(chem+"\t"+tchem+"\t0.0\n")
				op2.write(chem+"\t"+tchem+"\t0.0\n")
			else:
				m2    = Chem.MolFromSmiles(tsme)
				RDLogger.DisableLog('rdApp.*')
				f1    = AllChem.GetMorganFingerprintAsBitVect(m1,2)
				f2    = AllChem.GetMorganFingerprintAsBitVect(m2,2)
				mgsim = DataStructs.FingerprintSimilarity(f1,f2)
				op1.write(chem+"\t"+tchem+"\t"+str(mgsim)+"\n")

				t1 = Torsions.GetTopologicalTorsionFingerprintAsIntVect(m1)
				t2 = Torsions.GetTopologicalTorsionFingerprintAsIntVect(m2)
				ttsim = DataStructs.DiceSimilarity(t1,t2)
				op2.write(chem+"\t"+tchem+"\t"+str(ttsim)+"\n")
	op1.close()
	op2.close()
	print("chemical similarities generated!")


def get_genes_and_chem_from_total_list(testchem):
	genes = []
	chems = []
	genefile="/espresso/kaiwang/exposome/revision/gene.list.txt"
	chemfile="/espresso/kaiwang/exposome/revision/chem.list.txt"
	with open(genefile) as f:
		for line in f:
			line = line.strip()
			genes.append(line)
	with open(chemfile) as f:
		for line in f:
			line = line.strip()
			if line == testchem:
				next
			else:
				chems.append(line)
	return(chems,genes)
	
def check_path(path):
	if path.endswith("/"):
		out = path
	else:
		out = path+"/"
	return(out)

def get_idx_for_test_chem(chem,path):
	idx = -1
	if os.path.exists(path+"test_chem_idx.txt") is False:
		print("File is not exists!")
		exit()
	else:
		with open(path+"test_chem_idx.txt") as f:
			for i in f:
				i = i.strip()
				c = i.split("\t")[0]
				if c == chem:
					idx = int(i.split("\t")[1])
					break
	return(idx)

def generate_test_chem_idx(file,idx):
	out = open(file,"w")
	for k,v in idx.items():
		out.write(k+"\t"+str(v)+"\n")
	out.close()

def genereate_new_matrix(outf,idx,inf):
	mat = np.loadtxt(inf)
	mat[:,idx] = 0.0
	np.savetxt(outf,mat,delimiter="\t",fmt='%.6f')
	print("New Matrix generated!")

def get_testchem_idx(file,testchem):
	out = {}
	with open(file) as f:
		for line in f:
			line = line.strip()
			chem = line.split("\t")[0]
			idx  = line.split("\t")[1]
			if chem in testchem:
				out[chem] = idx
	return(out)


def generate_gene_and_chemical_list(seed,testchem):
	totalchem=[]
	with open("/espresso/kaiwang/exposome/data/CTD/chemical/homo.chemical.name.txt") as f:
		for line in f:
			line = line.strip()
			if line in testchem:
				next
			else:
				totalchem.append(line)
	numofchem = 5000 - len(testchem)
	outchem = shuffle_list(totalchem,numofchem,54645687)

	posgene = readchemgenefile("/espresso/kaiwang/exposome/case_study_for_paper/genes.txt")
	totgene = readchemgenefile("/espresso/kaiwang/exposome/case_study_for_paper/total.gene.txt")

	totalgene = []
	for i in totgene:
		if i in posgene:
			next
		else:
			totalgene.append(i)
	numofgene = 10000-len(posgene)
	outgene = posgene+shuffle_list(totalgene,numofgene,24543132)

	return(outchem,outgene,posgene)

def shuffle_list(inputlist, numofitem, seed):
	outlist = []
	ids = list(np.random.RandomState(seed).permutation(len(inputlist))[0:numofitem])
	for i in ids:
		outlist.append(inputlist[i])
	return(outlist)

def readchemgenefile(file):
	out = []
	with open(file) as f:
		for line in f:
			line = line.strip()
			out.append(line)
	return(out)

def assgin_list_with_numID(inputlist,seed):
	# will shuffle the list
	out = {}
	i = 0 
	inputlist = shuffle_list(inputlist,len(inputlist),seed)
	for item in inputlist:
		out[item] = i
		i+=1
	return(out)

def get_chem(file):
	out = []
	with open(file) as f:
		for line in f:
			line = line.strip()
			out.append(line)
	if len(out)>1:
		print("The file should only contain one chemical!")
		exit()
	else:
		return(out[0])

def generate_subfolder(ofpath):
	print("Generating subfolders ...")
	for i in range(0,200):
		cmd = "mkdir "+ofpath+"Set"+str(i)
		subprocess.run(cmd,shell=True)
	print("Subfolder generated!")

def obtain_chemicals_from_chemical_similarity(names,file1,file2):
	# very slow might need use parallel in future
	out = []
	with open(file1) as f:
		for line in f:
			line = line.strip()
			chem1 = line.split("\t")[0]
			chem2 = line.split("\t")[1]
			if chem1 in names and chem2 in names:
				out.append(line)
	with open(file2) as f:
		for line in f:
			line = line.strip()
			chem1 = line.split("\t")[0]
			chem2 = line.split("\t")[1]
			if chem1 in names and chem2 in names:
				out.append(line)
	return(out)

def generate_ccm_matrix(ccmsp,chemicals):
	#initialize matrix
	ccm_size = len(list(chemicals.keys()))
	ccm = np.identity(ccm_size)
	for line in ccmsp:
		chem1 = line.split("\t")[0]
		chem2 = line.split("\t")[1]
		idx1  = chemicals[chem1]
		idx2  = chemicals[chem2]
		ccm[idx1][idx2] = float(line.split("\t")[2])
		ccm[idx2][idx1] = float(line.split("\t")[2])
	return(ccm)	


def read_and_findidx_genesimfile(file,gene):
	# gene is a dict of gene id as key and number as value
	col = [] # save header for indexing
	row = [] # save the ids of select rows 
	tp_array = []
	colidx = [] # save col ids for subsetting the rows 
	with open(file) as f:
		header = f.readline().strip()
		for i in header.split():
			col.append(i.strip('"'))
		for line in f:
			line = line.strip()
			ids = line.split()[0].strip('"')
			if ids in gene:
				tp_array.append(line.split()[1:])
				row.append(ids)
	array = np.array(tp_array)
	for i in row:
		tp = col.index(i)
		colidx.append(tp)
	m = array[:,colidx]
	
	msize = len(list(gene.keys()))
	
	ggm = np.identity(msize)

	for i in itertools.combinations(gene.keys(), 2):
		rowid = i[0]
		colid = i[1]
		if(rowid in row and colid in row):
			subrowid = int(row.index(rowid))
			subcolid = int(row.index(colid))
			ggm[gene[rowid]][gene[colid]] = m[subrowid][subcolid] # upper triangle
			ggm[gene[colid]][gene[rowid]] = m[subrowid][subcolid] # lower triangle
	return(ggm)

def get_M(file,gene):
	# gene is a dict of gene id as key and number as value
	col = [] # save header for indexing
	row = [] # save the ids of select rows 
	tp_array = []
	colidx = [] # save col ids for subsetting the rows 
	with open(file) as f:
		header = f.readline().strip()
		for i in header.split():
			col.append(i.strip('"'))
		for line in f:
			line = line.strip()
			ids = line.split()[0].strip('"')
			if ids in gene:
				tp_array.append(line.split()[1:])
				row.append(ids)
	array = np.array(tp_array)
	for i in row:
		tp = col.index(i)
		colidx.append(tp)
	m = array[:,colidx]
	return(m,row)

def generate_pseudo_ggm(m,row,gene,seed,pseudo,path,ID):
	print("generate pseudo gene gene matrices ...")
	msize = len(list(gene.keys()))
	ggm = np.identity(msize)
	print("get pseudo content ...")
	mean = 0.1
	sd   = 0.4
	size = len(m[:,1])*len(m[:,1])-len(m[:,1])
	if pseudo == 1:
		pseudo_content = list(generate_truncnorm_sample(0,0.3,size,mean,sd,seed))
		pseudo_content = shuffle_list(pseudo_content,size,seed)
	elif pseudo == 2:
		pseudo_content = list(generate_truncnorm_sample(0.3,0.7,size,mean,sd,seed))
		pseudo_content = shuffle_list(pseudo_content,size,seed)
	elif pseudo == 3:
		pseudo_content = list(generate_truncnorm_sample(0.7,1,size,mean,sd,seed))
		pseudo_content = shuffle_list(pseudo_content,size,seed)
	j = 0
	for i in itertools.combinations(gene.keys(), 2):
		rowid = i[0]
		colid = i[1]
		if(rowid in row and colid in row):
			subrowid = int(row.index(rowid))
			subcolid = int(row.index(colid))
			ggm[gene[rowid]][gene[colid]] = pseudo_content[j] # upper triangle
			ggm[gene[colid]][gene[rowid]] = pseudo_content[j] # lower triangle
			j+=1
	print("writing to "+ str(ID) +" files ...")
	np.savetxt(path+"pseudo_"+str(ID)+"_gene_gene_matrix.txt",ggm,delimiter="\t",fmt='%.6f')

def obtain_main_content(chems,genes,N,file):
	ValueCol = N + 5
	out = []
	with open(file) as f:
		for line in f:
			line = line.strip()
			if line.startswith("ChemicalID"):
				next
			else:
				name = line.split("\t")[3]
				geneid = line.split("\t")[5]
				value  = line.split("\t")[ValueCol]
				if name in chems and geneid in genes:
					out.append(name + "\t" + geneid + "\t" + str(value))
	return(out)

def generate_main_matrix(genenames,chemicals,mainsp):
	row_size = len(list(genenames.keys()))
	col_size = len(list(chemicals.keys()))
	mm	   = np.zeros((row_size,col_size))
	for line in mainsp:
		chem = line.split("\t")[0]
		gene = line.split("\t")[1]
		value= line.split("\t")[2]
		row_idx = genenames[gene]
		col_idx = chemicals[chem]
		mm[row_idx][col_idx] = value
	return(mm)

def save_gene_chem_to_file(genenames,chemicals,path):
	gene = open(path+"gene_id.txt",'w')
	chem = open(path+"chem_id.txt",'w')
	for k,v in genenames.items():
		gene.write(k + "\t" + str(v) + "\n")
	for k,v in chemicals.items():
		chem.write(k + "\t" + str(v) + "\n")
	gene.close()
	chem.close()

def get_pos_neg_chems_from_gold_standard(pos,neg):
	chems = {}
	with open(pos) as f:
		for line in f:
			line = line.strip()
			chems[line] = 1
	with open(neg) as f:
		for line in f:
			line = line.strip()
			chems[line] = 1
	return(list(chems.keys()))

def get_pos_neg_genes_from_gold_standard(pos,neg,num,seed):
	genes = {}
	tp = []
	with open(pos) as f:
		for line in f:
			line = line.strip()
			genes[line] = 1
	with open(neg) as f:
		for line in f:
			line = line.strip()
			if line in genes:
				next
			else:
				tp.append(line)
	tp = shuffle_list(tp,len(tp),seed)
	for i in tp:
		genes[i] = 1
		if(len(genes.keys()) == num):
			break
	return(list(genes.keys()))

def get_extra_genes_from_cminputfile_casestudy(chems,genes,numofchem,numofgene,cminputfile,seed):
	geneout = []
	if(numofgene == 400):
		geneout = genes
		return(geneout)
	elif(numofgene>400):
		extra_gene_num = numofgene - 400
		tpgene = []
		with open(cminputfile) as f:
			next(f)
			for line in f:
				line = line.strip()
				genename = line.split("\t")[5]
				if genename not in genes:
					tpgene.append(genename)
		tpgene = list(dict.fromkeys(tpgene))
		tpgene = shuffle_list(tpgene,extra_gene_num,seed)
		geneout = tpgene+genes
		return(geneout)
	else:
		print('Set right number!')

def get_extra_chems_from_ACT_file(CaseorNot,seed,chems,numofchem):
	chemout = []
	if(numofchem == 200):
		print("Did not get any extra genes in this dataset.")
		print("167 gold standard chems are include, get extra 31 chems with BPA & BPF & BPS.")
		tp = []
		with open(CaseorNot) as f:
			for line in f:
				if line.startswith("bisphenol"):
					next
				elif line.strip() not in chems:
					tp.append(line.strip())
		tp = list(dict.fromkeys(tp))
		tp = shuffle_list(tp,31,seed)
		chemout = chems + tp + ["bisphenol A","bisphenol S"] #BPF is in the pos chem
		return(chemout)
	elif(numofchem>200):
		extra_chem_num = numofchem - len(chems) - 2
		tp = []
		with open(CaseorNot) as f:
			for line in f:
				if line.startswith("bisphenol"):
					next
				elif line.strip() not in chems:
					tp.append(line.strip())
		tp = list(dict.fromkeys(tp))
		tp = shuffle_list(tp,extra_chem_num,seed)
		chemout = chems + tp + ["bisphenol A","bisphenol S"] #BPF is in the pos chem
		return(chemout)


def get_extra_genes_chems_from_cminputfile(chems,genes,numofchem,numofgene,cminputfile,seed):
	geneout = []
	chemout = []
	if(numofgene == 400 and numofchem ==200):
		print("Did not get any extra genes in this dataset.")
		print("167 gold standard chems are include, get extra 33 chems.")
		tp = []
		with open(cminputfile) as f:
			next(f)
			for line in f:
				line = line.strip()
				chemical = line.split("\t")[3]
				if chemical	not in chems:
					tp.append(chemical)
		tp = list(set(tp))
		tp = shuffle_list(tp,33,seed)
		geneout = genes
		chemout = chems + tp
		return(chemout,geneout)
	elif(numofgene>400 and numofchem>200):
		extra_gene_num = numofgene - 400
		extra_chem_num = numofchem - len(chems)
		tpchem = []
		tpgene = []
		with open(cminputfile) as f:
			next(f)
			for line in f:
				line = line.strip()
				genename = line.split("\t")[5]
				chemical = line.split("\t")[3]
				if genename not in genes:
					tpgene.append(genename)
				if chemical	not in chems:
					tpchem.append(chemical)
		tpchem = list(set(tpchem))
		tpgene = list(set(tpgene))
		tpchem = shuffle_list(tpchem,extra_chem_num,seed)
		tpgene = shuffle_list(tpgene,extra_gene_num,seed)
		geneout = tpgene+genes
		chemout = tpchem+chems
		return(chemout,geneout)
	else:
		print("set right number.")

def get_extra_genes_chems_from_similar_files(chems,genes,numofchem,numofgene,cminputfile,seed,level):
	# the path of similar files are in the get_corresponding_level_files funtions
	# the function is return two list, i.e. chem list and gene list randomly within one level
	path = "/espresso/kaiwang/exposome/data/CTD/test.run.3/similar.file/"
	(chemfiles,genefiles) = get_corresponding_level_files(level,path)
	chemsover = list(set(get_overlap_item(path,chemfiles,"chem"))-set(chems))
	genesover = list(set(get_overlap_item(path,genefiles,"gene"))-set(genes))
	if(numofgene == 400 and numofchem ==200):
		print("Did not get any extra genes in this dataset.")
		print("167 gold standard chems are include, get extra 33 chems.")
		tp = shuffle_list(chemsover,33,seed)
		geneout = genes
		chemout = chems + tp
		return(chemout,geneout)
	elif(numofgene>400 and numofchem>200):
		extra_gene_num = numofgene - 400
		extra_chem_num = numofchem - len(chems)
		tpchem = shuffle_list(chemsover,extra_chem_num,seed)
		tpgene = shuffle_list(genesover,extra_gene_num,seed)
		geneout = tpgene+genes
		chemout = tpchem+chems
		return(chemout,geneout)
	else:
		print("set right number.")

def get_overlap_item(path,files,types):
	out = []
	if(types == "chem"):
		chem1 = get_content(path+files[0])
		chem2 = get_content(path+files[1])
		out = intersection(chem1,chem2)
	elif(types == "gene"):
		gene1 = get_content(path+files[0])
		gene2 = get_content(path+files[1])
		gene3 = get_content(path+files[2])
		out = intersection(intersection(gene1,gene2),gene3)
	return(out)

def intersection(lst1, lst2):
    return list(set(lst1) & set(lst2))

def get_content(file):
	out = []
	with open(file) as f:
		for line in f:
			line = line.strip()
			if line.startswith("Gene") or line.startswith("Chemical"):
				next
			else:
				con = line.split()[1]
				out.append(con)
	return(out)

def get_corresponding_level_files(level,path):
	if(level == 1):
		types = "low"
		print("generate ---low--- similarity data ...")
	elif(level == 2):
		types = "mid"
		print("generate ---mid--- similarity data ...")
	elif(level == 3):
		types = "high"
		print("generate ---high--- similarity data ...")
	files = os.listdir(path)
	chemfiles = []
	genefiles = []
	for i in files:
		if i.startswith("chem") and i.endswith(types+".txt"):
			chemfiles.append(i)
		elif i.startswith("gene") and i.endswith(types+".txt"):
			genefiles.append(i)
	return(chemfiles,genefiles)

def read_pos_neg_combination(pos,neg):
	poscom = getcom(pos)
	negcom = getcom(neg)
	return(poscom,negcom)

def getcom(file):
	tp = []
	with open(file) as f:
		for line in f:
			line = line.strip()
			if line.startswith("ChemicalID"):
				next
			else:
				tp.append(line.split("\t")[3] + "\t" + line.split("\t")[5])
	return(tp)

def read_gene_chem(gene,chem):
	gene = getsub(gene)
	chem = getsub(chem)
	return(gene,chem)

def getsub(file):
	tp = []
	with open(file) as f:
		for line in f:
			line = line.strip()
			sub  = line.split("\t")[0]
			tp.append(sub)
	return(tp)

def get_index(gene,chem,com,path,types):
	res = []
	if types == "pos":
		out = open(path + "pos.index.txt","w")
	if types == "neg":
		out = open(path + "neg.index.txt","w")
	out.write("chem_index\tgene_index\tchemName\tgeneid\n")
	for idxg, g in enumerate(gene):
		for idxc, c in enumerate(chem):
			tp = c + "\t" + g
			if tp in com:
				out.write(str(idxc) + "\t" + str(idxg) + "\t" + c + "\t" + g + "\n")
				res.append(str(idxc) + "\t" + str(idxg))
	out.close()
	return(res)

def get_random_test_set(posindex,negindex,path,rad,size,numofsets):
	out = []
	for i in range(numofsets):
		tpfile = path + "test.index." + str(i) + ".txt"
		test = generate_test_set(rad,i,size,posindex,negindex,tpfile)
		out.append(test)
	return(out)

def generate_test_set(rad,i,size,posindex,negindex,tpfile):
	outset = []
	out = open(tpfile,"w")
	out.write("chem_index(col)\tgene_index(row)\ttype\n")
	posids = list(np.random.RandomState(rad+i).permutation(len(posindex))[0:size])
	for i in posids:
		outset.append(posindex[i])
		out.write(posindex[i] + "\tpos\n")
	negids = list(np.random.RandomState(rad+i+2).permutation(len(negindex))[0:size])
	for i in negids:
		outset.append(negindex[i])
		out.write(negindex[i] + "\tneg\n")
	out.close()
	return(outset)

def generate_new_main(file,testset,path,val):
	for i in range(len(testset)):
		dat = np.loadtxt(file)
		for j in testset[i]:
			colidx = int(j.split("\t")[0])
			rowidx = int(j.split("\t")[1])
			dat[rowidx,colidx] = val
		np.savetxt(path+"gene_chem_matrix_replace_values_of_testdata."+str(i)+".txt",dat,delimiter="\t",fmt='%.6f')

def generate_train_test_data(mode,testsize,foldsize,posfile,negfile,geneidfile,chemidfile,ofpath,seed,mainfile):
	replace_value = 0
	if(mode == "test"):
		print("Test size will be 50x2 since this is in test mode...")
		ts = 50
		if(foldsize == None):
			print("Please set the foldsize.")
			exit()
	elif(mode == "gcm" or mode =="gct" or mode == "similar" or mode == "random"):
		if(testsize == None or foldsize == None):
			print("Please set the testsize and foldsize.")
			exit()
		else:
			ts = testsize

	print("start reading the interaction files ...")
	(poscom, negcom) = read_pos_neg_combination(posfile,negfile)
	print("start reading the gene and chem id file ...")
	(gene, chem) =read_gene_chem(geneidfile,chemidfile)
	print("start writing pos and neg index into files ...")
	posindex = get_index(gene,chem,poscom,ofpath,"pos")
	negindex = get_index(gene,chem,negcom,ofpath,"neg")
	print("generate random test set ...")
	testset = get_random_test_set(posindex,negindex,ofpath,seed,ts,foldsize)
	print("generate a new main matrix with only test item to " + str(replace_value) + " ...")
	generate_new_main(mainfile,testset,ofpath,replace_value)
	print("Done!")


def generate_test_data(path):
	print("generating files in "+path+"...")
	bp = ["bisphenol A","bisphenol S","bisphenol F"]
	idxfor = {}
	idxidx = []
	if os.path.exists(path+"chem_id.txt"):
		with open(path+"chem_id.txt") as f:
			for line in f:
				line = line.strip()
				name = line.split("\t")[0]
				idx  = line.split("\t")[1]
				if name == "bisphenol A":
					idxfor[name] = idx
				elif name == "bisphenol S":
					idxfor[name] = idx
				elif name == "bisphenol F":
					idxfor[name] = idx
	else:
		print("File is not here!")
		exit()
	idxfile = open(path+"BPA.idx.txt",'w')
	for i in bp:
		idxfile.write(i+"\t"+idxfor[i]+"\n")
		idxidx.append(int(idxfor[i]))
	idxfile.close()
	# idxidx = np.array(idxidx)
	if os.path.exists(path+"gene_chem_matrix.txt"):
		outmat = path+"BPA.original.txt"
		outtest= path+"gene_chem_matrix.test.txt"
		mat = np.loadtxt(path+"gene_chem_matrix.txt",delimiter="\t")
		tp  = np.transpose(mat)
		tp1 = tp[idxidx,]
		np.savetxt(outmat,tp1,delimiter="\t",fmt="%.6f")
		tp[idxidx[1:3],] = 0.0
		np.savetxt(outtest,np.transpose(tp),delimiter="\t",fmt="%.6f")
	else:
		print("matrix file is not here!")
		exit()
	print("Done")

