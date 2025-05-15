# -*- coding: utf-8 -*-
# @Author: KaiWang
# @Date:   2024-01-08 10:07:08
# @Last Modified by:   wangdaha
# @Last Modified time: 2024-01-08 14:29:11
import os,sys
import numpy as np
from functions import *


def generate_case_data(ifpath,ofpath,seed,testchem,numofchem,numofgene,cminputfile,chemfile1,chemfile2,genefile1,genefile2,genefile3,num_col,testchemfile1,testchemfile2):
	# get gene and chem
	print("get total 5000 chemicals with 10000 genes ...")
	chems,genes = get_genes_and_chem_from_total_list(testchem)
	print("get test chems and genes ...")
	chemicals = assgin_list_with_numID(shuffle_list([testchem]+shuffle_list(chems,numofchem-1,seed+2334),numofchem,seed+3312),seed+123)
	genenames = assgin_list_with_numID(shuffle_list(genes,numofgene,seed+213213),seed+123)
	print("save genes and chems into files ...")
	save_gene_chem_to_file(genenames,chemicals,ofpath)
	
	if(len(chemicals) != 0 and len(genenames) != 0):
		N = int(num_col)
		print("obtain chemfile1 content ...")
		chem1_content = obtain_chemicals_from_chemical_similarity(chemicals,chemfile1,testchemfile1)
		print("obtain chemfile2 content ...")
		chem2_content = obtain_chemicals_from_chemical_similarity(chemicals,chemfile2,testchemfile2)
		print("generate chemfile1 file ...")

		ccm1 = generate_ccm_matrix(chem1_content,chemicals)
		print("generate chemfile2 file ...")
		ccm2 = generate_ccm_matrix(chem2_content,chemicals)
		print("save to file ...")
		np.savetxt(ofpath+"chem_chem_1_matrix.txt",ccm1,delimiter="\t",fmt='%.6f')
		np.savetxt(ofpath+"chem_chem_2_matrix.txt",ccm2,delimiter="\t",fmt='%.6f')
		#-------------------
		#process gene
		#-------------------
		print("obtain genefile1 content ...")
		ggm1 = read_and_findidx_genesimfile(genefile1,genenames)
		print("obtain genefile2 content ...")
		ggm2 = read_and_findidx_genesimfile(genefile2,genenames)
		print("obtain genefile3 content ...")
		ggm3 = read_and_findidx_genesimfile(genefile3,genenames)
		print("save to file ...")
		np.savetxt(ofpath+"gene_gene_1_matrix.txt",ggm1,delimiter="\t",fmt='%.6f')
		np.savetxt(ofpath+"gene_gene_2_matrix.txt",ggm2,delimiter="\t",fmt='%.6f')
		np.savetxt(ofpath+"gene_gene_3_matrix.txt",ggm3,delimiter="\t",fmt='%.6f')
		print("obtain main content...")
		maincontent = obtain_main_content(chemicals,genenames,N,cminputfile)
		print("generate main matrix...")
		mm = generate_main_matrix(genenames,chemicals,maincontent)
		print("save to file...")
		np.savetxt(ofpath+"gene_chem_matrix.txt",mm,delimiter="\t",fmt='%.6f')
	else:
		print("no chemicals and genes selected!")

def generate_for_cmmc_run(path,testchem):
	if os.path.exists(path+"gene_chem_matrix.txt") is True and os.path.exists(path+"chem_id.txt") is True:
		print("generating data for cmmc ...")
		chemidx = get_testchem_idx(path+"chem_id.txt",testchem)
		print("generating test chem idx ...")
		generate_test_chem_idx(path+"test_chem_idx.txt",chemidx)
		idx     = list(map(int,list(chemidx.values())))
		print("generate new mat ...")
		genereate_new_matrix(path+"gene_chem_case_mat.txt",idx,path+"gene_chem_matrix.txt")
		print("Done!")
	else:
		print("Files are not exists!")
		exit()


