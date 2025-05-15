# -*- coding: utf-8 -*-
# @Author: KaiWang
# @Date:   2024-01-08 10:05:33
# @Last Modified by:   wangdaha
# @Last Modified time: 2024-01-08 14:29:24

# this code is to generate data for CMMC run.

#python 1.generate.data.py -C ./inputfolder/chem.txt -I ./inputfolder/ -O ./outputfolder/ -s 122 -N 1

import os,sys,getopt
from os import path
import argparse
import numpy as np
from gcd import *

def main(argv):
	parser = argparse.ArgumentParser(description="This code is used to generate the case matrix datasets from processed CTD data for CMMC run.")
	parser.add_argument("-C","--ChemicalsName",help="The text file that contains the chemical you would like to test")
	parser.add_argument("-I","--InputDir",help="folder contains the data that used to generate the matrices")
	parser.add_argument("-O","--OutputDir",help="folder for storing the matrices")
	parser.add_argument("-c","--numofchem",help="default is 300",type=int)# total chem num
	parser.add_argument("-g","--numofgene",help="default is 500",type=int)#total gene num
	parser.add_argument("-s","--seed",help="random seed for generate each matrix",type=int)
	parser.add_argument("-N","--num_col",help="main data column selected, N1 or N2, use 1 for N1 as recommended",type=int)
	args=parser.parse_args()


	if(args.ChemicalsName ==None):
		print("Please set the Name for the testing chemical!")
		exit()
	else:
		chem = get_chem(args.ChemicalsName)
	if(args.InputDir == None):
		print("Please set inputDir")
		exit()
	else:
		ifpath = check_path(args.InputDir)
	if(args.OutputDir == None):
		print("this code is used to generate case study data for the manuscript ...")
		ofpath="/espresso/kaiwang/exposome/revision/case_data/"
	else:
		ofpath = check_path(args.OutputDir)

	if(args.seed == None):
		print("please set a random number for data generation. This seed can be any number.")
		exit()
	else:
		seed = args.seed

	if(args.num_col == None):
		print("Please set column number for interaction values!")
		exit()

	if(args.numofchem == None):
		print("Set the number of chemicals use the default setting, i.e., 300")
		numofchem = 300
	else:
		numofchem = args.numofchem

	if(args.numofgene == None):
		print("Set the number of genes use the default setting, i.e., 500")
		numofgene = 500
	else:
		numofgene = args.numofgene

	cminputfile = ifpath + "homo.chemical.gene.interaction.txt"

	print("generating matrices data into " + ofpath + " for " + chem + " with" + " -N " + str(args.num_col)  + " -g "+ str(numofgene) + " -c " + str(numofchem) + "...")
	# generate folders to save matrix data
	#generate_subfolder(ofpath)
	# generate chemical similarities for the test chemical
	#generate_chemical_similarities(ifpath,chem)

	chemfile1 = ifpath + "homo.chemical.mg.sim.sparse.txt"
	chemfile2 = ifpath + "homo.chemical.tt.sim.sparse.txt"
	genefile1 = ifpath + "homo.gene.sim.MF.txt"
	genefile2 = ifpath + "homo.gene.sim.BP.txt"
	genefile3 = ifpath + "homo.gene.sim.CC.txt"

	testchemfile1 = ifpath + "chem.mg.txt"
	testchemfile2 = ifpath + "chem.tt.txt"

	for i in range(0,200):
		newpath = ofpath+"Set"+str(i)+"/"
		generate_case_data(ifpath,newpath,seed+i,chem,numofchem,numofgene,cminputfile,chemfile1,chemfile2,genefile1,genefile2,genefile3,args.num_col,testchemfile1,testchemfile2)
		generate_for_cmmc_run(newpath,chem)


if __name__ == "__main__":
	main(sys.argv[1:])

