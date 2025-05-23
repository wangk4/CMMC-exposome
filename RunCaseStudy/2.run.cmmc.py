# -*- coding: utf-8 -*-
# @Author: KaiWang
# @Date:   2024-01-08 14:33:16
# @Last Modified by:   wangdaha
# @Last Modified time: 2024-01-08 16:22:40
import os,sys,getopt
from os import path
import argparse
import subprocess
from functions import check_path

def main(argv):
	parser = argparse.ArgumentParser(description="This code is used to run CMMC on the data generated by generate.data.py. The compiled cmmc should be in the same folder with this code.")
	parser.add_argument("-D","--DataFolder",help="Folder for storing the matrices. This is the OutputFolder in generate.data.py.")
	args=parser.parse_args()

	if(args.DataFolder ==None):
		print("Please set the Name for the testing chemical!")
		exit()
	else:
		folder = check_path(args.DataFolder)

	for i in range(0,200):
		subfolder=folder+"Set"+str(i)+"/"
		out1=subfolder+"gene_chem_hat_mat.run1.csv"
		out2=subfolder+"gene_chem_hat_mat.run2.csv"
		out3=subfolder+"gene_chem_hat_mat.run3.csv"
		out4=subfolder+"gene_chem_hat_mat.run4.csv"
		out5=subfolder+"gene_chem_hat_mat.run5.csv"
		out6=subfolder+"gene_chem_hat_mat.run6.csv"
		cmin=subfolder+"gene_chem_case_mat.txt"
		c1=subfolder+"chem_chem_1_matrix.txt"
		c2=subfolder+"chem_chem_2_matrix.txt"
		g1=subfolder+"gene_gene_1_matrix.txt"
		g2=subfolder+"gene_gene_2_matrix.txt"
		g3=subfolder+"gene_gene_3_matrix.txt"
	
		CMD1="./cmmc"+" "+cmin+" "+g1+" "+c1+" "+out1
		CMD2="./cmmc"+" "+cmin+" "+g1+" "+c2+" "+out2
		CMD3="./cmmc"+" "+cmin+" "+g2+" "+c1+" "+out3
		CMD4="./cmmc"+" "+cmin+" "+g2+" "+c2+" "+out4
		CMD5="./cmmc"+" "+cmin+" "+g3+" "+c1+" "+out5
		CMD6="./cmmc"+" "+cmin+" "+g3+" "+c2+" "+out6
		
		subprocess.run(CMD1,shell=True)
		subprocess.run(CMD2,shell=True)
		subprocess.run(CMD3,shell=True)
		subprocess.run(CMD4,shell=True)
		subprocess.run(CMD5,shell=True)
		subprocess.run(CMD6,shell=True)

if __name__ == "__main__":
	main(sys.argv[1:])

	