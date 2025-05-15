# CMMC-exposome

This site contains the benchmark datasets and C++ source code of the coupled matrix-matrix completion algorithm for exposome target interaction prediction.

## C++ Source Code

[Armadillo](http://arma.sourceforge.net/) with [OpenBLAS](https://www.openblas.net/) are needed to compile the code.

Complaing example:

<pre><code> g++ CMMC.cpp -o CMMC -std=c++11 -O3 -I/path/to/armadillo/\<include\>/ -I/path/to/openblas/\<include\>/ -L/path/to/openblas/openblas/\<lib\>/ -lopenblas </code></pre>

## Benchmark Dataset

Three different sizes of the input matrices, i.e., 200 chemicals by 400 genes (200_400), 400 chemicals by 600 genes(400_600), and 600 chemicals by 800 genes (600_800), with corresponding coupled matrices can be found in [BenchMarkdatasetsN1](https://github.com/sartorlab/CMMC-on-Exposome-Prediction/tree/main/BenchMarkdatasetsN1) and [BenchMarkdatasetsN2](https://github.com/sartorlab/CMMC-on-Exposome-Prediction/tree/main/BenchMarkdatasetsN2).

For each dataset: 
- "gene_chem_matrix.txt" is the chemical gene interaction matrix. 

- "chem_chem_1_matrix.txt" is the chemical similarites calculated based on the Morgan fingerprints.
- "chem_chem_2_matrix.txt" is the chemical similarites calculated based on the Topological Torsion fingerprints.

- "gene_gene_1_matrix.txt" is the gene similarites calculated based on GOMF.
- "gene_gene_2_matrix.txt" is the gene similarites calculated based on GOBP.
- "gene_gene_3_matrix.txt" is the gene similarites calculated based on GOCC.

- "neg.index.txt" and "pos.index.txt" contain the position indexes of true positve and true negative interactions in the chemical gene interaction matrix. 

	- For the true positive dataset, only the most confident chemical-gene pairs, i.e., those that have more than 20 interaction records (not unique PubMed IDs) were selected, which resulted in 462 chemical-gene pairs in the true positive dataset. To generate the true negative dataset, chemicals were sorted by the number genes they interact with to select for the most well-studied chemicals. Then, the top 20 chemicals with the greatest number of interacting genes were selected as chemicals for the negative dataset. In this manner, we generated a chemical-gene pair list with the 20 chemicals and the genes having no interaction evidence with these 20 chemicals as a true negative dataset.  

 - 10 "gene_chem_matrix_replace_values_of_testdata.[X].txt" files are the chemical gene interaction matrix with the randomly selected true positive interactions replaced by 0. 
 - 10 "test.index.[X].txt" files are the indexes of randomly selected true positive and true negative interactions in corresponding chemical gene interaction matrices. 

 - "chem_id.txt" and "gene_id.txt" are the chemical names and gene Entrez IDs. 

## Run Case Study

To predict target genes with our data used in the case study for a chemical, code in [RunCaseStudy](https://github.com/sartorlab/CMMC-on-Exposome-Prediction/tree/main/RunCaseStudy) can be used to generate corresponding data and run CMMC.

Before run data generation, i.e., generate the 200 different dataset similar to our case study, users can use the code in [generateChemSimilarity](https://github.com/sartorlab/CMMC-on-Exposome-Prediction/tree/main/RunCaseStudy/generateChemSimilarity) to generate the overall chemical similarities and use the code in [generateGeneSimilarity](https://github.com/sartorlab/CMMC-on-Exposome-Prediction/tree/main/RunCaseStudy/generateGeneSimilarity) to generate the overall gene similarities. Then, run [1.generate.data.py](https://github.com/sartorlab/CMMC-on-Exposome-Prediction/blob/main/RunCaseStudy/1.generate.data.py) to generate the 200 different sets. 

The name of the chemical that used in this testing should be saved in the [chem.txt](https://github.com/sartorlab/CMMC-on-Exposome-Prediction/blob/main/RunCaseStudy/chem.txt) file. 

Example of running data generation:

<pre><code> python 1.generate.data.py -C ./inputfolder/chem.txt -I ./inputfolder/ -O ./outputfolder/ -s 122 -N 1 </code></pre>

Run 

<pre><code> python 1.generate.data.py -h 
            usage: generate.data.py [-h] [-C CHEMICALSNAME] [-I INPUTDIR] [-O OUTPUTDIR]
                        [-c NUMOFCHEM] [-g NUMOFGENE] [-s SEED] [-N NUM_COL]

            This code is used to generate the case matrix datasets from processed CTD data
            for CMMC run.

            optional arguments:
              -h, --help            show this help message and exit
              -C CHEMICALSNAME, --ChemicalsName CHEMICALSNAME
                                    The text file that contains the chemical you would
                                    like to test
              -I INPUTDIR, --InputDir INPUTDIR
                                    folder contains the data that used to generate the
                                    matrices
              -O OUTPUTDIR, --OutputDir OUTPUTDIR
                                    folder for storing the matrices
              -c NUMOFCHEM, --numofchem NUMOFCHEM
                                    default is 300
              -g NUMOFGENE, --numofgene NUMOFGENE
                                    default is 500
              -s SEED, --seed SEED  random seed for generate each matrix
              -N NUM_COL, --num_col NUM_COL
                                    main data column selected, N1 or N2, use 1 for N1 as
                                    recommended

</code></pre>

for detils.

After the data generated, users can run 

<pre><code> python 2.run.cmmc.py -D ./outputfolder/ </code></pre>

to run CMMC to generate the prediction results.
