import os
import pandas as pd
import numpy as np
import subprocess

'''
    This method takes the vcf file which has not been annotated and runs through the GATK's genotype refinement workflow. 
    PLEASE KEEP IN MIND THAT THE DIRECTORY CONTAINING PEDIGREE FILES SHOULD BE IN THE MAIN DIRECTORY, AND THE NAME
    CONVENTION OF THE PEDIGREE FILE SHOULD BE SIMILAR AS VCF FILE EXCEPT THE EXTENSION.

'''

gatk="/Users/bowcock_lab/Desktop/Analysis_Softwares/gatk-4.1.5.0/gatk"
reference_genome="/Users/bowcock_lab/Desktop/Analysis_Softwares/Reference_Fastas/Homo_sapiens_assembly38.fasta"
high_snp = "/Users/bowcock_lab/Desktop/Analysis_Softwares/Databases/1000G_phase1.snps.high_confidence.hg38_BIALLELIC_ONLY.vcf.gz"

try:

    def index_vcf(vcf_file):
        cmd4 = gatk, "IndexFeatureFile", "-I", vcf_file
        process4 = subprocess.Popen(cmd4).wait()

    def gatk_denovo(main_dir, file_extension):
        file = os.listdir(main_dir)
        for f in file:
            if f.endswith(file_extension):
                input_file = main_dir + "/" + f
                index_file = input_file.replace(".vcf.gz",".vcf.gz.tbi")
                if not os.path.exists(index_file):
                    index_vcf(input_file)
                ped_file = main_dir + "/Pedigree_Files/" + f.replace(".vcf.gz",".ped.txt")
                cgp_file = main_dir + "/" +f.replace(".vcf.gz","_recalibratedVariants.postCGP.vcf")
                g_filtered_file = main_dir + "/" + f.replace(".vcf.gz","_recalibratedVariants.postCGP.Gfilteredcall.vcf")
                denovo_file = main_dir + "/" + f.replace(".vcf.gz","recalibratedVariants.postCGP.Gfiltered.deNovos.vcf")
                filter = "GQ < 20.0"
                cmd1 = [gatk, "CalculateGenotypePosteriors", "-R",reference_genome, "-supporting", high_snp,"-ped",ped_file,"-V",input_file ,"-O",cgp_file]
                cmd2 = [gatk, "VariantFiltration", "-R", reference_genome, "-V" , cgp_file , "-G-filter",filter,"-G-filter-name",'"lowGQ"',"-O", g_filtered_file]
                cmd3 = [gatk , "VariantAnnotator", "-R", reference_genome, "-V", g_filtered_file, "-A","PossibleDeNovo","-ped",ped_file,"-O", denovo_file]

            process1 = subprocess.Popen(cmd1).wait()
            process2 = subprocess.Popen(cmd2).wait()
            process3 = subprocess.Popen(cmd3)
            process3.communicate()


except Exception as e:
    print("An Error has occurred while running the method gatk_denovo. An error is : {}".format(e))