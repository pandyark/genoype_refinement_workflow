import os
import subprocess

import pandas as pd

vcftools="/Users/bowcock_lab/Desktop/Analysis_Softwares/vcftools-master/src/cpp/vcftools"

try:
    def _hwe_cal_(main_dir,file_extension):
        if file_extension == ".vcf.gz":
            file = os.listdir(main_dir)
            for f in file:
                if f.endswith(file_extension):
                    hwe_input = main_dir + "/" + f
                    hwe_test_file = hwe_input.replace(".vcf", "_HWE")
                    logfile = hwe_test_file.replace(".hwe", "_logfile.txt")
                    hwe_cmd = [vcftools,"--gzvcf",hwe_input,"--hardy","--out",hwe_test_file]
                    hwe_process = subprocess.Popen(hwe_cmd, stdout=open(logfile,'w'))
                    hwe_process.wait()

        if file_extension == ".vcf":
            file = os.listdir(main_dir)
            for f in file:
                if f.endswith(file_extension):
                    hwe_input = main_dir + "/" + f
                    hwe_test_file = hwe_input.replace(".vcf", "_HWE")
                    logfile = hwe_test_file.replace(".hwe", "_logfile.txt")
                    hwe_cmd = [vcftools, "--vcf", hwe_input, "--hardy", "--out", hwe_test_file]
                    hwe_process = subprocess.Popen(hwe_cmd, stdout=open(logfile, 'w'))
                    hwe_process.wait()


    def _hwe_vcf_merge_gnomAD_exome(main_dir):
        file = os.listdir(main_dir)
        for f in file:
            if f.endswith("_VTT.txt"):
                vtt_file = main_dir + "/" + f
            if f.endswith(".hwe"):
                hwe_file = main_dir + "/" + f
        vtt_df1 = pd.read_csv(vtt_file, sep='\t', low_memory=False, keep_default_na=False)
        hwe_df2 = pd.read_csv(hwe_file, sep='\t', low_memory=False, keep_default_na=False)
        hwe_df2.drop(hwe_df2.columns[[0, 2, 3, 5, 6, 7]], axis=1, inplace=True)
        merge_file = (pd.merge(vtt_df1, hwe_df2, on='POS', how='left'))
        merge_file['AF'] = pd.to_numeric(merge_file.AF, errors='coerce')
        # This statement should be used for the original Biome file that was provided in rg_bowcoa01 directory
        #merge_file['gnomAD_genome_NFE'] = pd.to_numeric(merge_file.gnomAD_genome_NFE, errors='coerce')

        # This statement should be used when the vcf files have been annotated with new annovar exome database
        merge_file['gnomAD_exome_AF_nfe'] = pd.to_numeric(merge_file.gnomAD_exome_AF_nfe, errors='coerce')
        #merge_file['OR'] = merge_file['AF'] / merge_file['gnomAD_exome_AF_nfe']
        merge_file['Re-Calculated_AF'] = (2 * merge_file['HOM-VAR'] + 1 * merge_file['HET']) / (2 * merge_file['NCALLED'])
        merge_file['OR'] = merge_file['Re-Calculated_AF'] / merge_file['gnomAD_exome_AF_nfe']
        merge_output = vtt_file.replace(".txt", "_final_Merge.xlsx")
        merge_file.to_excel(merge_output)

    def _hwe_vcf_merge_gnomAD_genome(main_dir):
        file = os.listdir(main_dir)
        for f in file:
            if f.endswith("_VTT.txt"):
                vtt_file = main_dir + "/" + f
            if f.endswith(".hwe"):
                hwe_file = main_dir + "/" + f
        vtt_df1 = pd.read_csv(vtt_file, sep='\t', low_memory=True, keep_default_na=False)
        hwe_df2 = pd.read_csv(hwe_file, sep='\t', low_memory=True, keep_default_na=False)
        hwe_df2.drop(hwe_df2.columns[[0, 2, 3, 5, 6, 7]], axis=1, inplace=True)
        merge_file = (pd.merge(vtt_df1, hwe_df2, on='POS', how='left'))
        merge_file['AF'] = pd.to_numeric(merge_file.AF, errors='coerce')
        # This statement should be used for the original Biome file that was provided in rg_bowcoa01 directory
        #merge_file['gnomAD_genome_NFE'] = pd.to_numeric(merge_file.gnomAD_genome_NFE, errors='coerce')

        # This statement should be used when the vcf files have been annotated with new annovar exome database
        merge_file['gnomAD_genome_NFE'] = pd.to_numeric(merge_file.gnomAD_genome_NFE, errors='coerce')
        merge_file['Re-calculated_AF'] = (2 * merge_file['HOM-VAR'] + 1 * merge_file['HET']) / (2 * merge_file['NCALLED'])
        merge_file['OR'] = merge_file['Re-calculated_AF'] / merge_file['gnomAD_genome_NFE']
        #merge_file['Re-calculated_AF'] = (2 * merge_file['HOM-VAR'] + 1 * merge_file['HET']) / (2 * merge_file['NCALLED'])
        merge_output = vtt_file.replace(".txt", "_final_Merge.xlsx")
        merge_file.to_excel(merge_output)

except Exception as e:
    print("An error has occurred while running the hew_cal method. The error is: ".format(e))
