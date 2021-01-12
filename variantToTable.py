import os
import subprocess

gatk="/Users/bowcock_lab/Desktop/Analysis_Softwares/gatk-4.1.5.0/gatk"


try:
    def variantToTable(main_dir, file_extension):
        file = os.listdir(main_dir)
        for f in file:
            if f.endswith(file_extension):
                anno_file = main_dir + "/" + f
                vtt_file = main_dir + "/" + f.replace(file_extension,"_VTT.txt")
            # Doesn't have genetic model annotation
                cmd7 = [gatk,"VariantsToTable","-V",anno_file,"-F","CHROM" ,"-F","POS","-F","REF","-F","ALT","-F","TYPE","-F","Gene.refGene" ,"-F","GeneDetail.refGene",
                              "-F","AAChange.refGene","-F" ,"Func.refGene","-F","ExonicFunc.refGene","-F","hiConfDeNovo","-F","loConfDeNovo","-F" ,"AC","-F","DP","-F","MQ","-F","QD", "-F","AF",
                              "-F","CLNDN","-F","CLNSIG","-F","GTEx_V6p_tissue","-F","gnomAD_exome_AF_ALL","-F","gnomAD_exome_AF_afr","-F","gnomAD_exome_AF_amr","-F","gnomAD_exome_AF_asj",
                              "-F","gnomAD_exome_AF_eas","-F","gnomAD_exome_AF_female","-F","gnomAD_exome_AF_fin","-F","gnomAD_exome_AF_male","-F","gnomAD_exome_AF_nfe",
                              "-F","gnomAD_exome_AF_oth","-F","gnomAD_exome_AF_popmax","-F","gnomAD_exome_AF_raw","-F","gnomAD_exome_AF_sas","-F","AS_FS","-F","ExAC_AFR","-F","ExAC_ALL","-F","ExAC_AMR",
                              "-F","ExAC_EAS","-F","ExAC_FIN","-F","ExAC_NFE","-F","ExAC_OTH","-F","ExAC_SAS","-F","Kaviar_AF","-F","ExcessHet","-F","CADD_phred","-F","DANN_score","-F","FATHMM_pred",
                              "-F","FATHMM_score","-F","MetaLR_pred","-F","MetaLR_score","-F","MetaSVM_pred","-F","MetaSVM_score","-F","MutationTaster_pred",
                              "-F","MutationTaster_score","-F","PROVEAN_pred","-F","PROVEAN_score","-F","Polyphen2_HDIV_pred","-F","Polyphen2_HDIV_score","-F","Polyphen2_HVAR_pred",
                              "-F","Polyphen2_HVAR_score","-F","REVEL_rankscore","-F","REVEL_score","-F","SIFT_pred","-F","SIFT_score","-F","avsnp150","-F","fathmm-MKL_coding_pred",
                              "-F","fathmm-MKL_coding_score","-GF","GT","-O",vtt_file]

                process7 = subprocess.Popen(cmd7).wait()


    def variantToTable_with_geneModel(main_dir, file_extension):
        file = os.listdir(main_dir)
        for f in file:
            if f.endswith(file_extension):
                anno_file = main_dir + "/" + f
                vtt_file = main_dir + "/" + f.replace(file_extension, "_geneMod_VTT.txt")
                # Includes genetic model annotation
                cmd7 = [gatk,"VariantsToTable","-V",anno_file,"-F","CHROM" ,"-F","POS","-F","REF","-F","ALT","-F","TYPE","-F","Gene.refGene" ,"-F","GeneDetail.refGene",
                              "-F","AAChange.refGene","-F" ,"Func.refGene","-F","ExonicFunc.refGene","-F","hiConfDeNovo","-F","loConfDeNovo","-F", "GeneticModels","-F" ,"AC","-F","DP","-F","MQ","-F","QD", "-F","AF",
                              "-F","CLNDN","-F","CLNSIG","-F","GTEx_V6p_tissue","-F","gnomAD_exome_AF_ALL","-F","gnomAD_exome_AF_afr","-F","gnomAD_exome_AF_amr","-F","gnomAD_exome_AF_asj",
                              "-F","gnomAD_exome_AF_eas","-F","gnomAD_exome_AF_female","-F","gnomAD_exome_AF_fin","-F","gnomAD_exome_AF_male","-F","gnomAD_exome_AF_nfe",
                              "-F","gnomAD_exome_AF_oth","-F","gnomAD_exome_AF_popmax","-F","gnomAD_exome_AF_raw","-F","gnomAD_exome_AF_sas","-F","AS_FS","-F","ExAC_AFR","-F","ExAC_ALL","-F","ExAC_AMR",
                              "-F","ExAC_EAS","-F","ExAC_FIN","-F","ExAC_NFE","-F","ExAC_OTH","-F","ExAC_SAS","-F","Kaviar_AF","-F","ExcessHet","-F","CADD_phred","-F","DANN_score","-F","FATHMM_pred",
                              "-F","FATHMM_score","-F","MetaLR_pred","-F","MetaLR_score","-F","MetaSVM_pred","-F","MetaSVM_score","-F","MutationTaster_pred",
                              "-F","MutationTaster_score","-F","PROVEAN_pred","-F","PROVEAN_score","-F","Polyphen2_HDIV_pred","-F","Polyphen2_HDIV_score","-F","Polyphen2_HVAR_pred",
                              "-F","Polyphen2_HVAR_score","-F","REVEL_rankscore","-F","REVEL_score","-F","SIFT_pred","-F","SIFT_score","-F","avsnp150","-F","fathmm-MKL_coding_pred",
                              "-F","fathmm-MKL_coding_score","-GF","GT","-O",vtt_file]
                process7 = subprocess.Popen(cmd7).wait()

    def gnomad_ALL_edit (main_dir, file_extension):
        file = os.listdir(main_dir)
        for f in file:
            if f.endswith(file_extension):
                non_editate_file = main_dir + "/" + f
                cmd_gno = ["sed", "-i" ,'"s/[[:<:]]gnomAD_exome_AF[[:>:]]/gnomAD_exome_AF_ALL/g"', non_editate_file]
                process_gnom = subprocess.Popen(cmd_gno).wait()


    def variantToTable_biome(main_dir, file_extension):
        file = os.listdir(main_dir)
        for f in file:
            if f.endswith(file_extension):
                anno_file = main_dir + "/" + f
                vtt_file = main_dir + "/" + f.replace(file_extension, "_VTT.txt")
                biome_cmd = [gatk, "VariantsToTable", "-V", anno_file, "-F", "CHROM", "-F", "POS", "-F", "REF", "-F", "ALT","-F", "TYPE", "-F", "Gene.refGene", "-F", "GeneDetail.refGene",
                             "-F", "AAChange.refGene", "-F", "Func.refGene", "-F", "ExonicFunc.refGene","-F", "GeneticModels", "-F", "AC", "-F", "DP", "-F", "MQ",
                             "-F", "QD", "-F", "AF","-F", "CLNDN", "-F", "CLNSIG", "-F", "GTEx_V6p_tissue", "-F","1000g2015aug_ALL.sites","-F", "1000g2015aug_AFR.sites","-F","1000g2015aug_AMR.sites",
                             "-F", "1000g2015aug_EAS.sites", "-F", "1000g2015aug_EUR.sites","-F","1000g2015aug_SAS.sites", "-F", "gnomAD_exome_AF_ALL", "-F","gnomAD_exome_AF_afr", "-F", "gnomAD_exome_AF_amr",
                             "-F", "gnomAD_exome_AF_asj","-F", "gnomAD_exome_AF_eas", "-F", "gnomAD_exome_AF_female", "-F", "gnomAD_exome_AF_fin", "-F","gnomAD_exome_AF_male", "-F", "gnomAD_exome_AF_nfe",
                             "-F", "gnomAD_exome_AF_oth", "-F", "gnomAD_exome_AF_popmax", "-F", "gnomAD_exome_AF_raw", "-F","gnomAD_exome_AF_sas", "-F", "AS_FS", "-F", "ExAC_AFR", "-F", "ExAC_ALL",
                             "-F", "ExAC_AMR","-F", "ExAC_EAS", "-F", "ExAC_FIN", "-F", "ExAC_NFE", "-F", "ExAC_OTH", "-F", "ExAC_SAS", "-F","Kaviar_AF", "-F", "ExcessHet", "-F", "CADD_phred",
                             "-F", "DANN_score", "-F", "FATHMM_pred","-F", "FATHMM_score", "-F", "MetaLR_pred", "-F", "MetaLR_score", "-F", "MetaSVM_pred", "-F","MetaSVM_score", "-F", "MutationTaster_pred",
                             "-F", "MutationTaster_score", "-F", "PROVEAN_pred", "-F", "PROVEAN_score", "-F","Polyphen2_HDIV_pred", "-F", "Polyphen2_HDIV_score", "-F", "Polyphen2_HVAR_pred",
                             "-F", "Polyphen2_HVAR_score", "-F", "REVEL_rankscore", "-F", "REVEL_score", "-F", "SIFT_pred","-F", "SIFT_score", "-F", "avsnp150", "-F", "fathmm-MKL_coding_pred",
                             "-F", "fathmm-MKL_coding_score", "-F","HET","-F","HOM-REF","-F","HOM-VAR","-F","NO-CALL","-F","VAR","-F","NSAMPLE","-F","NCALLED","-O",vtt_file]

                biome_process = subprocess.Popen(biome_cmd).wait()

    def main_bioME_VariantTOTable(main_dir, file_exntension):
        file = os.listdir(main_dir)
        for f in file:
            if f.endswith(file_exntension):
                anno_file = main_dir + "/" + f
                vtt_file = main_dir + "/" + f.replace(file_exntension, "_VTT.txt")
                cmd6 = [gatk, "VariantsToTable", "-V", anno_file, "-F", "CHROM", "-F", "POS", "-F", "REF", "-F", "ALT","-F", "TYPE", "-F", "Gene.refGene",
                        "-F", "GeneDetail.refGene", "-F", "AAChange.refGene", "-F", "Func.refGene", "-F","ExonicFunc.refGene", "-F", "AC", "-F", "DP",
                        "-F", "MQ", "-F", "QD", "-F", "AF", "-F", "CLNDN", "-F", "CLNSIG", "-F","1000g2015aug_ALL.sites","-F", "1000g2015aug_AFR.sites",
                        "-F","1000g2015aug_AMR.sites","-F", "1000g2015aug_EAS.sites", "-F", "1000g2015aug_EUR.sites","-F","1000g2015aug_SAS.sites",
                        "-F", "gnomAD_genome_ALL", "-F", "gnomAD_genome_AFR", "-F", "gnomAD_genome_AMR","-F", "gnomAD_genome_ASJ", "-F", "gnomAD_genome_EAS",
                        "-F", "AF_female", "-F","gnomAD_genome_FIN", "-F", "AF_male", "-F", "gnomAD_genome_NFE", "-F", "gnomAD_genome_OTH", "-F", "AF_popmax",
                        "-F", "AF_raw", "-F", "AF_sas", "-F", "AS_FS","-F", "ExAC_AFR", "-F", "ExAC_ALL", "-F", "ExAC_AMR",
                        "-F", "ExAC_EAS", "-F", "ExAC_FIN", "-F", "ExAC_NFE", "-F", "ExAC_OTH", "-F", "ExAC_SAS", "-F","Kaviar_AF", "-F", "ExcessHet",
                        "-F", "CADD_phred", "-F", "DANN_score", "-F", "FATHMM_pred", "-F", "FATHMM_score", "-F","MetaLR_pred", "-F", "MetaLR_score",
                        "-F", "MetaSVM_pred", "-F", "MetaSVM_score", "-F", "MutationTaster_pred", "-F","MutationTaster_score", "-F", "PROVEAN_pred",
                        "-F", "PROVEAN_score", "-F", "Polyphen2_HDIV_pred", "-F", "Polyphen2_HDIV_score", "-F","Polyphen2_HVAR_pred", "-F", "Polyphen2_HVAR_score",
                        "-F", "SIFT_pred", "-F", "SIFT_score", "-F", "avsnp150", "-F", "fathmm-MKL_coding_pred", "-F","fathmm-MKL_coding_score",
                        "-F", "HET","-F", "HOM-REF", "-F", "HOM-VAR", "-F", "NO-CALL", "-F", "VAR", "-F", "NSAMPLE", "-F","NCALLED", "-O", vtt_file]
                process6 = subprocess.Popen(cmd6)
                process6.communicate()

except Exception as e:
    print("An Error has occurred while running the method variantToTable. An Error is : {}".format(e))