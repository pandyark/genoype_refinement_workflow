import sys
import os
import subprocess
import pandas as pd
import numpy as np
from variantToTable import variantToTable,variantToTable_with_geneModel, gnomad_ALL_edit, variantToTable_biome, main_bioME_VariantTOTable
from hardy_W_calc import _hwe_cal_ , _hwe_vcf_merge_gnomAD_exome, _hwe_vcf_merge_gnomAD_genome
from annovar_annotation import annovar_annotation
from gatk_denovo import gatk_denovo
from vcf_filter import vcf_filter
from merge_dup_FPKM import merge_dup_FPKM

main_dir = sys.argv[1]

gatk="/Users/bowcock_lab/Desktop/Analysis_Softwares/gatk-4.1.5.0/gatk"
high_snp = "/Users/bowcock_lab/Desktop/Analysis_Softwares/Databases/1000G_phase1.snps.high_confidence.hg38_BIALLELIC_ONLY.vcf.gz"
reference_genome="/Users/bowcock_lab/Desktop/Analysis_Softwares/Reference_Fastas/Homo_sapiens_assembly38.fasta"
bcftools="/Users/bowcock_lab/Desktop/Analysis_Softwares/bcftools-1.10.2/bcftools"
table_annovar="/Users/bowcock_lab/Desktop/Analysis_Softwares/annovar/table_annovar.pl"
humandb="/Users/bowcock_lab/Desktop/Analysis_Softwares/humandb_annovar"

file = os.listdir(main_dir)


def index_vcf(vcf_file):
    cmd4 = gatk, "IndexFeatureFile","-I", vcf_file
    process4 = subprocess.Popen(cmd4).wait()


denovo_file = None
current_file = None

#vcf_filter(main_dir,"VTT.txt")
#gnomad_ALL_edit(main_dir,"hg38_multianno.vcf")
#annovar_annotation(main_dir,"_Left_norm.ano.input.step2.vcf")
#variantToTable(main_dir,"_anno.hg38_multianno.vcf")
#vcf_filter(main_dir,"_VTT.txt")
#variantToTable_with_geneModel(main_dir,".vcf")
variantToTable_biome(main_dir,".vcf")
_hwe_cal_(main_dir,".vcf")
_hwe_vcf_merge_gnomAD_exome(main_dir)




# try:
#     for f in file:
#         if f.endswith(".vcf.gz"):
#             input_file = main_dir + "/" + f
#             index_file = input_file.replace(".vcf.gz",".vcf.gz.tbi")
#             if not os.path.exists(index_file):
#                 index_vcf(input_file)
#             ped_file = main_dir + "/Pedigree_Files/" + f.replace(".vcf.gz",".ped.txt")
#             cgp_file = main_dir + "/" +f.replace(".vcf.gz","_recalibratedVariants.postCGP.vcf")
#             g_filtered_file = main_dir + "/" + f.replace(".vcf.gz","_recalibratedVariants.postCGP.Gfilteredcall.vcf")
#             denovo_file = main_dir + "/" + f.replace(".vcf.gz","recalibratedVariants.postCGP.Gfiltered.deNovos.vcf")
#             filter = "GQ < 20.0"
#             cmd1 = [gatk, "CalculateGenotypePosteriors", "-R",reference_genome, "-supporting", high_snp,"-ped",ped_file,"-V",input_file ,"-O",cgp_file]
#             cmd2 = [gatk, "VariantFiltration", "-R", reference_genome, "-V" , cgp_file , "-G-filter",filter,"-G-filter-name",'"lowGQ"',"-O", g_filtered_file]
#             #cmd3 = [gatk , "VariantAnnotator", "-R", reference_genome, "-V", g_filtered_file, "-A","PossibleDeNovo","-ped",ped_file,"-O", denovo_file]
#
#         process1 = subprocess.Popen(cmd1).wait()
#         process2 = subprocess.Popen(cmd2).wait()
#         #process3 = subprocess.Popen(cmd3)
#         process2.communicate()
#
#         if f.endswith(".Gfilteredcall.vcf"):
#             input_file = main_dir + "/" + f
#             preanno_step_1 = main_dir + "/" + f.replace("recalibratedVariants.postCGP.Gfilteredcall.vcf",".ann.step1.vcf")
#             preanno_step_2 = main_dir + "/" + f.replace("recalibratedVariants.postCGP.Gfilteredcall.vcf",".Left_norm.ano.input.step2.vcf")
#             cmd4 = [bcftools,"norm","-m-",input_file,"-o",preanno_step_1,input_file]
#             cmd5 = [bcftools,"norm","-f",reference_genome,"-o",preanno_step_2,preanno_step_1]
#
#             process4 = subprocess.Popen(cmd4).wait()
#             process5 = subprocess.Popen(cmd5)
#             process5.communicate()
#
#         # For Denovo vcf file:
#         if f.endswith(".Gfiltered.deNovos.vcf"):
#             input_file = main_dir + "/" + f
#             preanno_step_1 = main_dir + "/" + f.replace("recalibratedVariants.postCGP.Gfiltered.deNovos.vcf",".ann.step1.vcf")
#             preanno_step_2 = main_dir + "/" + f.replace("recalibratedVariants.postCGP.Gfiltered.deNovos.vcf",".Left_norm.ano.input.step2.vcf")
#             cmd4 = [bcftools,"norm","-m-both","-o",preanno_step_1,input_file]
#             cmd5 = [bcftools,"norm","-f",reference_genome,"-o",preanno_step_2,preanno_step_1]
#
#             process4 = subprocess.Popen(cmd4).wait()
#             process5 = subprocess.Popen(cmd5)
#             process5.communicate()
#
#         if f.endswith(".Left_norm.ano.input.step2.vcf"):
#             input_file = main_dir + "/" + f
#             anno_output_file = main_dir + "/" + f.replace(".Left_norm.ano.input.step2.vcf","_anno")
#             cmd6 = [table_annovar,input_file,humandb,"--buildver","hg38","--outfile",anno_output_file,"--remove","--protocol","refGene,ensGene,cytoBand,clinvar_20190305,dbnsfp35a,exac03,gnomad211_exome,kaviar_20150923,avsnp150",
#                     "--operation","g,g,r,f,f,f,f,f,f","--nastring",".","--polish", "--vcfinput"]
#
#             process6 = subprocess.Popen(cmd6)
#             process6.communicate()
#
#
#         if f.endswith("_anno.hg38_multianno.vcf"):
#             anno_file = main_dir + "/" + f
#             vtt_file = main_dir + "/" + f.replace("_anno.hg38_multianno.vcf","_VTT.txt")
#         # Doesn't have genetic model annotation
#             cmd7 = [gatk,"VariantsToTable","-V",anno_file,"-F","CHROM" ,"-F","POS","-F","REF","-F","ALT","-F","TYPE","-F","Gene.refGene" ,"-F","GeneDetail.refGene",
#                           "-F","AAChange.refGene","-F" ,"Func.refGene","-F","ExonicFunc.refGene","-F","hiConfDeNovo","-F","loConfDeNovo","-F" ,"AC","-F","DP","-F","MQ","-F","QD", "-F","AF",
#                           "-F","CLNDN","-F","CLNSIG","-F","GTEx_V6p_tissue","-F","gnomAD_exome_AF_ALL","-F","gnomAD_exome_AF_afr","-F","gnomAD_exome_AF_amr","-F","gnomAD_exome_AF_asj",
#                           "-F","gnomAD_exome_AF_eas","-F","gnomAD_exome_AF_female","-F","gnomAD_exome_AF_fin","-F","gnomAD_exome_AF_male","-F","gnomAD_exome_AF_nfe",
#                           "-F","gnomAD_exome_AF_oth","-F","gnomAD_exome_AF_popmax","-F","gnomAD_exome_AF_raw","-F","gnomAD_exome_AF_sas","-F","AS_FS","-F","ExAC_AFR","-F","ExAC_ALL","-F","ExAC_AMR",
#                           "-F","ExAC_EAS","-F","ExAC_FIN","-F","ExAC_NFE","-F","ExAC_OTH","-F","ExAC_SAS","-F","Kaviar_AF","-F","ExcessHet","-F","CADD_phred","-F","DANN_score","-F","FATHMM_pred",
#                           "-F","FATHMM_score","-F","MetaLR_pred","-F","MetaLR_score","-F","MetaSVM_pred","-F","MetaSVM_score","-F","MutationTaster_pred",
#                           "-F","MutationTaster_score","-F","PROVEAN_pred","-F","PROVEAN_score","-F","Polyphen2_HDIV_pred","-F","Polyphen2_HDIV_score","-F","Polyphen2_HVAR_pred",
#                           "-F","Polyphen2_HVAR_score","-F","REVEL_rankscore","-F","REVEL_score","-F","SIFT_pred","-F","SIFT_score","-F","avsnp150","-F","fathmm-MKL_coding_pred",
#                           "-F","fathmm-MKL_coding_score","-GF","GT","-O",vtt_file]
#
#             process7 = subprocess.Popen(cmd7).wait()
#
#
#         if f.endswith("_VTT.txt"):
#             current_file = f
#             panda_file = main_dir + "/" + f
#             patient_ID = f.split('.')[0]
#             data = pd.read_csv(panda_file, sep='\t',low_memory=False,keep_default_na=False)
#             print(data.info(verbose=False, memory_usage="deep"))
#             data.loc[data['hiConfDeNovo'] != "NA",'hiConfDeNovo'] = 'hiConfDeNovo'
#             data.loc[data['loConfDeNovo'] != "NA", 'loConfDeNovo'] = 'loConfDeNovo'
#             data['denovos'] = data['hiConfDeNovo'].str.cat(data['loConfDeNovo'],sep='')
#             data.loc[data['denovos'] == "NANA", 'denovos'] = 'NA'
#             data.loc[data['denovos'] == "NAloConfDeNovo", 'denovos'] = 'loConfDeNovo'
#             data.loc[data['denovos'] == "hiConfDeNovoNA", 'denovos'] = 'hiConfDeNovo'
#             out_file = "/Users/bowcock_lab/Desktop/Test.csv"
#             # data['ExAC_ALL'] = pd.to_numeric(data.ExAC_ALL, errors='coerce')
#             # data['CADD_phred'] = pd.to_numeric(data.CADD_phred, errors='coerce')
#             # data['DANN_score'] = pd.to_numeric(data.DANN_score, errors='coerce')
#             # data['AF'] = pd.to_numeric(data.AF, errors='coerce')
#
#             # #denovo_file = data[(data['Func.refGene'].isin(['exonic','exonic\splicing','splicing'])) & (data['ExonicFunc.refGene'] == "nonsynonymous_SNV") & (data['AF'] < 0.001) &  (data['ExAC_ALL'] < 0.001 ) & (data['DANN_score'] > 0.9 ) & (data['CADD_phred'] > 20) & (data['DP'] > 10) &
#             # (data['CLNSIG'].isin(['Pathogenic','Likely_benign','Likely_pathogenic','Uncertain_significance','Conflicting_interpretations_of_pathogenicity']))]
#
#             #denovo_file = data[(data['Func.refGene'].isin(['exonic', 'exonic;splicing', 'splicing','nc_RNA_exonic','nc_RNA_exonic;splicing','nc_RNA_splicing'])) & (data['ExonicFunc.refGene'] == "nonsynonymous_SNV") & (data['AF'] < 0.001) & (data['ExAC_ALL'] < 0.001) & (data['DANN_score'] > 0.9) & (data['CADD_phred'] > 20) & (data['DP'] > 10)]
#
#             gatk_file = (data[(data['denovos'].isin(['hiConfDeNovo','loConfDeNovo']))])
#
#             # length = (len(denovo_file.columns) - 1)
#             # denovo_file['alt_gt'] = denovo_file['REF'] + "/" + denovo_file['ALT']
#             # child_ID = (denovo_file.columns[length - 4])
#             # mother_ID = (denovo_file.columns[length - 3])
#             # father_ID = (denovo_file.columns[length - 2])
#             # child_gt = denovo_file[[child_ID]]
#             # mother_gt = denovo_file[[mother_ID]]
#             # father_gt = denovo_file[[father_ID]]
#             # denovo_file['manual_denovo'] = np.where((denovo_file['alt_gt'] == denovo_file[child_ID]) & (denovo_file[child_ID] != denovo_file[mother_ID]) & (denovo_file[child_ID] != denovo_file[father_ID]), 'Match_ALT', 'Match_REF')
#             # denovo_file.drop(['hiConfDeNovo', 'loConfDeNovo'], axis=1)
#             gatk_file.drop(['hiConfDeNovo'], axis=1)
#             gatk_file.drop(['loConfDeNovo'], axis=1)
#             #out_csv = main_dir+ "/" + current_file.replace("_VTT.txt", "_Sig_var.csv")
#             gatk_out_file = main_dir+ "/" + current_file.replace("_VTT.txt", "_gatk_anno_var.xlsx")
#             #denovo_file.to_csv(out_csv)
#             gatk_file.to_excel(gatk_out_file)
#             #out_excel = main_dir + "/" + current_file.replace("_VTT.txt", "_Sig_var.xlsx")
#             #print(out_csv)
#             # dropped_col.to_excel(out_excel)
#
# except Exception as e:
#     print(e)


# anno_file = "/Users/bowcock_lab/Desktop/Denovo_VCF_Files/Family_VCFs/fam_PP028_Mend_Vio_output.vcf"
# vtt_file =  anno_file.replace(".vcf","_VTT.txt")
# # Includes genetic model annotation
# cmd7 = [gatk,"VariantsToTable","-V",anno_file,"-F","CHROM" ,"-F","POS","-F","REF","-F","ALT","-F","TYPE","-F","Gene.refGene" ,"-F","GeneDetail.refGene",
#               "-F","AAChange.refGene","-F" ,"Func.refGene","-F","ExonicFunc.refGene","-F","hiConfDeNovo","-F","loConfDeNovo","-F" ,"AC","-F","DP","-F","MQ","-F","QD", "-F","AF",
#               "-F","CLNDN","-F","CLNSIG","-F","GTEx_V6p_tissue","-F","gnomAD_exome_AF_ALL","-F","gnomAD_exome_AF_afr","-F","gnomAD_exome_AF_amr","-F","gnomAD_exome_AF_asj",
#               "-F","gnomAD_exome_AF_eas","-F","gnomAD_exome_AF_female","-F","gnomAD_exome_AF_fin","-F","gnomAD_exome_AF_male","-F","gnomAD_exome_AF_nfe",
#               "-F","gnomAD_exome_AF_oth","-F","gnomAD_exome_AF_popmax","-F","gnomAD_exome_AF_raw","-F","gnomAD_exome_AF_sas","-F","AS_FS","-F","ExAC_AFR","-F","ExAC_ALL","-F","ExAC_AMR",
#               "-F","ExAC_EAS","-F","ExAC_FIN","-F","ExAC_NFE","-F","ExAC_OTH","-F","ExAC_SAS","-F","Kaviar_AF","-F","ExcessHet","-F","CADD_phred","-F","DANN_score","-F","FATHMM_pred",
#               "-F","FATHMM_score","-F","MetaLR_pred","-F","MetaLR_score","-F","MetaSVM_pred","-F","MetaSVM_score","-F","MutationTaster_pred",
#               "-F","MutationTaster_score","-F","PROVEAN_pred","-F","PROVEAN_score","-F","Polyphen2_HDIV_pred","-F","Polyphen2_HDIV_score","-F","Polyphen2_HVAR_pred",
#               "-F","Polyphen2_HVAR_score","-F","REVEL_rankscore","-F","REVEL_score","-F","SIFT_pred","-F","SIFT_score","-F","avsnp150","-F","fathmm-MKL_coding_pred",
#               "-F","fathmm-MKL_coding_score","-GF","GT","-O",vtt_file]
#
# # cmd7 = [gatk,"VariantsToTable","-V",anno_file,"-F","CHROM" ,"-F","POS","-F","REF","-F","ALT","-F","TYPE","-F","Gene.refGene" ,"-F","GeneDetail.refGene",
# #               "-F","AAChange.refGene","-F" ,"Func.refGene","-F","ExonicFunc.refGene","-F","hiConfDeNovo","-F","loConfDeNovo", "-F", "GeneticModels","-F" ,"AC","-F","DP","-F","MQ","-F","QD", "-F","AF",
# #               "-F","CLNDN","-F","CLNSIG","-F","GTEx_V6p_tissue","-F","gnomAD_exome_AF_ALL","-F","gnomAD_exome_AF_afr","-F","gnomAD_exome_AF_amr","-F","gnomAD_exome_AF_asj",
# #               "-F","gnomAD_exome_AF_eas","-F","gnomAD_exome_AF_female","-F","gnomAD_exome_AF_fin","-F","gnomAD_exome_AF_male","-F","gnomAD_exome_AF_nfe",
# #               "-F","gnomAD_exome_AF_oth","-F","gnomAD_exome_AF_popmax","-F","gnomAD_exome_AF_raw","-F","gnomAD_exome_AF_sas","-F","AS_FS","-F","ExAC_AFR","-F","ExAC_ALL","-F","ExAC_AMR",
# #               "-F","ExAC_EAS","-F","ExAC_FIN","-F","ExAC_NFE","-F","ExAC_OTH","-F","ExAC_SAS","-F","Kaviar_AF","-F","ExcessHet","-F","CADD_phred","-F","DANN_score","-F","FATHMM_pred",
# #               "-F","FATHMM_score","-F","MetaLR_pred","-F","MetaLR_score","-F","MetaSVM_pred","-F","MetaSVM_score","-F","MutationTaster_pred",
# #               "-F","MutationTaster_score","-F","PROVEAN_pred","-F","PROVEAN_score","-F","Polyphen2_HDIV_pred","-F","Polyphen2_HDIV_score","-F","Polyphen2_HVAR_pred",
# #               "-F","Polyphen2_HVAR_score","-F","REVEL_rankscore","-F","REVEL_score","-F","SIFT_pred","-F","SIFT_score","-F","avsnp150","-F","fathmm-MKL_coding_pred",
# #               "-F","fathmm-MKL_coding_score","-GF","GT","-O",vtt_file]
#
# process7 = subprocess.Popen(cmd7).wait()



# panda_file = "/Users/bowcock_lab/Desktop/enGenome_Denovo/fam_PP028_genmod_Affected_Dad_VTT.txt"
# data = pd.read_csv(panda_file, sep='\t',low_memory=False,keep_default_na=False)
# print(data.info(verbose=False, memory_usage="deep"))
# data.loc[data['hiConfDeNovo'] != "NA",'hiConfDeNovo'] = 'hiConfDeNovo'
# data.loc[data['loConfDeNovo'] != "NA", 'loConfDeNovo'] = 'loConfDeNovo'
# data['denovos'] = data['hiConfDeNovo'].str.cat(data['loConfDeNovo'],sep='')
# data.loc[data['denovos'] == "NAloConfDeNovo", 'denovos'] = 'loConfDeNovo'
# data.loc[data['denovos'] == "hiConfDeNovoNA", 'denovos'] = 'hiConfDeNovo'
# data.loc[data['denovos'] == "NANA", 'denovos'] = 'NA'
# #out_file = "/Users/bowcock_lab/Desktop/Denovo_VCF_Files/Family_VCFs/MSPP001_Mend_Vio_gatk_anno_var.txt"
# data['ExAC_ALL'] = pd.to_numeric(data.ExAC_ALL, errors='coerce')
# data['CADD_phred'] = pd.to_numeric(data.CADD_phred, errors='coerce')
# data['DANN_score'] = pd.to_numeric(data.DANN_score, errors='coerce')
# data['AF'] = pd.to_numeric(data.AF, errors='coerce')
#
# # #denovo_file = data[(data['Func.refGene'].isin(['exonic','exonic\splicing','splicing'])) & (data['ExonicFunc.refGene'] == "nonsynonymous_SNV") & (data['AF'] < 0.001) &  (data['ExAC_ALL'] < 0.001 ) & (data['DANN_score'] > 0.9 ) & (data['CADD_phred'] > 20) & (data['DP'] > 10) &
# # (data['CLNSIG'].isin(['Pathogenic','Likely_benign','Likely_pathogenic','Uncertain_significance','Conflicting_interpretations_of_pathogenicity']))]
#
# #denovo_file = data[(data['Func.refGene'].isin(['exonic', 'exonic;splicing', 'splicing','nc_RNA_exonic','nc_RNA_exonic;splicing','nc_RNA_splicing'])) & (data['ExonicFunc.refGene'] == "nonsynonymous_SNV") & (data['AF'] < 0.001) & (data['ExAC_ALL'] < 0.001) & (data['DANN_score'] > 0.9) & (data['CADD_phred'] > 20) & (data['DP'] > 10)]
#
# #gatk_file = (data[(data['denovos'].isin(['hiConfDeNovo','loConfDeNovo']))])
#
# # length = (len(denovo_file.columns) - 1)
# # denovo_file['alt_gt'] = denovo_file['REF'] + "/" + denovo_file['ALT']
# # child_ID = (denovo_file.columns[length - 4])
# # mother_ID = (denovo_file.columns[length - 3])
# # father_ID = (denovo_file.columns[length - 2])
# # child_gt = denovo_file[[child_ID]]
# # mother_gt = denovo_file[[mother_ID]]
# # father_gt = denovo_file[[father_ID]]
# # denovo_file['manual_denovo'] = np.where((denovo_file['alt_gt'] == denovo_file[child_ID]) & (denovo_file[child_ID] != denovo_file[mother_ID]) & (denovo_file[child_ID] != denovo_file[father_ID]), 'Match_ALT', 'Match_REF')
# # denovo_file.drop(['hiConfDeNovo', 'loConfDeNovo'], axis=1)
# data.drop(['hiConfDeNovo'], axis=1)
# data.drop(['loConfDeNovo'], axis=1)
# #out_csv = main_dir+ "/" + current_file.replace("_VTT.txt", "_Sig_var.csv")
# data_out_file = panda_file.replace("_VTT.txt","_Inh_Pat.xlsx")
# #denovo_file.to_csv(out_csv)
# data.to_excel(data_out_file)
# #out_excel = main_dir + "/" + current_file.replace("_VTT.txt", "_Sig_var.xlsx")
# #print(out_csv)
# # dropped_col.to_excel(out_excel)










