import os
import pandas as pd
import numpy as np
'''
    This method takes the text file that has been generated from GATK VariantToTable tool and reads in Pandas.
    It merges two columns "hiConfDeNovo" and "loConfDeNovo" into third column "denovos". Next, it makes some of the columns
    numeric because they contain both string and integers, i.e. "ExAC_ALL" column contains an allele frequency as well as "." 
    where ExAC_ALL doesn't have allele frequency for the variant.  
'''


try:
    def vcf_filter(main_dir, file_extension):
        file = os.listdir(main_dir)
        for f in file:
            if f.endswith(file_extension):
                input_file = main_dir + "/" + f
                data = pd.read_csv(input_file, sep='\t',low_memory=False,keep_default_na=False)
                print(data.info(verbose=False, memory_usage="deep"))
                data.loc[data['hiConfDeNovo'] != "NA",'hiConfDeNovo'] = 'hiConfDeNovo'
                data.loc[data['loConfDeNovo'] != "NA", 'loConfDeNovo'] = 'loConfDeNovo'
                data['denovos'] = data['hiConfDeNovo'].str.cat(data['loConfDeNovo'],sep='')
                data.loc[data['denovos'] == "NAloConfDeNovo", 'denovos'] = 'loConfDeNovo'
                data.loc[data['denovos'] == "hiConfDeNovoNA", 'denovos'] = 'hiConfDeNovo'
                data.loc[data['denovos'] == "NANA", 'denovos'] = 'NA'

                # These lines of the code converts the columns to numeric
                data['ExAC_ALL'] = pd.to_numeric(data.ExAC_ALL, errors='coerce')
                data['CADD_phred'] = pd.to_numeric(data.CADD_phred, errors='coerce')
                data['DANN_score'] = pd.to_numeric(data.DANN_score, errors='coerce')
                data['gnomAD_exome_AF_ALL'] = pd.to_numeric(data.gnomAD_exome_AF_ALL, errors='coerce')

                #denovo_file = data[(data['Func.refGene'].isin(['exonic','exonic\splicing','splicing'])) & (data['ExonicFunc.refGene'] == "nonsynonymous_SNV") & (data['AF'] < 0.001) &  (data['ExAC_ALL'] < 0.001 ) & (data['DANN_score'] > 0.9 ) & (data['CADD_phred'] > 20) & (data['DP'] > 10) &
                #(data['CLNSIG'].isin(['Pathogenic','Likely_benign','Likely_pathogenic','Uncertain_significance','Conflicting_interpretations_of_pathogenicity']))]

                # It filters the pandas data frame using the following criteria
                #denovo_file = data[(data['Func.refGene'].isin(['exonic','exonic\x3bsplicing','splicing','nc_RNA_exonic','nc_RNA_exonic\x3bsplicing','nc_RNA_splicing'])) & (data['ExonicFunc.refGene'].isin(['nonsynonymous_SNV','unknown', '.'])) & (data['AF'] < 0.001) & (data['ExAC_ALL'] < 0.001) & (data['DANN_score'] > 0.9) & (data['CADD_phred'] > 20) & (data['DP'] > 10)]
                #denovo_file = data[(data['Func.refGene'].isin(['exonic', 'exonic\x3bsplicing', 'splicing', 'nc_RNA_exonic', 'nc_RNA_exonic\x3bsplicing','nc_RNA_splicing'])) & (data['ExonicFunc.refGene'].isin(['nonsynonymous_SNV', 'unknown', '.'])) & (data['gnomAD_exome_AF_ALL'] < 0.001) & (data['ExAC_ALL'] < 0.001) & (data['DANN_score'] > 0.9) & (data['CADD_phred'] > 20) & (data['DP'] > 10)]

                denovo_file = data.loc[((data.Func_refGene == 'exonic') | (data.Func_refGene == 'exonic\x3bsplicing') | ( data.Func_refGene == 'splicing') | (data.Func_refGene == 'nc_RNA_exonic') | (data.Func_refGene =='nc_RNA_exonic\x3bsplicing') | (data.Func_refGene =='nc_RNA_splicing')) & ((data.ExonicFunc_refGene != "synonymous_SNV"))  & ((data.DANN_score > 0.9) | (data.DANN_score == '.')) & ((data.CADD_phred > 20) | (data.CADD_phred == '.')) & ((data.DP > 10))]
                data_out_file_1 = input_file.replace(file_extension, "gatk_anno_var.xlsx")
                denovo_file.to_excel(data_out_file_1)

                # Please run the following line if you want filter GATK annotated de-novos
                # gatk_file = (data[(data['denovos'].isin(['hiConfDeNovo','loConfDeNovo']))])
                # data_out_file_1 = input_file.replace(file_extension, "gatk_anno_var.xlsx")
                # gatk_file.to_excel(data_out_file_1)

                # length = (len(denovo_file.columns) - 1)
                # denovo_file['alt_gt'] = denovo_file['REF'] + "/" + denovo_file['ALT']
                # child_ID = (denovo_file.columns[length - 4])
                # mother_ID = (denovo_file.columns[length - 3])
                # father_ID = (denovo_file.columns[length - 2])
                # child_gt = denovo_file[[child_ID]]
                # mother_gt = denovo_file[[mother_ID]]
                # father_gt = denovo_file[[father_ID]]
                # denovo_file['manual_denovo'] = np.where((denovo_file['alt_gt'] == denovo_file[child_ID]) & ((denovo_file[child_ID] != denovo_file[mother_ID]) & (denovo_file[child_ID] != denovo_file[father_ID])), 'Match_ALT', 'Match_REF')
                #
                # denovo_file.drop(['hiConfDeNovo'], axis=1)
                # denovo_file.drop(['loConfDeNovo'], axis=1)
                # #out_csv = main_dir+ "/" + current_file.replace("_VTT.txt", "_Sig_var.csv")
                # data_out_file_2 = input_file.replace(file_extension,"exonic_damaging_manual_denovo.xlsx")
                # #denovo_file.to_csv(out_csv)
                # denovo_file.to_excel(data_out_file_2)


except Exception as e:
    print("An Error has occurred while running vcf_filter method. And the error message is: {}".format(e))