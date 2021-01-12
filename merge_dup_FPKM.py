
import os
import pandas as pd
import numpy as np
import xlrd
'''
    This method takes the file with the provided extension and reads in Pandas. It then selects the columns with the 
    duplicate values and sum up the row values. At the end, the program writes back the unique rows into new dataframe.   
'''


try:
    def merge_dup_FPKM(main_dir, file_extension):
        file = os.listdir(main_dir)
        for f in file:
            if f.endswith(file_extension):
                input_file = main_dir + "/" + f
                #data = pd.read_excel(input_file)
                data = pd.read_csv(input_file,sep=",")
                print(data.info(verbose=False, memory_usage="deep"))
                #unique_df = data.groupby(['tracking_id', 'GENE']).sum().reset_index()
                data_col = list(data.columns)
                #d = {data_col[2]: 'sum', data_col[3]: 'sum', data_col[4]: 'sum', data_col[5]: 'sum',data_col[6]:'sum',data_col[7]:'sum', data_col[8]:'sum', data_col[9]:'sum', data_col[10]:'sum',data_col[11]:'sum', data_col[12]:'sum', data_col[13]:'sum',data_col[14]:'sum'}
                d = {data_col[2]: 'sum', data_col[3]: 'sum', data_col[4]: 'sum', data_col[5]: 'sum',data_col[6]:'sum',data_col[7]:'sum', data_col[8]:'sum', data_col[9]:'sum', data_col[10]:'sum',data_col[11]:'sum'}
                unique_df = (data.groupby(['tracking_id', 'GENE'], sort=False, as_index=False).agg(d).reindex(columns=data.columns))
                print(unique_df)
                #unique_df_out = input_file.replace(".xlsx","_unique.xlsx")
                unique_df_out = input_file.replace(".csv", "_unique.xlsx")
                unique_df.to_excel(unique_df_out)
                # print(unique_df)


except Exception as e:
    print("An Error has occurred while running vcf_filter method. And the error message is: {}".format(e))