import os
import subprocess



reference_genome="/Users/bowcock_lab/Desktop/Analysis_Softwares/Reference_Fastas/Homo_sapiens_assembly38.fasta"
bcftools="/Users/bowcock_lab/Desktop/Analysis_Softwares/bcftools-1.10.2/bcftools"
table_annovar="/Users/bowcock_lab/Desktop/Analysis_Softwares/annovar/table_annovar.pl"
humandb="/Users/bowcock_lab/Desktop/Analysis_Softwares/humandb_annovar"


try:
    def annovar_annotation(main_dir, file_extension):
        file = os.listdir(main_dir)
        for f in file:
            if f.endswith(file_extension):
                # Please note that the ideal file extension should be "recalibratedVariants.postCGP.Gfiltered.deNovos.vcf"
                input_file = main_dir + "/" + f
                preanno_step_1 = main_dir + "/" + f.replace(file_extension,"_preanno.step1.vcf")
                preanno_step_2 = main_dir + "/" + f.replace(file_extension,"_Left_norm.ano.input.step2.vcf")
                cmd4 = [bcftools,"norm","-m-both","-o",preanno_step_1,input_file]
                cmd5 = [bcftools,"norm","-f",reference_genome,"-o",preanno_step_2,preanno_step_1]

                # process4 = subprocess.Popen(cmd4).wait()
                # process5 = subprocess.Popen(cmd5)
                # process5.communicate()

            if f.endswith("_Left_norm.ano.input.step2.vcf"):
                input_file = main_dir + "/" + f
                anno_output_file = main_dir + "/" + f.replace("_Left_norm.ano.input.step2.vcf","_anno")
                cmd6 = [table_annovar,input_file,humandb,"--buildver","hg38","--outfile",anno_output_file,"--remove","--protocol","refGene,ensGene,cytoBand,clinvar_20190305,dbnsfp35a,exac03,gnomad211_exome,kaviar_20150923,avsnp150",
                        "--operation","g,g,r,f,f,f,f,f,f","--nastring",".","--polish", "--vcfinput"]

                process6 = subprocess.Popen(cmd6)
                process6.communicate()

except Exception as e:
    print("An Error has occurred while running the method annovar_annotation. The error is :  {}".format(e))
