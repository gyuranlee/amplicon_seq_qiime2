#!/home/grlee/anaconda3/bin/python
#!/bin/bash
#!/home/grlee/anaconda3/envs/qiime2-2022.2/bin/qiime


####### [[ mixed_laterprocess.py ]] #######
### mixed colony 분석과정 전체 진행 코드 정리 ###
#########(abundance table 이후 과정)#########
###########################################
### ngs data mixed colony pre-processing ###
############################################


#n은 이후 과정이 이전 과정과 다른 날짜에 진행되어도 이전 날짜 기준으로 맞추는 것으로 함



import os, glob, argparse, natsort,sys
import pandas as pd
import numpy as np
import seaborn as sns
import random
import plotly as plt
import matplotlib.pyplot as plt 

parser = argparse.ArgumentParser(description = 'mixed_colony')
parser.add_argument('-i','--input', help = 'input path for fastq file', type = str)
parser.add_argument('-o','--output', help = 'output path', type = str)
parser.add_argument('-p', '--percentage', help = 'cut percent', type = str)
parser.add_argument('-m-', '--metadata', help = 'metadata file', type = str)
parser.add_argument('-d', '--directory', help = 'make input directory file folder', type = str) 
parser.add_argument('-n', '--name', help = 'make file name with date', type = str)
parser.add_argument('-v', '--value', help = 'asv filtering value', type = str)
args = parser.parse_args()



###########################################################################

# 변수지정
i = args.input             #input raw data 저장경로(fastq파일)
o = args.output            #output 파일 저장할 경로
m = args.metadata          #metadata 파일 저장경로 (tsv)
p = args.percentage        #otu clustering percentages
v = args.value             #asv filtering value
n = args.name              #output file, directory name
d = args.directory

###########################################################################


'''
현재 입력 변수

-o /home/grlee/output/분석할 데이터 디렉토리(-ex.mixed)(마지막에 슬래쉬 붙이지 않기 !!!주의!!!)
-n 날짜 디렉토리명에 붙어있는 날짜(-ex.mix_NGS_1)
-v asv filtering 할 수치 (그래프상 확인하고 정하기)
'''




###########################cut percent 정하는 단계#############################

########################## 퍼센트 정한 후에 수정된 파일로 저장하기 ########################### 
########################## 0.01보다 작은 asv 값 0으로 교체하기 ############################    
#이 퍼센트값 파라미터로 둘지 고민해보기 (그럼 디폴트값 줘야할듯)
print()
print("--------------change table with cut asv percent--------------")

out_path = o + '/' + n + '/3-dada2/' + n + '_abundancetable_format' 
df = pd.read_csv(out_path + '/asv-abundance-table.tsv', sep = '\t')
df = df.set_index("Unnamed: 0")

for i in df:
    df.loc[df[i] < float(v), i] = 0.000000
    
for i in df.columns:
    if df[i].sum() == 0 :                    #컬럼별 합이 0이 나오는거
        df.drop(columns = i, inplace = True)


changed_asv_t = out_path + '/asv_percent' + v + '.tsv'
df.to_csv(changed_asv_t,sep = '\t')


print()
print("--------------finish--------------")



########################## 파일 format convert 후에 qiime에 다시 import 시키기 ########################### 
print()
print("--------------import qiime asv percent cut table--------------")

#file format convert 하기 전에 우선 qiime에 import 할 수 있는 형태로 데이터프레임 형태 변경해야함
df = pd.read_csv(changed_asv_t , sep ='\t')
df = df.set_index(['Unnamed: 0'])
df = df.transpose()
final_asv_t = out_path + '/asv_percent' + v + '.tsv'
df.to_csv(final_asv_t, sep = '\t')


#biom convert
out_biom = out_path + '/asv_percent' + v + '.biom'
cmd =  "/home/grlee/anaconda3/envs/qiime2-2022.2/bin/biom convert \
-i {} \
-o {} \
--to-hdf5".format(final_asv_t, out_biom)
os.system(cmd)

#biom to qza
out_qza = out_path + '/asv_percent' + v + '.qza'
cmd = "/home/grlee/anaconda3/envs/qiime2-2022.2/bin/qiime tools import \
--input-path {} \
--type 'FeatureTable[Frequency]' \
--input-format 'BIOMV210Format' \
--output-path {}".format(out_biom, out_qza)
os.system(cmd)

#변환 과정에서 생기는 biom 파일 삭제 
os.remove(out_biom)

print()
print("--------------finish--------------")




########################## table visualization ########################### 
print()
print("--------------visualization asv percent filtered table--------------")

asv_vis = out_path + '/asv_percent' + v + '.qzv'
cmd = "/home/grlee/anaconda3/envs/qiime2-2022.2/bin/qiime feature-table summarize \
--i-table {} \
--o-visualization {}".format(out_qza, asv_vis)
os.system(cmd)

print()
print("--------------finish--------------")



############################################ filter seq (for tax assign) #################################################
### table을 absolute한 수치로 만들고 > 그 후에 퍼센트 설정해서 filtering한 테이블을 이용해 seq도 filtering해줌 (taxonomy assign과정에 필요)
print()
print("-----------filtering sequence using asv percent filtered table----------")

filtering_outseq = out_path + '/asv_percent' + v + '_seq.qza'
dada2_out_r = o + '/' + n + '/3-dada2/' + n + '_rep-seq.qza'
out_qza = out_path + '/asv_percent' + v + '.qza'

cmd =  "/home/grlee/anaconda3/envs/qiime2-2022.2/bin/qiime feature-table filter-seqs \
--i-data {} \
--i-table {} \
--o-filtered-data {}".format(dada2_out_r, out_qza, filtering_outseq)
os.system(cmd)
print()
print("-------------------filtering seq finish------------------")



########################## taxonomy assignment ##########################
print()
print("--------------start taxonomy assignment--------------")
taxa_out = o + '/' + n + '/6-taxonomy/asv' + '_taxonomy_assignment95' + '.qza'

cmd = "/home/grlee/anaconda3/envs/qiime2-2022.2/bin/qiime feature-classifier classify-consensus-blast \
--i-query {} \
--i-reference-reads '/home/grlee/db/NCBI-classify-consensus/ncbi_16S_seq.qza' \
--i-reference-taxonomy '/home/grlee/db/NCBI-classify-consensus/ncbi_16S_tax.qza' \
--p-perc-identity 0.95 \
--o-classification {}".format(filtering_outseq, taxa_out)
os.system(cmd)
print()
print("--------------finish--------------")


print()
print("--------------visualization taxonomy assignment--------------")
vis_out = o + '/' + n + '/6-taxonomy/asv' + '_taxonomy_assignment95' + '.qzv'
cmd = "/home/grlee/anaconda3/envs/qiime2-2022.2/bin/qiime metadata tabulate \
--m-input-file {} \
--o-visualization {}".format(taxa_out, vis_out)
os.system(cmd)
print()
print("--------------finish--------------")



########################## taxonomy abundance table to taxa barplot ##########################
#차임에서 barplot그리는데 이용된 table을 살펴보면 abundance table형태로 바꿔줘야함 
print()
print("--------------taxonomy abundance table start--------------")

taxa_abundancet = o + '/' + n + '/6-taxonomy/asv' + '_taxonomy_assignment95_abundancetable' + '.qza'
taxa_collapsed =  o + '/' + n + '/6-taxonomy/asv' + '_taxonomy_assignment95_collapsed' + '.qza'

cmd = "/home/grlee/anaconda3/envs/qiime2-2022.2/bin/qiime taxa collapse \
--i-table {} \
--i-taxonomy {} \
--p-level 7 \
--o-collapsed-table {}".format(out_qza, taxa_out, taxa_collapsed)
os.system(cmd)


cmd = "/home/grlee/anaconda3/envs/qiime2-2022.2/bin/qiime feature-table relative-frequency \
--i-table {} \
--o-relative-frequency-table {}".format(taxa_collapsed, taxa_abundancet)
os.system(cmd)

print()
print("--------------finish--------------")



print()
print("--------------abundance table to tsv format--------------")

#biom convert
out_path = o + '/' + n + '/6-taxonomy/asv' + '_taxonomy_assignment95_abundancetable_format'
cmd =  "/home/grlee/anaconda3/envs/qiime2-2022.2/bin/qiime tools export \
--input-path {} \
--output-path {}".format(taxa_abundancet, out_path)
os.system(cmd)



#biom to tsv
out_tsv = out_path + '/feature-table.tsv'
cmd =  "/home/grlee/anaconda3/envs/qiime2-2022.2/bin/biom convert \
-i {} \
-o {} \
--to-tsv".format(out_path + '/feature-table.biom', out_tsv)
os.system(cmd)

#tsv 파일 새로 저장하기
df = pd.read_csv(out_tsv, sep='\t')
# df1 = df.set_index(keys='#OTU ID')
# df1 = df1.transpose()  # 행,열 전환 
df.to_csv(out_path + '/taxonomy_abundancetable.tsv' , sep = '\t')

#파일 변환 과정에서 생긴 중간파일 삭제
os.remove(out_path + '/feature-table.tsv')
os.remove(out_path + '/feature-table.biom')

print()
print("--------------finish--------------")










# ##################################### asv cut percent table except zero ( 필요에따라 사용하면됨 ) ######################################
# #표에서 전체 값이 0인 값들로만 이뤄진 컬럼 삭제하기 위함
# print()
# print("-------------- table filtering except zero -------------")

# out_path = o + '/' + n + '/6-taxonomy'
# changed_asv_filtered = out_path + '/asv_percent_except_zero.tsv'
# df = pd.read_csv(changed_asv_t, sep = '\t')
# for i in df.columns:
#     if df[i].sum() == 0 :                    #컬럼별 합이 0이 나오는거
#         df.drop(columns = i, inplace = True)

# df.rename(columns={'Unnamed: 0' : 'id'}, inplace = True)            
# df.to_csv(changed_asv_filtered , sep = '\t')
# print()
# print("--------------finish--------------")







###################################### taxonomy barplot with matplotlib ######################################
print()
print("--------------taxonomy barplot with matplotlib-------------")

out_path = o + '/' + n + '/6-taxonomy/'

taxa_abundancet = out_path + '/asv_taxonomy_assignment95_abundancetable_format/taxonomy_abundancetable.tsv'
df = pd.read_csv(taxa_abundancet, sep ='\t', header = 1)   # plot 만들어지는 형식에 맞게끔 df 수정하기
df = df.set_index("#OTU ID")
df = df.transpose()
#0인값 표에서 없애기
for i in df.columns:
    if df[i].sum() == 0 :                    #컬럼별 합이 0이 나오는거
        df.drop(columns = i, inplace = True)

df.plot(kind="bar", stacked=True)
plt.rcParams['figure.figsize'] = [30,30]
plt.legend(ncol=2, loc = (1.0, 1.0), bbox_to_anchor = (1, 0.95))
plt.title("Taxonomy barplot")
plt.rc('legend', fontsize = 5)
plt.rc('axes', labelsize = 13)                           #x,y축 label폰트
plt.tick_params(axis = 'x', labelsize =10, width = 25)



# plt.tight_layout()
plt.savefig(out_path + 'barplot_savefig_200dpi.png', dpi=200, bbox_inches ='tight', pad_inches = 0.2)  #dpi는 해상도(디폴트100) #pad_inches 숫자 커질수록 여백 많아짐
print()
print("--------------finish--------------")




#########################################barplot이루는 균주 taxonomy 정보#########################################
tax = []
df = pd.read_csv(o + "/" + n + "/6-taxonomy/asv_taxonomy_assignment95_abundancetable_format/taxonomy_abundancetable.tsv", sep = '\t', header = 1 )

df = df.transpose()
df = df.rename(columns = df.iloc[0])
df = df.drop(df.index[0])
df = df.reset_index()

id = df['index'].values.tolist()


for i in df.columns:
    if df[i].sum() == 0 :                    #컬럼별 합이 0이 나오는거
        df.drop(columns = i, inplace = True)


del df['index']
for i in range (len(df)):
    tax.append(df.columns[df.loc[i] != 0])

df = pd.DataFrame(tax, index = id)
df.to_csv(o + "/" + n + "/6-taxonomy/taxonomybarplot_table.tsv", sep = '\t')
