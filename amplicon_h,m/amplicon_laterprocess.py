#!/home/grlee/anaconda3/bin/python
#!/bin/bash
#!/home/grlee/anaconda3/envs/qiime2-2022.2/bin/qiime


############# amplicon data 분석과정 전체 진행 코드 정리 #############
###################(abundance table 이후 과정)###################
##############################################################
### ngs data pre-processing(percent filtering 필요 없는 경우) ###
#############################################################


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
parser.add_argument('-p', '--percentage', help = 'cut percent', type = float, default = 0)
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
-p otu clustering 할 퍼센티지 (otu clustering 안할거면 0입력, 할거면 0.99같은 identity 퍼센티지 입력하기-앞선 확인 과정으로 99퍼센트로 실행해왔었는데 추후 상황에 따라 identity threshold 바뀔 수 있음)
'''


########################## asv abundance table 파일(.tsv) format convert 후에 qiime에 다시 import 시키기 ########################### 
print()
print("--------------import qiime asv absolute table--------------")

changed_asv_t = o + '/' + n + '/3-dada2/'+ n + '_abundancetable_format/asv-abundance-table.tsv'
out_path = o + '/' + n + '/3-dada2/'+ n + '_abundancetable_format'

#file format convert 하기 전에 우선 qiime에 import 할 수 있는 형태로 데이터프레임 형태 변경해야함
df = pd.read_csv(changed_asv_t , sep ='\t')
df = df.set_index(['Unnamed: 0'])
df = df.transpose()
final_asv_t = out_path + '/asv-table-absolute.tsv'
df.to_csv(final_asv_t, sep = '\t')


#biom convert
out_biom = out_path + '/asv_percent.biom'
cmd =  "/home/grlee/anaconda3/envs/qiime2-2022.2/bin/biom convert \
-i {} \
-o {} \
--to-hdf5".format(final_asv_t, out_biom)
os.system(cmd)

#biom to qza
out_qza = out_path + '/asv-table-absolute.qza'
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


#p값을 float로 받기 때문에 이를 파일명에 쓰기 위해서 str로 전환해주는 과정 필요
p2 = str(p)


##################### otu clustering directory 형성 (p값을 입력했을때만 -> 0이 아닐때)########################
if p != 0:
    otu_path = o + '/' + n + '/' + 'otu' + p2
    os.mkdir(otu_path)
    
    dir_list = [ o + '/' + n + '/' + 'otu' + p2 + '/3-dada2', o + '/' + n + '/' + 'otu' + p2 + '/4-diversity', o + '/' + n + '/' + 'otu' + p2 + '/5-phylogeny', o + '/' + n + '/' + 'otu' + p2 + '/6-taxonomy']
    for path in dir_list:
        os.mkdir(path)

    
    diversity = [ o + '/' + n + '/' + 'otu' + p2 +'/4-diversity/alpha', o + '/' + n + '/' + 'otu' + p2 +'/4-diversity/beta',] 
    for path in diversity:
        os.mkdir(path)             
 


    ############################# otu clustering ##############################
    print()
    print("--------------filtering sequence for otu clustering--------------")

    dada2_out_r = o + '/' + n + '/3-dada2/' + n + '_rep-seq.qza'
    dada2_out_t = o + '/' + n + '/3-dada2/' +  n + '_table.qza'


    ######## otu clustering 시작 ########
    print()
    print("--------------start otu clustering--------------")

    #usearch 프로그램 사용해 clustering하는 방식은 sequence 기반 clustering 이라 clustering이 반영된 table 정보가  차후 필요하기 때문에 qiime 내 vsearch 기능 사용함

    cluster_out_t = o + '/' + n + '/otu' + p2 + '/' + '3-dada2/' +  n + '_otucluster_table' + '.qza'
    cluster_out_r = o + '/' + n + '/otu' + p2 + '/' + '3-dada2/' +  n + '_otucluster_seq' + '.qza'


    cmd =  "/home/grlee/anaconda3/envs/qiime2-2022.2/bin/qiime vsearch cluster-features-de-novo \
    --i-table {} \
    --i-sequences {} \
    --p-perc-identity {} \
    --o-clustered-table {} \
    --o-clustered-sequences {}".format(dada2_out_t, dada2_out_r, p, cluster_out_t, cluster_out_r)
    os.system(cmd)


    print()
    print("--------------finish--------------")



    ######### otu clustering abundance table 만들기 #########
    print()
    print("--------------abundance table for otu clustering--------------")

    qza_out_2 = o + '/' + n + '/otu' + p2 + '/' + '3-dada2/' +  n + '_otu_abundancetable' + '.qza'
    cmd =  "/home/grlee/anaconda3/envs/qiime2-2022.2/bin/qiime feature-table relative-frequency \
    --i-table {} \
    --o-relative-frequency-table {}".format(cluster_out_t, qza_out_2)
    os.system(cmd)

    
    print()
    print("--------------finish--------------")
    
    
    
    ######### otu clustering abundance table > to tsv format ######### 
    print()
    print("--------------abundance table to tsv format--------------")

    #biom convert
    out_path = o + '/' + n + '/otu0.99/3-dada2/' + n + '_abundancetable_format'
    cmd =  "/home/grlee/anaconda3/envs/qiime2-2022.2/bin/qiime tools export \
    --input-path {} \
    --output-path {}".format(qza_out_2, out_path)
    os.system(cmd)


    #biom to tsv
    out_tsv = out_path + '/feature-table.tsv'
    cmd =  "/home/grlee/anaconda3/envs/qiime2-2022.2/bin/biom convert \
    -i {} \
    -o {} \
    --to-tsv".format(out_path + '/feature-table.biom', out_tsv)
    os.system(cmd)


    #tsv 파일 새로 저장하기(합이 1이 되도록 absolute 한 수치로 변환하는 과정 필요)
    df = pd.read_csv(out_tsv, sep='\t', header = 1)
    df = df.set_index("#OTU ID")
    df = df.div(df.sum(axis=1), axis = 0)
    df.to_csv(out_path + '/otu-abundance-table.tsv' , sep = '\t')


    #파일 변환 과정에서 생긴 중간파일 삭제
    os.remove(out_path + '/feature-table.tsv')
    os.remove(out_path + '/feature-table.biom')

    print()
    print("--------------finish--------------")
    
    
    #####다시 qza로 변환 (taxonomy abundance table형성시 필요)#####
    #biom convert
    out_biom = out_path + '/otu-abundance-table.biom'
    final_otu_t = out_path + '/otu-abundance-table.tsv'
    
    cmd =  "/home/grlee/anaconda3/envs/qiime2-2022.2/bin/biom convert \
    -i {} \
    -o {} \
    --to-hdf5".format(final_otu_t, out_biom)
    os.system(cmd)

    #biom to qza
    out_qza = out_path + '/otu-abundance-table.qza'
    cmd = "/home/grlee/anaconda3/envs/qiime2-2022.2/bin/qiime tools import \
    --input-path {} \
    --type 'FeatureTable[Frequency]' \
    --input-format 'BIOMV210Format' \
    --output-path {}".format(out_biom, out_qza)
    os.system(cmd)

    #변환 과정에서 생기는 biom 파일 삭제 
    os.remove(out_biom)



########################## taxonomy assignment ##########################
print()
print("--------------start taxonomy assignment--------------")
taxa_out = o + '/' + n + '/6-taxonomy/asv' + '_taxonomy_assignment95' + '.qza'
taxa_out_2 = o + '/' + n + '/otu' + p2 + '/' + '6-taxonomy/otu' + p2 + '_taxonomy_assignment95' + '.qza'
dada2_out_r = o + '/' + n + '/3-dada2/' + n + '_rep-seq.qza'

if p!= 0:
    cmd = "/home/grlee/anaconda3/envs/qiime2-2022.2/bin/qiime feature-classifier classify-consensus-blast \
    --i-query {} \
    --i-reference-reads '/home/grlee/db/NCBI-classify-consensus/ncbi_16S_seq.qza' \
    --i-reference-taxonomy '/home/grlee/db/NCBI-classify-consensus/ncbi_16S_tax.qza' \
    --p-perc-identity 0.95 \
    --o-classification {}".format(cluster_out_r, taxa_out_2)
    os.system(cmd)
    print()
    print("--------------finish--------------")
    
    print()
    print("--------------visualization taxonomy assignment--------------")
    vis_out_2 = o + '/' + n + '/otu' + p2 + '/' + '6-taxonomy' + '/otu' + p2 + '_taxonomy_assignment95' + '.qzv'
    cmd = "/home/grlee/anaconda3/envs/qiime2-2022.2/bin/qiime metadata tabulate \
    --m-input-file {} \
    --o-visualization {}".format(taxa_out_2, vis_out_2)
    os.system(cmd)
    print()
    print("--------------finish--------------")
        



cmd = "/home/grlee/anaconda3/envs/qiime2-2022.2/bin/qiime feature-classifier classify-consensus-blast \
--i-query {} \
--i-reference-reads '/home/grlee/db/NCBI-classify-consensus/ncbi_16S_seq.qza' \
--i-reference-taxonomy '/home/grlee/db/NCBI-classify-consensus/ncbi_16S_tax.qza' \
--p-perc-identity 0.95 \
--o-classification {}".format(dada2_out_r, taxa_out)
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
if p !=  0:
    print()
    print("--------------taxonomy abundance table start--------------")

    taxa_abundancet_2 = o + '/' + n + '/otu' + p2 + '/' + '6-taxonomy/otu' + p2 + '_taxonomy_assignment95_abundancetable' + '.qza'
    taxa_collapsed_2 =  o + '/' + n + '/otu' + p2 + '/' + '6-taxonomy/otu' + p2 + '_taxonomy_assignment95_collapsed' + '.qza'

    cmd = "/home/grlee/anaconda3/envs/qiime2-2022.2/bin/qiime taxa collapse \
    --i-table {} \
    --i-taxonomy {} \
    --p-level 7 \
    --o-collapsed-table {}".format(out_qza, taxa_out_2, taxa_collapsed_2)
    os.system(cmd)
    
    cmd = "/home/grlee/anaconda3/envs/qiime2-2022.2/bin/qiime feature-table relative-frequency \
    --i-table {} \
    --o-relative-frequency-table {}".format(taxa_collapsed_2, taxa_abundancet_2)
    os.system(cmd)

    print()
    print("--------------finish--------------")
    
    
    
    #biom convert
    out_path = o + '/' + n + '/otu' + p2 + '/' + '6-taxonomy/otu' + p2 + '_taxonomy_assignment95_abundancetable_format'
    cmd =  "/home/grlee/anaconda3/envs/qiime2-2022.2/bin/qiime tools export \
    --input-path {} \
    --output-path {}".format(taxa_abundancet_2, out_path)
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





##################################### asv cut percent table except zero ( 필요에따라 사용하면됨 ) ######################################
#표에서 전체 값이 0인 값들로만 이뤄진 컬럼 삭제하기 위함
print()
print("-------------- table filtering except zero -------------")

out_path = o + '/' + n + '/6-taxonomy'
changed_asv_filtered = out_path + '/asv_percent_except_zero.tsv'
df = pd.read_csv(changed_asv_t, sep = '\t')
for i in df.columns:
    if df[i].sum() == 0 :                    #컬럼별 합이 0이 나오는거
        df.drop(columns = i, inplace = True)

df.rename(columns={'Unnamed: 0' : 'id'}, inplace = True)            
df.to_csv(changed_asv_filtered , sep = '\t')
print()
print("--------------finish--------------")







###################################### taxonomy barplot with matplotlib ######################################
if p != 0 :
    print()
    print("--------------taxonomy barplot with matplotlib-------------")

    out_path = o + '/' + n + '/otu' + p2 + '/' + '6-taxonomy'

    taxa_abundancet = out_path + '/otu0.99_taxonomy_assignment95_abundancetable_format/taxonomy_abundancetable.tsv'
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
if p != 0 :
    tax = []
    df = pd.read_csv(o + '/' + n + '/otu' + p2 + '/' + '6-taxonomy/otu0.99_taxonomy_assignment95_abundancetable_format/taxonomy_abundancetable.tsv', sep = '\t', header = 1 )

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
    df.to_csv(o + '/' + n + '/otu' + p2 + '/' + '6-taxonomy/taxonomybarplot_table.tsv', sep = '\t')
    
    
     
    
    
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

