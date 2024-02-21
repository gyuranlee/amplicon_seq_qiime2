#!/home/grlee/anaconda3/bin/python
#!/bin/bash
#!/home/grlee/anaconda3/envs/qiime2-2022.2/bin/qiime


### amplicon sequencing data qiime처리 process ### (""" + binary 한 수치로 변환하는 과정까지 """ ) ####
###############################################################################################
################################### ngs data pre-processing ###################################
###############################################################################################


# [amplicon sequencing data를 qiime2(현재 스크립트에선 2022.2 버전이용)를 이용해 abundance table 만드는 과정까지 정리한 python script file]
# 작성자 - 이규란
# 최초 작성 일 - 2022.11.15
# version - 1.0


#########################################입력 변수(사용법)##########################################
#모든 경로는 /home/grlee/~ 부터 입력 해야함

'''
현재 입력 해야할 변수들
i = '/home/grlee/raw/mixed~/ngsdata/H~/'(fastq파일저장된경로)
m = '/home/grlee/raw/mixed~/mixedcolony_metadata_1.tsv' (metadata 파일 저장경로)
n = mix_NGS_1(파일명 - 현재 날짜로 지정하고 있음)   >> 저장될 output경로에 n이라는 폴더를 만들기 위함 
o = '/home/grlee/output/~ 저장할 경로(마지막에 슬래쉬 붙이지 않기 !!!주의!!!)
f = '0' 몇개의 feature를 가지는 sample까지 filtering 할건지 ( 샘플 개수가 너무 적다면 filtering하지 않거나(-f 0입력하기) 샘플 개수를 고려해서 filtering 할 개수를 정해야함. mixed의 경우 우선 다 0으로 했음  )
fl = 200
rl = 140
'''
###############################################################################################




import os, glob, argparse, natsort, sys
import pandas as pd
import numpy as np

parser = argparse.ArgumentParser(description = 'mixed_colony')
parser.add_argument('-i','--input', help = 'input path for fastq file', type = str)
parser.add_argument('-o','--output', help = 'output path', type = str)
parser.add_argument('-p', '--percentage', help = 'cut percent', type = str)
parser.add_argument('-type', help = 'qiime import type', type = str, default = 'SampleData[PairedEndSequencesWithQuality]')
parser.add_argument('-if', '--input-format', help = 'qiime input format', default = 'PairedEndFastqManifestPhred33V2')
parser.add_argument('-p-front-f','--front', help = 'qiime cutadapt primer trimmed front f', type = str)
parser.add_argument('-p-front-r', '--reverse', help = 'qiime cutadapt primer trimmed front r', type = str)
parser.add_argument('-m-', '--metadata', help = 'metadata file', type = str)
parser.add_argument('-d', '--directory', help = 'make input directory file folder', type = str) 
parser.add_argument('-n', '--name', help = 'make file name with date', type = str)
parser.add_argument('-f', '--filtering', help = 'how many sample filtering when feature table filter-features', type = int, default = 0)
parser.add_argument('-fl', '--forwardlength', help = 'trimming length for front end in dada2', type = int, default = 200)
parser.add_argument('-rl', '--reverselength', help = 'trimming length for reverse end in dada2', type = int, default = 140)
args = parser.parse_args()

###########################################################################

# 변수지정
i = args.input             #input raw data 저장경로(fastq파일)
o = args.output            #output 파일 저장할 경로
m = args.metadata          #metadata 파일 저장경로 (tsv)
p = args.percentage        #otu clustering percentages
t = args.type
n = args.name              #output file, directory name
fr = args.front
rv = args.reverse
d = args.directory
f = args.filtering
fl = args.forwardlength
rl = args.reverselength

print("#################################### argument ######################################")
print("##################  input fastq raw data path: {}".format(i))
print("##################  output directory path: {}".format(o + '/' + n))
print("##################  output directory name: {}".format(n))
print("####################################################################################")



# 하위 폴더 생성하기

###################################################################################################################
# output 디렉토리에 이미 name으로 된 하위 디렉토리가 존재할 시, 에러 코드 출력

out_path = o + '/' + n
if os.path.exists(out_path) == True:
    print()
    print()
    sys.exit("******* FileExistsError: 이미 해당 name으로 생성된 디렉토리가 존재합니다 *******")

else:
    os.mkdir(out_path)
 
###################################################################################################################
## directory file folder input 날짜에 맞춰서 만들기(하위 파일 여러개 만들때 리스트 사용) ##

dir_list = [ o + '/' + n +'/0-raw', o + '/' + n +'/1-import', o + '/' + n +'/2-cutadapt', o + '/' + n +'/3-dada2', o + '/' + n +'/4-diversity', o + '/' + n +'/5-phylogeny',  o + '/' + n +'/6-taxonomy']
for path in dir_list:
    os.mkdir(path)

    
diversity = [ o + '/' + n +'/4-diversity/alpha', o + '/' + n +'/4-diversity/beta',] 
for path in diversity:
    os.mkdir(path)             


###################################################################################################################
##############################################fastq 파일 to manifest################################################
# index = False 쓰면 저장될때 숫자열 인덱스 사라지고 저장됨. # 파일 포맷은 tsv로 지정해야함

print()
print("---------------------fastq to manifest---------------------")


lst = []
# input data file (fastq.gz) 의 경로를 dir로 지정해줘야함
#처음 input으로 fastqc 데이터 들어있는 디렉토리 지정해주면되고 / output으로 앞으로 생길 파일명에 들어갈 고유 표시(ex.날짜)입력하기
target_dir = i


# forward, reverse 파일 샘플명이 같아서 둘중 하나만 불러와서 샘플명 따서 첫 컬럼에 넣기위해
for i in os.listdir(target_dir):
    if not ('.ipynb_checkpoints' in i) and ('_1.fastq.gz' in i):    # 스크립트를 로데이터 있는 경로에서 만들면 생기는 파일임 만약 있다면 확인해서 제거하기
        lst.append(i.split('_1.fastq.gz')[0])
        lst = natsort.natsorted(lst)

# forward 파일 경로만 불러오기 위해서 조건 생성   
lst2 = []    
for i in glob.glob(os.path.join(target_dir, '*1.fastq.gz')):
    lst2.append(i)
    lst2 = natsort.natsorted(lst2)

# reverse 파일 경로만 불러오기 위해서 조건 생성    
lst3 = []
for i in glob.glob(os.path.join(target_dir, '*2.fastq.gz')):
    lst3.append(i)
    lst3 = natsort.natsorted(lst3)
    

df = pd.DataFrame([lst, lst2, lst3])
# 형식에 맞추기 위해 행, 열 전환해 줘야함 / 컬럼명 지정하기
df = df.transpose()
df.columns = ['sample-id', 'forward-absolute-filepath', 'reverse-absolute-filepath']
manifest_path = o + '/' + n + '/1-import/' + n + '_mani.tsv'
df
df.to_csv(manifest_path, sep='\t', index = False)  #파일 저장시에 쓰기
# index = False 쓰면 저장될때 숫자열 인덱스 사라지고 저장됨. # 파일 포맷은 tsv로 지정해야함

print()
print("-----------------fastq to manifest finish------------------")






###################### 여기서 부터 qiime 실행이 되어야함 #########################
###########################################################################
###########################################################################




#################################################qiime import#######################################################
# manifest로 나온 output이 import가 되어야함
print()
print("----------------------qiime import-------------------------")
            


qiime_out = o + '/' + n + '/1-import/' + n + '_mani.qza'


cmd =  "/home/grlee/anaconda3/envs/qiime2-2022.2/bin/qiime tools import \
--type 'SampleData[PairedEndSequencesWithQuality]' \
--input-path {} \
--output-path {} \
--input-format PairedEndFastqManifestPhred33V2".format(manifest_path, qiime_out)
os.system(cmd)

print()
print("--------------------qiime import finish--------------------")


########################################### visualization (quality check) ##########################################
# import가 제대로 되었는지, quality plot보고 퀄리티 형태 확인하기

print()
print("-------------------import file visualizaion----------------")

mani_visual = o + '/' + n + '/1-import/' + n + '_mani.qzv'
cmd =  "/home/grlee/anaconda3/envs/qiime2-2022.2/bin/qiime demux summarize \
--i-data {} \
--o-visualization {}".format(qiime_out, mani_visual)
os.system(cmd)

print()
print("---------------import file visualizaion finish-------------")


################################################### cutadapt #######################################################
print()
print("-----------------cutadapt primer trimming------------------")

# py에서는 입력받을 변수로 지정하기
pf = 'GTGCCAGCMGCCGCGGTAA'
pr = 'GGACTACHVGGGTWTCTAAT'

cutadapt_out = o + '/' + n + '/2-cutadapt/' + n + '_primer_trimmed.qza'
cmd =  "/home/grlee/anaconda3/envs/qiime2-2022.2/bin/qiime cutadapt trim-paired \
--i-demultiplexed-sequences {} \
--p-front-f {} \
--p-front-r {} \
--p-discard-untrimmed \
--o-trimmed-sequences {} \
--p-cores 32".format(qiime_out, pf, pr, cutadapt_out) 
os.system(cmd)
 # primer 없던 sequence 는 버리기 (discard-untrimmed)
print()
print("---------------------cutadapt finish-----------------------")


############################################# cutadapt visualization ###############################################
print()
print("-----------------primmer trimmed visualizaion--------------")

cutadapt_visual = o + '/' + n + '/2-cutadapt/' + n + '_primer_trimmed.qzv'
cmd =  "/home/grlee/anaconda3/envs/qiime2-2022.2/bin/qiime demux summarize \
--i-data {} \
--o-visualization {}".format(cutadapt_out, cutadapt_visual)
os.system(cmd)

print()
print("-------------primmer trimmed visualizaion finish-----------")




##################################################### dada2 ########################################################
print()
print("-------------------denoising with dada2--------------------")

# denoising 과정에서 merge 시키고 다듬는 과정 진행됨 

dada2_out_t = o + '/' + n + '/3-dada2/' + n + '_table.qza'
dada2_out_s = o + '/' + n + '/3-dada2/' + n + '_stats.qza'
dada2_out_r = o + '/' + n + '/3-dada2/' + n + '_rep-seq.qza'

cmd =  "/home/grlee/anaconda3/envs/qiime2-2022.2/bin/qiime dada2 denoise-paired \
--i-demultiplexed-seqs {} \
--p-trunc-len-f {} \
--p-trunc-len-r {} \
--o-table {} \
--o-representative-sequences {} \
--o-denoising-stats {}".format(cutadapt_out, fl, rl, dada2_out_t, dada2_out_r, dada2_out_s)

os.system(cmd)
print()
print("------------------------dada2 finish-----------------------")


################################################ dada2 visualization ###############################################
print()
print("-------------------visualization dada2---------------------")

# 이 과정부터 metadata file 필요함
# table, seq, stat 각각 visualization 하는 명령어 다르니 필요한 목적에 따라 선택하면 됨
# chimera percent , merge percent 정보가 필요하면 stat // feature count나 asv정보가 필요하면 table을 봐야함

dada2_visual_t = o + '/' + n + '/3-dada2/' + n + '_table.qzv'
dada2_visual_s = o + '/' + n + '/3-dada2/' + n + '_stats.qzv'
dada2_visual_r = o + '/' + n + '/3-dada2/' + n + '_rep-seq.qzv'


#table
metadata = m
cmd =  "/home/grlee/anaconda3/envs/qiime2-2022.2/bin/qiime feature-table summarize \
--i-table {} \
--m-sample-metadata-file {} \
--o-visualization {}".format(dada2_out_t, metadata , dada2_visual_t)
os.system(cmd)


#visualization 단계는 table, stat, seq 필요에 따라 주석 해제하기
#stat
cmd =  "/home/grlee/anaconda3/envs/qiime2-2022.2/bin/qiime metadata tabulate \
--m-input-file {} \
--o-visualization {}".format(dada2_out_s, dada2_visual_s)

# #seq
# cmd =  "/home/grlee/anaconda3/envs/qiime2-2022.2/bin/qiime feature-table tabulate-seqs \
# --i-data {} \
# --o-visualization {}".format(dada2_out_r, dada2_visual_r)

os.system(cmd)
print()
print("--------------------------finish---------------------------")




#이 과정을 f값에 따라서 넣을지 말지를 선택하게끔 (왜냐면 f = 0 인채로 돌면 오류가 남)

if f != 0 :
    
################################################# filter table (샘플개수에 따라 다르게 해야함) #####################################################
    print()
    print("-----------filtering table feature count under f----------")

# feature count 적은 asv의 경우 기술적인 오류로 생겨난 경우일 가능성이 높아서 filtering 하는 단계

    filtering_out =  o + '/' + n + '/3-dada2/' + n + '_filtered_table_' + f + '.qza'



    cmd =  "/home/grlee/anaconda3/envs/qiime2-2022.2/bin/qiime feature-table filter-features \
    --i-table {} \
    --p-min-samples f \
    --o-filtered-table {}".format(dada2_out_t, filtering_out)

    os.system(cmd)
    print()
    print("-------------------filtering table finish------------------")
    




###################################### abundance table ##########################################
print()
print("-----------------abundance table each asvs-----------------")

if f != 0:
    filtering_out =  o + '/' + n + '/3-dada2/' + n + '_filtered_table_' + f + '.qza'
else:
    filtering_out = dada2_out_t


abundancetable_out = o + '/' + n + '/3-dada2/' + n + '_abundancetable.qza'


cmd =  "/home/grlee/anaconda3/envs/qiime2-2022.2/bin/qiime feature-table relative-frequency \
--i-table {} \
--o-relative-frequency-table {}".format(filtering_out, abundancetable_out)
os.system(cmd)
print()
print("-------------------abundance table finish------------------")


########################################## abundancee table visualization ##########################################
print()
print("---------------abundance table visualization---------------")

abundancetable_vis = o + '/' + n + '/3-dada2/' + n + '_abundance_table.qzv'


cmd =  "/home/grlee/anaconda3/envs/qiime2-2022.2/bin/qiime metadata tabulate \
--m-input-file {} \
--o-visualization {}".format(abundancetable_out ,abundancetable_vis)
os.system(cmd)

print()
print("--------------------------finish---------------------------")





########################################## abundancee table to tsv format ##########################################
print()
print("--------------abundance table to tsv format--------------")

#biom convert
out_path = o + '/' + n + '/3-dada2/' + n + '_abundancetable_format'
cmd =  "/home/grlee/anaconda3/envs/qiime2-2022.2/bin/qiime tools export \
--input-path {} \
--output-path {}".format(abundancetable_out, out_path)
os.system(cmd)


#biom to tsv
out_tsv = out_path + '/feature-table.tsv'
cmd =  "/home/grlee/anaconda3/envs/qiime2-2022.2/bin/biom convert \
-i {} \
-o {} \
--to-tsv".format(out_path + '/feature-table.biom', out_tsv)
os.system(cmd)


##!!!여기서 나온 tsv파일이 tap 이 안된 상태로 나와짐..그걸 고치기 위해서 tsv파일 읽어와서 tap 넣어서 수정하는 과정 한번 더 거침 
#tsv 파일 새로 저장하기
df = pd.read_csv(out_tsv, sep='\t', header = 1)
df = df.transpose()
df = df.rename(columns = df.iloc[0])
df = df.drop(df.index[0])
df.to_csv(out_path + '/asv-abundance-table.tsv' , sep = '\t')


#파일 변환 과정에서 생긴 중간파일 삭제
os.remove(out_path + '/feature-table.tsv')
os.remove(out_path + '/feature-table.biom')

print()
print("--------------finish--------------")







# ######################################### abundancee table to binary value (효림님 드린 정보)##########################################
# print()
# print("--------------abundance table with binary value--------------")
# df = pd.read_csv(out_path + '/feature-table-tap.tsv' , sep = '\t')

# pd.set_option('display.max_rows', None)
# pd.set_option('display.max_columns', None)

# #필요에 따라 데이터프레임 형식 바꿔줘야함 (feature가 컬럼으로 와야함 / 첫째 행 컬럼으로 지정해야하고 첫 컬럼 네임 'id' 되도록 변경하기)
# df = df.transpose()                           # 행,열 전환
# df = df.rename(columns = df.iloc[0])          # 첫번째 행 정보 컬럼으로 만들기
# df = df.drop(df.index[0])                     # 첫번째 행 남아있는거 삭제
# df = df.rename(columns = {'#OTU ID' : 'id'})  # 첫 컬럼 이름 id로 바꿔주기(만약 df에 otu id가 컬럼 네임이 아니라면 저 부분 !!!!!변경해야함!!!!!)

# id_list = df['id'].to_list()
# df.index = id_list
# df.drop(['id'], axis = 1, inplace = True)     # id 내용들이 인덱스로 오게끔 바꿔주기


# #함수설정
# def f(x):
#     if x != '0.0':                            # 이 데이터프레임은 0이 아닌 0.0이었음 이거 확인하고 다르면 !!!!!변경해줘야함!!!!!
#         return 1
#     else:
#         return 0


# pd.set_option('display.max_rows', None)
# pd.set_option('display.max_columns', None)

# # df

# df = df.applymap(f)
# df.to_csv(out_path + '/feature-table-binary.tsv' , sep = '\t')

# print()
# print("--------------finish--------------")


