#!/home/grlee/anaconda3/bin/python
#!/bin/bash
#!/home/grlee/anaconda3/envs/qiime2-2022.2/bin/qiime

####################################################
########### amplicon data diversity 분석 ############
####################################################

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
parser.add_argument('-de', '--depth', help = 'alpha rarefaction sampling depth', type = str)
parser.add_argument('-e', '--vector', help = 'vector for alpha diversity', type = str)
parser.add_argument('-c','--column', help = 'metadata column 선택하기 for beta group significance', type = str)
parser.add_argument('-a', '--distancemetrics', help = 'distance metrics 선택 for alpha group significance', type = str)
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
a = args.distancemetrics
c = args.column
de = args.depth
e = args.vector

###########################################################################


'''
현재 입력 변수

-o /home/grlee/output/분석할 데이터 디렉토리(-ex.mixed)(마지막에 슬래쉬 붙이지 않기 !!!주의!!!)
-n 날짜 디렉토리명에 붙어있는 날짜(-ex.mix_NGS_1)
-de alpha rarefaction에서 max-depth결정 (dada2결과 stat확인 --> non-chimeric수 가장 적은 숫자로 결정)
-m /home/grlee/raw/human/metadatafile~
-e vector 선택하기 (ex.faith_pd_vector)
-a distance matrix 선택하기 core metrix중에서 고르는거 (ex.weighted_unifrac_distance_matrix)
-c metadata column선택하기 (beta과정에서 필요 / categorical 컬럼만 선택 가능하니까 numeric 컬럼을 선택하고 싶다면 메타데이터 수정필요)
'''



################################# alpha rarefaction ###############################
print()
print("--------------alpha rarefaction--------------")

dada2_out_t = o + '/' + n + '/3-dada2/' + n + '_table.qza'
out_rarefacition = o + '/' + n + '/4-diversity/' + de + '_alphararefaction.qzv'


#depth설정시에 table qzv 파일 열어서 확인해봐야함 
cmd = "/home/grlee/anaconda3/envs/qiime2-2022.2/bin/qiime diversity alpha-rarefaction \
--i-table {} \
--p-max-depth {} \
--m-metadata-file {} \
--o-visualization {}".format(dada2_out_t, de, m, out_rarefacition) ## filtering table 사용하는 경우에는 dada2_out_t쓰면 안됨.. 
os.system(cmd)

print()
print("--------------finish--------------")



################################# phylogeny ###############################
print()
print("--------------phylogeny--------------")

dada2_out_r = o + '/' + n + '/3-dada2/' + n + '_rep-seq.qza'
out_alignment = o + '/' + n + '/4-diversity/' + n + '_aligned-seq.qza'
out_maskedalignment = o + '/' + n +  '/4-diversity/' + n + '_masked_aligned-seqs.qza'
out_tree = o + '/' + n + '/4-diversity/' + n + '_unrooted-tree.qza'
out_rootedtree = o + '/' + n + '/4-diversity/' + n + '_rooted-tree.qza'

cmd = "/home/grlee/anaconda3/envs/qiime2-2022.2/bin/qiime phylogeny align-to-tree-mafft-fasttree \
--i-sequences {} \
--o-alignment {} \
--o-masked-alignment {} \
--o-tree {} \
--o-rooted-tree {}".format(dada2_out_r, out_alignment, out_maskedalignment, out_tree, out_rootedtree)
os.system(cmd)

print()
print("--------------finish--------------")



################################# core metrics(diversity분석에 쓰일 파일들) ###############################
print()
print("--------------core metrics--------------")

out_dir = o + '/' + n + '/4-diversity/' + n + '_core-diversity-results'
cmd = "/home/grlee/anaconda3/envs/qiime2-2022.2/bin/qiime diversity core-metrics-phylogenetic \
--i-phylogeny {} \
--i-table {} \
--p-sampling-depth {} \
--m-metadata-file {} \
--output-dir {}".format(out_rootedtree, dada2_out_t, de, m, out_dir)
os.system(cmd)

print()
print("--------------finish--------------")


################################# alpha diversity ###############################
##### 그룹간의 얼마나 다른지 확인하고 싶은 경우 #####
print()
print("--------------alpha group significance--------------")

out_alphagroup = o + '/' + n + '/4-diversity/' + n + '_' + e + '-group-sig.qzv'

cmd = "/home/grlee/anaconda3/envs/qiime2-2022.2/bin/qiime diversity alpha-group-significance \
--i-alpha-diversity {} \
--m-metadata-file {} \
--o-visualization {}".format(out_dir + '/' + e + '.qza', m, out_alphagroup)
os.system(cmd)

print()
print("--------------finish--------------")




################################# beta diversity ###############################
print()
print("--------------beta group significance--------------")

out_beta = o + '/' + n + '/4-diversity/' + n + '_' + a + '-group-sig.qzv'

cmd = "/home/grlee/anaconda3/envs/qiime2-2022.2/bin/qiime diversity beta-group-significance \
--i-distance-matrix {} \
--m-metadata-file {} \
--m-metadata-column {} \
--p-method permanova \
--p-pairwise \
--o-visualization {}".format(out_dir + '/' + a + '.qza', m, c, out_beta)
os.system(cmd)

print()
print("--------------finish--------------")
