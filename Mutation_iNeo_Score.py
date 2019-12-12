# -*- coding: utf-8 -*-
import pandas as pd 
import matplotlib.pyplot as plt
import math
import numpy as np
from decimal import Decimal
from openpyxl import load_workbook
import warnings
import sys
import argparse
import copy
import xlrd
import xlwt
import yaml
import os

try:
    parser = argparse.ArgumentParser(description='Parameters for this formula')
    parser.add_argument('-i',type=str)
    parser.add_argument('-o',type=str)
    args = parser.parse_args()
    Filename = args.i
    OutFilename = args.o
except:
    print('input or output error')
    sys.exit(1)

params_file = open('/export/apps/iNeo_score_formula/params.test.yaml','r',encoding='utf-8')
params_content = params_file.read()
params = yaml.load(params_content)
#归一化函数
def normalization(data):
    result = []
    for i in data:
        norm = (i - min(data))/(max(data) -min(data))
        result.append(norm)
    return result
#激励函数
def Sigmoid(x):
    y = 1/(1+ np.exp(-10*x))
    y = y - 0.5
    return y
#生长曲线函数
def Grow(x,L,a,b):
    y = L/(1 + a*np.exp(-b*x))
    return y
#转化ABC
def Mut_lev_ABC(x):
    if(x == 'A'):
        return 0.8802
    if(x == 'B'):
        return 0.0855
    if(x == 'C'):
        return 0.001
    if(x == '.'):
        return 0.001
    if(x not in 'ABC'):
        return 0.001
#二选一
def choose_exist(x,y):
    if(x == '.' and y != '.'):
        return y
    elif(y == '.' and x != '.'):
        return x
    elif(x == '.' and y == '.'):
        return 0
    elif(x != '.' and y!= '.'):
        return max(x,y)
#计算排名
def Rank(x,array):
    array = sorted(array,reverse=False)
    x_index = array.index(x)
    rank = x_index/len(array)
    return rank
#将缺失值转化为列最小值
def missing_val(col):
    temp = [i for i in col if i != '.']
    col_min = min(temp)
    col_Copy = copy.copy(col)
    for i in range(len(Mut_level)):
        if(col_Copy[i] == '.'):
            col_Copy[i] = col_min
        else:
            pass
    return col_Copy
#pd中打印出整行
def get_line(x):
    y = ''
    for i in x:
        y = str(y) + '\t' + str(i)
    return(y)

pre_neo_score = []
pre_neo_score2 = []
bonus = []

input = pd.read_csv(Filename,sep='\t',header=0,index_col=False)
I_AC = list(input['I_AC'])
II_AC = list(input['II_AC'])
AC = []
for i in range(len(I_AC)):
    AC.append(float(I_AC[i]) + float(II_AC[i]))
k = 3/(max(AC)-min(AC))

DNA_AF_array = []
RNA_AF_array = []
Exp_array = []

for i in range(len(input)):
    DNA_AF_array.append(float(choose_exist(input['DNA_AF'][i],input['DNA_bamcount_AF'][i])))
    RNA_AF_array.append(float(choose_exist(input['RNA_AF'][i],input['RNA_bamcount_AF'][i])))
    Exp_array.append(float(choose_exist(input['RPKM'][i],input['TCGA_exp'][i])))

output_data = pd.DataFrame(columns=['sample_num','Chr','Position','Ref', 'Alt','Mut_type','AA_change','Gene_name','ENSG',
                     'DNA_AF', 'DNA_AF_Rank','AF_calculate','RNA_AF','RNA_AF_Rank','RNAF_calculate', 
                     'M', 'M_calculate','DNA_Mut_level','RNA_Mutation_level','RPKM', 'RPKM_Rank',
                     'RPKM_calculate','MHC_I_epitope','MHC_I_Rank','MHC_I_calculate','MHC_II_epitope', 'MHC_II_Rank','MHC_II_calculate',
                     'DG', 'I_AC','I_AC_rank','II_AC','II_AC_rank','I_SB','I_SB_rank','II_SB','II_SB_rank','pre_Neo_Score_2(DNA_AF_cal*RNA_AF_cal*RPKM_cal*MHC_I_calculate)','AC(I_AC+II_AC)','AC_norm','Bonus(AC_norm+MHC_II_calculate*0.5)'])

for i in range(len(input)):
    line = input.iloc[i][:]
    sample_num = line.sample_num
    Chr = line.Chr
    Position = line.Position
    Ref = line.Ref
    Alt = line.Alt
    Mut_type = line.Mut_type
    AA_change = line.AA_change
    Gene_name = line.Gene_name
    ENSG = line.ENSG
    Mut_AA_num = line.Mut_AA_num
    DNA_AF = line.DNA_AF
    DNA_bamcount_AF = line.DNA_bamcount_AF
    RNA_AF = line.RNA_AF
    RNA_bamcount_AF = line.RNA_bamcount_AF
    DNA_Mut_lev = line.DNA_Mut_lev
    RNA_Mut_lev = line.RNA_Mut_lev
    RPKM = line.RPKM
    TCGA_exp = line.TCGA_exp
    Is_driver_gene = 0
    num_mut_I_SB_WB_epitope = line.num_mut_I_SB_WB_epitope
    num_mut_II_SB_WB_epitope = line.num_mut_II_SB_WB_epitope
    I_AC = line.I_AC
    II_AC = line.II_AC
    num_mut_I_SB_epitope = line.num_mut_I_SB_epitope
    num_mut_II_SB_epitope = line.num_mut_II_SB_epitope
#transform value
    DNA_AF_one = float(choose_exist(DNA_AF,DNA_bamcount_AF))
    RNA_AF_one = float(choose_exist(RNA_AF,RNA_bamcount_AF))
    Exp_one = float(choose_exist(RPKM,TCGA_exp))
    DNA_AF_trans = Grow(float(DNA_AF_one),params['DNA_AF_paras']['L'],params['DNA_AF_paras']['a'],params['DNA_AF_paras']['b'])
    RNA_AF_trans = Grow(float(RNA_AF_one),params['RNA_AF_paras']['L'],params['RNA_AF_paras']['a'],params['RNA_AF_paras']['b'])
    AF_trans = DNA_AF_trans * RNA_AF_trans
    AF_trans = AF_trans 
    Exp_trans = Grow(float(Exp_one)/100,params['RPKM_paras']['L'],params['RPKM_paras']['a'],params['RPKM_paras']['b'])
    if(float(Exp_one) <= 2):
        Exp_trans = 0.4*Exp_trans
    MHC_I_trans = Grow(num_mut_I_SB_WB_epitope,params['MHC_I_paras']['L'],params['MHC_I_paras']['a'],params['MHC_I_paras']['b'])
    Mut_num_trans = math.tanh(Mut_AA_num/10)
    AC_norm = k*((I_AC+II_AC)-min(AC))
    AC_sum = I_AC+II_AC
    MHC_II_trans = Grow(num_mut_II_SB_WB_epitope,params['MHC_II_paras']['L'],params['MHC_II_paras']['a'],params['MHC_II_paras']['b'])
#    DNA_Mut_lev_trans = Mut_lev_ABC(DNA_Mut_lev)
#    RNA_Mut_lev_trans = Mut_lev_ABC(RNA_Mut_lev)
    DNA_Mut_lev_trans = DNA_Mut_lev
    RNA_Mut_lev_trans = RNA_Mut_lev
#calculate rank
    DNA_AF_rank = Rank(float(DNA_AF_one),DNA_AF_array)
    RNA_AF_rank = Rank(float(RNA_AF_one),RNA_AF_array)
    Exp_rank = Rank(float(Exp_one),Exp_array)
    MHC_I_rank = Rank(num_mut_I_SB_WB_epitope,list(input['num_mut_I_SB_WB_epitope']))
    MHC_II_rank = Rank(num_mut_II_SB_WB_epitope,list(input['num_mut_II_SB_WB_epitope']))
    I_AC_rank = Rank(I_AC,list(input['I_AC']))
    II_AC_rank = Rank(II_AC,list(input['II_AC']))
    I_SB_rank = Rank(num_mut_I_SB_epitope,list(input['num_mut_I_SB_epitope']))
    II_SB_rank = Rank(num_mut_II_SB_epitope,list(input['num_mut_II_SB_epitope']))
#formula calculation    
#    temp_neo_score = MHC_I_weight*MHC_I_trans + Exp_weight*Exp_trans + AF_weight*AF_trans + MT_weight*Mut_num_trans + Mut_lev_weight*DNA_Mut_lev_trans*RNA_Mut_lev_trans
    temp_neo_score2 = MHC_I_trans*Exp_trans*AF_trans*Mut_num_trans
#    print(Chr+'\t'+str(Position)+'\t'+str(MHC_I_trans)+'\t'+str(Exp_trans)+'\t'+str(AF_trans))
    if(num_mut_II_SB_WB_epitope == 0):
        MHC_II_trans = 0
    temp_bonus = MHC_II_trans*params['MHC_II_weight'] + AC_norm
#filter
    if(float(DNA_AF_one) < 0.05 or float(RNA_AF_one) < 0.05 or DNA_AF_rank <0.1 or RNA_AF_rank<0.1):
        temp_neo_score = 0
        temp_neo_score2 = 0
        temp_bonus = 0
    if(RPKM<1):
        temp_neo_score = 0
        temp_neo_score2 = 0
        temp_bonus = 0
    if(float(Exp_one) < 0.5):
        temp_neo_score = 0
        temp_neo_score2 = 0
        temp_bonus = 0
    if(num_mut_I_SB_WB_epitope==0):        
        temp_neo_score = 0
        temp_neo_score2 = 0
        temp_bonus = 0
#put data into dataframe
    new = [sample_num,Chr,Position,Ref,Alt,Mut_type,AA_change,Gene_name,ENSG,DNA_AF_one,DNA_AF_rank,DNA_AF_trans,RNA_AF_one,RNA_AF_rank,RNA_AF_trans,Mut_AA_num,
          Mut_num_trans,DNA_Mut_lev_trans,RNA_Mut_lev_trans,Exp_one,Exp_rank,Exp_trans,num_mut_I_SB_WB_epitope,MHC_I_rank,MHC_I_trans,num_mut_II_SB_WB_epitope,
          MHC_II_rank,MHC_II_trans,Is_driver_gene,I_AC,I_AC_rank,II_AC,II_AC_rank,num_mut_I_SB_epitope,I_SB_rank,num_mut_II_SB_epitope,
            II_SB_rank,temp_neo_score2,AC_sum,AC_norm,temp_bonus]
    
    output_data.loc[i] = new
    
#    pre_neo_score.append(temp_neo_score)
    pre_neo_score2.append(temp_neo_score2)
    bonus.append(temp_bonus)

#pre_neo_score_range = max(pre_neo_score) - min(pre_neo_score)
pre_neo_score_range2 = max(pre_neo_score2) - min(pre_neo_score2)
pre_neo_score_norm = normalization(pre_neo_score)
pre_neo_score_norm2 = normalization(pre_neo_score2)
bonus_norm = normalization(bonus)
output_data['range(max(pre_Neo_Score_2) - min(pre_Neo_Score_2))'] = pre_neo_score_range2
#output_data['pre_score_norm'] = pre_neo_score_norm
output_data['pre_score_norm_2'] = pre_neo_score_norm2
output_data['bonus_norm'] = bonus_norm
#output_data['Neo_score'] = output_data['pre_Neo_Score'] + pre_neo_score_range*0.1*output_data['bonus_norm']
output_data['Neo_score_2(pre_Neo_Score_2 + range*0.1*bonus_norm)'] = output_data['pre_score_norm_2'] + pre_neo_score_range2*params['bonus_weight']*output_data['bonus_norm']
#output_data['Neo_score_norm'] = normalization(output_data['Neo_score'])
output_data['Neo_score_norm_2'] = normalization(output_data['Neo_score_2(pre_Neo_Score_2 + range*0.1*bonus_norm)'])

for i in range(len(output_data['Neo_score_norm_2'])):
#    if(output_data.loc[i,'DNA_AF_Rank'] >= 0.5 and output_data.loc[i,'RNA_AF_Rank'] >= 0.5 and output_data.loc[i,'RPKM_Rank'] >= 0.5
#      and output_data.loc[i,'MHC_I_Rank'] >= 0.5 and output_data.loc[i,'MHC_II_Rank'] >= 0.5):
#        output_data.loc[i,'Neo_score_norm'] = 10*output_data.loc[i,'Neo_score_norm']
    if(output_data.loc[i,'I_SB'] > 0 or output_data.loc[i,'II_SB'] > 0):
#        output_data.loc[i,'Neo_score_norm'] = 2*output_data.loc[i,'Neo_score_norm']
        I_SB_w = output_data.loc[i,'I_SB']
        if output_data.loc[i,'I_SB'] > 20:
            I_SB_w = 20
        II_SB_w = output_data.loc[i,'II_SB']
        if output_data.loc[i,'II_SB'] > 20:
            II_SB_w = 20
        I_AC_w = output_data.loc[i,'I_AC']
        II_AC_w = output_data.loc[i,'II_AC']
        output_data.loc[i,'Neo_score_norm_2(Neo_score_norm_2*SB_bonus)'] = output_data.loc[i,'Neo_score_norm_2']*(1+params['I_SB_weight']*I_SB_w + params['II_SB_weight']*II_SB_w)
        output_data.loc[i,'SB_bonus'] = 1+params['I_SB_weight']*I_SB_w + params['II_SB_weight']*II_SB_w
    else:
        output_data.loc[i,'Neo_score_norm_2(Neo_score_norm_2*SB_bonus)'] = output_data.loc[i,'Neo_score_norm_2']
        output_data.loc[i,'SB_bonus'] = 0
#output_data['Neo_score_norm_3(Neo_score_norm_2*SB_bonus)'] = normalization(output_data['Neo_score_norm_2'])
for i in range(len(output_data['Neo_score_norm_2(Neo_score_norm_2*SB_bonus)'])):
    zero_num = 0
    if(output_data.loc[i,'I_SB']==0):
        zero_num = zero_num + 1
    if(output_data.loc[i,'II_SB']==0):
        zero_num = zero_num + 1
    if(output_data.loc[i,'I_AC']==0):
        zero_num = zero_num + 1
    if(output_data.loc[i,'II_AC']==0):
        zero_num = zero_num + 1
    if(zero_num>=3):
        output_data.loc[i,'Neo_score_norm_2(Neo_score_norm_2*SB_bonus)'] = params['diminish']*output_data.loc[i,'Neo_score_norm_2(Neo_score_norm_2*SB_bonus)']
output_data.sort_values('Neo_score_norm_2(Neo_score_norm_2*SB_bonus)',inplace=True,ascending=False)
result = pd.ExcelWriter(OutFilename)
output_data.to_excel(excel_writer=result,sheet_name='Neo_score',startcol=0,index=False)
result.save()
