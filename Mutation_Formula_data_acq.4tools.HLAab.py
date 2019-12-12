#!/usr/bin/env python
# -*- coding: utf-8 -*-
from __future__ import print_function
import os
import sys
import pandas as pd
import numpy as np
import re

sys.path.append('/home/chenrc')
import python_module_sub as pm

mutation_objects = pm.creat_object_for_mutations(sys.argv[1],sys.argv[2])
HLA_I_ab = {}
HLA_I_ab_result = open(sys.argv[3],'r')
for line in HLA_I_ab_result:
	line = line.strip()
	HLA_type = line.split('\t')[0]
	ab = line.split('\t')[2]
	HLA_I_ab[HLA_type] = float(ab)
tumor_purity = float(sys.argv[4])
output_file_name = str(sys.argv[5])
result = pd.ExcelWriter(output_file_name)
output_df = pd.DataFrame(columns=['Chr','Position','Ref','Alt','Mut_type','Mut_AA_num','AA_change','Gene_name','ENSG','DNA_AF','DNA_bamcount_AF','RNA_AF','RNA_bamcount_AF','DNA_Mut_lev','DNA_RF_Mut_lev','RNA_Mut_lev','RPKM','TCGA_exp','Neo_qual','wild_IG','mut_IG','HLA_AB','Is_driver_gene','num_mut_I_SB_WB_epitope','num_mut_II_SB_WB_epitope','I_AC','II_AC','num_mut_I_SB_epitope','num_mut_I_WB_epitope','num_mut_II_SB_epitope','num_mut_II_WB_epitope','I_AC_A','I_AC_B','I_AC_C','II_AC_A','II_AC_B','II_AC_C'])
for i in range(len(mutation_objects)):
	mut = mutation_objects[i]
	DNA_AF = mut.DNA_AF
	RNA_AF = mut.RNA_AF
	DNA_bamcount_AF = mut.DNA_bamcount_AF
	RNA_bamcount_AF = mut.RNA_bamcount_AF
	if(DNA_AF != '.'):
		DNA_AF = float(mut.DNA_AF)*(1/tumor_purity)
	if(mut.DNA_bamcount_AF != '.'):
		DNA_bamcount_AF = float(mut.DNA_bamcount_AF)*(1/tumor_purity)
	if(RNA_AF != '.'):
		RNA_AF = float(mut.RNA_AF)*(1/tumor_purity)
	if(mut.RNA_bamcount_AF != '.'):
		RNA_bamcount_AF = float(mut.RNA_bamcount_AF)*(1/tumor_purity)
	MHC_I_AC_val = mut.I_AC.count('A')*5 + mut.I_AC.count('B')*3 + 2*mut.I_AC.count('C')
	MHC_II_AC_val = mut.II_AC.count('A')*2.5 + mut.II_AC.count('B')*1.5 + mut.II_AC.count('C')
	
	I_mut_no = [i for i in mut.mut_I_aff if 'WB' not in i and 'SB' not in i]	
	I_mut_WB = [i for i in mut.mut_I_aff if 'WB' in i]
	I_mut_SB = [i for i in mut.mut_I_aff if 'SB' in i]
	II_mut_no = [i for i in mut.mut_II_aff if 'WB' not in i and 'SB' not in i]
	II_mut_WB = [i for i in mut.mut_II_aff if 'WB' in i]
	II_mut_SB = [i for i in mut.mut_II_aff if 'SB' in i]

	MHC_I_HLAab_val = 0
	for item in mut.HLA_I:
		if(item == '.'):
			continue
		try:
			MHC_I_HLAab_val = MHC_I_HLAab_val +float(HLA_I_ab[item])
		except KeyError:
			MHC_I_HLAab_val = MHC_I_HLAab_val + 0.1
	mut_list = [mut.Chr,mut.Pos,mut.Ref,mut.Alt,mut.Mut_type,mut.Mut_AA_num,mut.aa_change,mut.gene_name,mut.ENSG,
                        DNA_AF,DNA_bamcount_AF,RNA_AF,RNA_bamcount_AF,mut.DNA_Mut_lev,mut.DNA_RF_Mut_lev,mut.RNA_Mut_lev,
                        mut.RPKM,mut.TCGA_exp,mut.Neo_qual,mut.wild_IG,mut.mut_IG,mut.HLA_AB,mut.Is_driver_gene,MHC_I_HLAab_val,str(len(II_mut_no)+len(II_mut_WB)+len(II_mut_SB)),MHC_I_AC_val,MHC_II_AC_val,str(len(I_mut_SB)),str(len(I_mut_WB)),str(len(II_mut_SB)),str(len(II_mut_WB)),str(mut.I_AC.count('A')),str(mut.I_AC.count('B')),str(mut.I_AC.count('C')),str(mut.II_AC.count('A')),str(mut.II_AC.count('B')),str(mut.II_AC.count('C'))]
	output_df.loc[i] = mut_list
output_df.to_excel(excel_writer=result,sheet_name='Neo_score',startcol=0,index=False)
result.save()
