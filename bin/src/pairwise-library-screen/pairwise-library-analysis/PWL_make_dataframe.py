'''
Inputs: outputs from PWL_screen_analysis.py, annotated library design file
Outputs: pandas dataframe 'indel_df' with computed off:on ratios

Josh Tycko
'''

import pandas as pd
import os
import datetime
import numpy as np

def var(df):
	if np.isnan(df['WIR StDev BR1'])or np.isnan(df['WIR StDev BR2']):
		return np.nanmax([df['WIR StDev BR1'],df['WIR StDev BR2']])
	return np.sqrt(float(df['WIR StDev BR1']**2 + df['WIR StDev BR2']**2)/4)

def offonratio(dfrow):
	ontarget = df[(df['GuideID']==dfrow['GuideID']) & (df['On Target']==1)]

	#if guide did not work on target > 2%, do not calculate specificity ratio
	if float(ontarget['Avg Weighted Indel Rate'].mean()) < 0.02: 
		return ''
	try:
		return float(dfrow['Avg Weighted Indel Rate'])/float(ontarget['Avg Weighted Indel Rate'].mean())
	except ZeroDivisionError:
		return ''

def offonratio2(dfrow):
	ontarget = df2[(df2['GuideID']==dfrow['GuideID']) & (df2['On Target']==1)]

	#if guide did not work on target > 2%, do not calculate specificity ratio
	if float(ontarget['Avg Weighted Indel Rate'].mean()) < 0.02: 
		return ''
	try:
		return float(dfrow['Avg Weighted Indel Rate'])/float(ontarget['Avg Weighted Indel Rate'].mean()) 
	except ZeroDivisionError:
		return ''

libdf = pd.read_csv('library_features_annotated.csv')
libdf = libdf.set_index(['Hammingbarcode'])

indel_dfBR3D03 = pd.read_pickle('Sample_BR3D03/BR3D03_indel_df')
indel_dfBR4D03 = pd.read_pickle('Sample_BR4D03/BR4D03_indel_df')
indel_dfBR3D14 = pd.read_pickle('Sample_BR3D14/BR3D14_indel_df')
indel_dfBR4D14 = pd.read_pickle('Sample_BR4D14/BR4D14_indel_df')

indel_dfBR3D03.index.name = 'Hammingbarcode'
indel_dfBR4D03.index.name = 'Hammingbarcode'
indel_dfBR3D14.index.name = 'Hammingbarcode'
indel_dfBR4D14.index.name = 'Hammingbarcode'

#merge features and indel rates on Hamming Codes
df = pd.merge(libdf,indel_dfBR3D03,left_index = True, right_index = True, how='left') 
df = pd.merge(df,indel_dfBR4D03,left_index = True, right_index = True, how='left', suffixes=(' BR1',' BR2')) #add next biorep
df['Day'] = 3

#to add timepoint as more rows
df2 = pd.merge(libdf,indel_dfBR3D14,left_index = True, right_index = True, how='left') #merge features and indel rates on Hamming Codes
df2 = pd.merge(df2,indel_dfBR4D14,left_index = True, right_index = True, how='left', suffixes=(' BR1',' BR2')) #add next biorep
df2['Day'] = 14

for dfx in [df, df2]:
	dfx['Avg Weighted Indel Rate'] = dfx[['Weighted Indel Rate BR1','Weighted Indel Rate BR2']].mean(axis = 1)
	dfx['StDev Weighted Indel Rate'] = dfx.apply(var,axis = 1)
df['Off:On Target Ratio'] = df.apply(offonratio, axis = 1)
df2['Off:On Target Ratio'] = df2.apply(offonratio2, axis = 1)

#merge along rows
df = pd.concat([df,df2]) 

df[['Good Reads BR1', 'Indel Reads BR1', 'rBCs pre-validation BR1', 'rBCs BR1', 'Good Reads BR2', 'Indel Reads BR2', 'rBCs pre-validation BR2', 'rBCs BR2']].replace('',0)
df['Total Reads'] = df[['Good Reads BR1','Good Reads BR2']].sum(axis = 1)
df['Total rBCs'] = df[['rBCs BR1','rBCs BR2']].sum(axis = 1)
print df.describe()
print df['Off:On Target Ratio'].describe()
df.to_csv('D03+D14_'+str(datetime.date.today())+'_database.csv')

