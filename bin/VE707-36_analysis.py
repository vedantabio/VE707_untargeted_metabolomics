# -*- coding: utf-8 -*-
"""
Created on Tue Feb 20 21:55:03 2024

@author: lcornacchione
"""

import pandas as pd 

#read in ratio and p value dataframes
ratio_df = pd.read_excel('../data/VEDA-01-23VS HEAT MAP.XLSX',sheet_name='Heat_map_for_analysis',header=2,index_col=0)
pval_df = pd.read_excel('../data/VEDA-01-23VS HEAT MAP.XLSX',sheet_name='p_values',header=3,index_col=0)

#get strain alias information and rename columns of dataframes to VE707-XX format
strain_df = pd.read_excel('../data/707_strainalias_strainID_key.xlsx')
strains = strain_df['strainID'].values
strain_dict = {}
for strain in strain_df['strain alias'].values:
    strain_dict[strain] = strain_df.loc[strain_df['strain alias']==strain,'strainID'].values[0]
#add new e coli info to strain dict 
strain_dict['P88D4_1'] = 'VE707-43'

#rename df columns to be the VE707-XX format 
ratio_df.rename(mapper=strain_dict,axis=1,inplace=True)
pval_df.rename(mapper=strain_dict,axis=1,inplace=True)

#%%

def identify_unique_metabolites(ratio_df, pval_df, strain_name, consumption_thresh=0.5, production_thresh=1.5, pval_level=0.05, rarity_level=3):
    """""
    Takes a dataframe comprised of ratios of metabolites relative to control medium for a given strain,
    identifies the metabolites that are produced and consumed by strain of interest (strain_name arg),
    then compares those lists to all other strains in the dataframe to identify unique and rare 
    (arg defined) consumed and produced metabolites. Take a pval_df argument which is a 
    matched dataframe to ratio_df, but contains the p value for the metabolite strain ratio. 
    also takes pval_level which is default set to 0.05 and rarity level which is default set to 3.
    rarity level specifies the number of other strains that produce/consume for a metabolite to be called
    as rare. It is inclusive (i.e. default 3 or fewer strains).
    Consumption and production thresholds default to 0.5 and 1.5 respectively, but can be adjusted.
    
    """    
    # Filter out the specified strain from the DataFrame
    other_strains = ratio_df.drop(columns=[strain_name])
        
    # Get the ratios for the specified strain
    strain_ratios = ratio_df[strain_name]
        
    # build a list of the consumed and produced metabolites for the strain
    strain_produced_metabolites = []
    strain_consumed_metabolites = []
    for metabolite, ratio in strain_ratios.items():
        if ratio >= production_thresh:
            p_val = pval_df.at[metabolite,strain_name]
            if p_val <= pval_level:
                strain_produced_metabolites.append(metabolite)
        if ratio <= consumption_thresh:
            p_val = pval_df.at[metabolite,strain_name]
            if p_val <= pval_level:
                strain_consumed_metabolites.append(metabolite)
                
    #identify unique and rare produced metabolites
    unique_produced = []
    rare_produced = pd.DataFrame(data=None,columns=['Other_Producing_Strains'])
    rare_produced.index.name = 'Metabolite'
    for metabolite in strain_produced_metabolites:
        other_strain_list = []
        for strain in other_strains.columns:
            if other_strains.at[metabolite,strain] >= production_thresh and pval_df.at[metabolite,strain] <= pval_level:
                other_strain_list.append(strain)
        if len(other_strain_list) < 1:
            unique_produced.append(metabolite)
        elif len(other_strain_list) > 1 and len(other_strain_list) <= rarity_level:
            rare_produced.loc[metabolite,'Other_Producing_Strains'] = other_strain_list
    unique_produced_series = pd.Series(unique_produced)
    unique_produced_series.name = 'Metabolite'
    
    #identify unique and rare consumed metabolites
    unique_consumed = []
    rare_consumed = pd.DataFrame(data=None,columns=['Other_Consuming_Strains'])
    rare_consumed.index.name = 'Metabolite'
    for metabolite in strain_consumed_metabolites:
        other_strain_list = []
        for strain in other_strains.columns:
            if other_strains.at[metabolite,strain] <= consumption_thresh and pval_df.at[metabolite,strain] <= pval_level:
                other_strain_list.append(strain)
        if len(other_strain_list) < 1:
            unique_consumed.append(metabolite)
        elif len(other_strain_list) > 1 and len(other_strain_list) <= rarity_level:
            rare_consumed.loc[metabolite,'Other_Consuming_Strains'] = other_strain_list
    unique_consumed_series = pd.Series(unique_consumed)
    unique_consumed_series.name = 'Metabolite'
    
    #capture parameters used to store in first sheet of the exported excel
    run_info = pd.DataFrame({'set_point':[consumption_thresh,production_thresh,pval_level,rarity_level]}, 
                            index=['consumption_threshold','production_threshold','p_value_threshold','rarity_threshold'])
    with pd.ExcelWriter('../results/unique_and_rare_metabolites_'+strain_name+'.xlsx') as writer:
        run_info.to_excel(writer, sheet_name='run_info',index=True)
        rare_produced.to_excel(writer, sheet_name='rare_produced',index=True)
        rare_consumed.to_excel(writer, sheet_name='rare_consumed',index=True)
        unique_produced_series.to_excel(writer, sheet_name='unique_produced',index=True)
        unique_consumed_series.to_excel(writer, sheet_name='unique_consumed',index=True)
#%%
identify_unique_metabolites(ratio_df, pval_df, 'VE707-36')

