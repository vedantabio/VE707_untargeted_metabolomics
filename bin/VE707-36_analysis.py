# -*- coding: utf-8 -*-
"""
Created on Tue Feb 20 21:55:03 2024

@author: lcornacchione
"""

import pandas as pd 
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

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
# Function to calculate -log2(p-value)
def calculate_log_pvalue(pval):
    return -np.log2(pval)

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
        elif len(other_strain_list) >= 1 and len(other_strain_list) <= rarity_level:
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
        elif len(other_strain_list) >= 1 and len(other_strain_list) <= rarity_level:
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
        
    #plot volcano plot of rare metabolites with other strains as well
    #build dataframe to plot
    data_list = []
    for metabolite in rare_produced.index:
        strains = rare_produced.at[metabolite,'Other_Producing_Strains']
        strains.append(strain_name)
        for strain in strains:
            data_list.append({
                'ratio': ratio_df.at[metabolite, strain],
                'pval': pval_df.at[metabolite, strain],
                'strain': strain,
                'metabolite': metabolite})
    for metabolite in rare_consumed.index:
        strains = rare_consumed.at[metabolite,'Other_Consuming_Strains']
        strains.append(strain_name)
        for strain in strains:
            data_list.append({
                'ratio': ratio_df.at[metabolite, strain],
                'pval': pval_df.at[metabolite, strain],
                'strain': strain,
                'metabolite': metabolite})
    #transform ratio and pval to log2 and -log2 respectively 
    plot_data = pd.DataFrame(data_list)
    plot_data['ratio'] = plot_data['ratio'].replace(0, 0.01)#handle -inf values when transforming 
    plot_data['ratio'] = np.log2(plot_data['ratio'])
    plot_data['pval'] = -(np.log2(plot_data['pval']))
    plot_data.to_excel('rare_plot_data.xlsx')
    #plot for each unique metabolite
    for metabolite in plot_data['metabolite'].unique():
        fig,ax = plt.subplots()
        fig.set_size_inches(6,5)
        metabolite_data = plot_data.loc[plot_data['metabolite']==metabolite,:]
        sns.scatterplot(data=metabolite_data, x='ratio',y='pval', hue='strain', legend='full', alpha=0.5)
        plt.xlabel('log2(Rare VE707-36 Significant Consumption-Production Metabolite Ratios)')
        plt.ylabel('-log2(p-value)')
        ax.set_xticks([-8,-6,-4,-2, 0, 2, 4])
        ax.set_yticks([0,5,10,15,20,25])
        plt.legend(loc='center left',bbox_to_anchor=(1,0.5))
        # Add text annotations for data points where the number of strains is below 4
        for index, row in metabolite_data.iterrows():
            ax.text(row['ratio']+0.1, row['pval'], row['metabolite'], fontsize=6, ha='left', va='center',alpha=0.8)
        fig.tight_layout()
        plt.savefig('../results/VE707-36_rare_plots/'+metabolite+'_volcano_plot.png',dpi=300)
    
#%%
identify_unique_metabolites(ratio_df, pval_df, 'VE707-36')
#%%
#Volcano plot 

# Calculate -log2(p-value) for each metabolite
log_pval_df = pval_df.applymap(calculate_log_pvalue)

# Calculate the -log2(p-value) for VE707-36
log_pval_ve707_36 = log_pval_df['VE707-36']

# Get the ratios for the specified strain
strain_name = 'VE707-36'
production_thresh = 1.5
consumption_thresh = 0.5
pval_level = 0.05
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

strain_produced_metabolites = strain_produced_metabolites + strain_consumed_metabolites
#dataframe for plotting 
plot_data = pd.DataFrame(data=None)
for metabolite in strain_produced_metabolites:
    hunga_ratio = strain_ratios.loc[metabolite]
    plot_data.loc[metabolite,'hungatella_ratio'] = hunga_ratio
    plot_data.loc[metabolite,'log_pval'] = log_pval_ve707_36.loc[metabolite]
    #calculate the number of strains that consume or produce given metabolite to significant level 
    strain_count = 0
    for strain in ratio_df.columns:
        if hunga_ratio >= 1.5:
            if ratio_df.at[metabolite,strain] >= 1.5 and pval_df.at[metabolite,strain] <= pval_level:
                strain_count = strain_count + 1
        if hunga_ratio <= 0.5: 
            if ratio_df.at[metabolite,strain] <= 0.5 and pval_df.at[metabolite,strain] <= pval_level:
                strain_count = strain_count + 1
    plot_data.at[metabolite,'num_strains'] = int(strain_count)

#plot the volcano plot 

#log2 transform the hungatella ratios for better visualization
plot_data['hungatella_ratio'] = plot_data['hungatella_ratio'].replace(0, 0.01)#handle -inf values to smallest value measureable by assay
plot_data['hungatella_ratio'] = np.log2(plot_data['hungatella_ratio'])

# Create a new column to represent the condition
plot_data['strain_condition'] = np.where(plot_data['num_strains'] <= 4, '<= 4 strains', '>= 4 strains')
plot_data.loc[plot_data['num_strains']==1,'strain_condition'] = 'VE707-36 only'

fig,ax = plt.subplots()
fig.set_size_inches(10,9)

sns.scatterplot(data=plot_data, x='hungatella_ratio',y='log_pval', size='num_strains', sizes=(20,200),
                hue='strain_condition', palette={'<= 4 strains': 'green', '>= 4 strains': 'blue', 'VE707-36 only':'red'}, legend='full', alpha=0.5)
plt.xlabel('log2(VE707-36 Significant Consumption-Production Metabolite Ratios)')
plt.ylabel('-log2(p-value)')
plt.legend(loc='center left',bbox_to_anchor=(1,0.5))
# Add text annotations for data points where the number of strains is below 4
for index, row in plot_data.iterrows():
    if row['num_strains'] <= 4:
        ax.text(row['hungatella_ratio']+0.1, row['log_pval'], index, fontsize=6, ha='left', va='center',alpha=0.8)
fig.tight_layout()
plt.savefig('../results/VE707-36_volcano_plot.png',dpi=300)

