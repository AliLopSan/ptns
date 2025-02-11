import pandas as pd
import matplotlib.pyplot as plt
from pySankey.sankey import sankey

plt.rcParams.update({
    "text.usetex": True,
    "font.family": "Helvetica"
})

def plot_sankey(species_data,hw_info,name,t_threshold):
    filtered = hw_info[hw_info['weight'] > t_threshold]
    KEGG_2_species = dict(zip(species_data['kegg_id'],species_data['species']))
    
    # Identify entries matching the pattern 'T' followed by numbers
    a_series = filtered['a']
    pattern = a_series.str.match(r'^T\d+$')
    # Replace matching entries with their values from the KEGG dictionary
    a_series[pattern] = a_series[pattern].map(KEGG_2_species)
    filtered['a'] = a_series
    
    # Identify entries matching the pattern 'T' followed by numbers
    b_series = filtered['b']
    pattern  = b_series.str.match(r'^T\d+$')
    # Replace matching entries with their values from the KEGG dictionary
    b_series[pattern] = b_series[pattern].map(KEGG_2_species)
    filtered['b'] = b_series
    
    
    sankey(left=filtered["a"], right=filtered["b"], aspect=25, fontsize=20,
           rightWeight=filtered['weight'],leftWeight=filtered['weight'])
    
    # Get current figure
    fig = plt.gcf()
    # Set size in inches
    fig.set_size_inches(12, 22)
    # Set the color of the background to white
    fig.set_facecolor("w")
    # Save the figure
    fig.savefig(name, bbox_inches="tight", dpi=300)
        


method = ['sankoff','genesis']
penaliz_type = ["equal","hgt_half","hgt_quarter"]
sdata = pd.read_csv("./real_data/interphylum_species_50.csv")

t_thresholds = [5,9,18]

for p in penaliz_type:
    data_name = "./results/sankoff_"+ p +"_hw_info.csv"
    data = pd.read_csv(data_name)
    for t in t_thresholds:
        res_name = "./results/plots/sankey/sankoff_"+p+"_sankey_"+"t_"+str(t)+".pdf"
        plot_sankey(sdata,data,res_name,t)

for p in penaliz_type:
    data_name = "./results/genesis_"+ p +"_hw_info.csv"
    data = pd.read_csv(data_name)
    for t in t_thresholds:
        res_name = "./results/plots/sankey/genesis_"+p+"_sankey_"+"t_"+str(t)+".pdf"
        plot_sankey(sdata,data,res_name,t)

        
        
