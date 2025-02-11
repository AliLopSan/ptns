import pandas as pd
import sys
from ast import literal_eval
from itertools import combinations_with_replacement


method = ['sankoff','genesis']
penaliz_type = ["equal","loss_half","hgt_half","loss_quarter","hgt_quarter"]
sdata = pd.read_csv("./real_data/interphylum_species_50.csv")

for m in method:
    for p in penaliz_type:
        data_name = "./real_data/"+m+"/penalizations/" + p + "/FA_"+m+"_"+ p +".csv"
        data = pd.read_csv(data_name,converters={'fa_list': literal_eval})

        pairs = dict()
        fa_lists = list(data['fa_list'])

        for fas in fa_lists:
            comparisons = list(combinations_with_replacement(fas, 2))
            for a,b in comparisons:
                if a != b :
                    temp = [a,b]
                    temp = sorted(temp)
                    key = '-'.join(temp)
                    if key not in pairs:
                        pairs[key] = 0
                    else:
                        pairs[key] = pairs[key] + 1

        hw = pd.DataFrame()
        a = []
        b = []
        w = []

        for key in pairs:
            if pairs[key] > 0:
                vals = key.split("-")
                a.append(vals[0])
                b.append(vals[1])
                w.append(pairs[key])
        hw['a'] = a
        hw['b'] = b
        hw['weight'] = w

        res = "./results/"+m+"_"+p+"_hw_info.csv"
        
        hw.to_csv(res,index=False)

