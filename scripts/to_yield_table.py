# Convert the `allinone.C` output text file to a table in .csv format
# python ./to_yield_table.py /path/to/input/text/file /path/to/output/csv/file

import csv
import sys
import numpy as np

text_path = sys.argv[1]
output_path = sys.argv[2]

print("reading:", text_path)
with open(text_path) as f:
    contents = f.read()
    contents_split = contents.split("STORE: ")[1:]
    contents_split = [c.split(" END STORE")[0] for c in contents_split]


header = ["Type", "Energy", "Mass", "Efficiency"]
print("writing:", output_path)
with open(output_path, 'w', encoding='UTF8', newline='') as fo:
    writer = csv.writer(fo)
    writer.writerow(header)
    for cont in contents_split:
        store_lt = []
        cont = cont.split('/')[-1]
        info, eff = cont.split('; e=')
        info = info.split('_reco.root')[0]
        info = info.split('_')

        if info[0][0] == "s":
            store_lt.append(info[0][0] + "_" + info[1][0])
            store_lt.append(info[2].split('E-')[-1])
            store_lt.append(info[3].split('m-')[-1])
        
        elif info[0][0] == "b":
            store_lt.append(info[0][0] + "_" + info[1])
            store_lt.append(info[2].split('E-')[-1])
            store_lt.append(np.nan)
        

        store_lt.append(eff)
        writer.writerow(store_lt)
