# Convert the `get_limit.py` output text file to a table in .csv format
# python ./to_V2_table.py /path/to/input/text/file /path/to/output/csv/file

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

# print(contents_split)


header = ["Energy", "Mass", "Maj", "Dir"]
print("writing:", output_path)
with open(output_path, 'w', encoding='UTF8', newline='') as fo:
    writer = csv.writer(fo)
    writer.writerow(header)
    for cont in contents_split:
        store_lt = []
        # print(cont)
        energy, cont = cont.split('e=')[1].split("; m=")
        mass, cont = cont.split("; V2_M=")
        V2_m, V2_d= cont.split("; V2_D=")

        store_lt.append(energy)
        store_lt.append(mass)
        store_lt.append(V2_m)
        store_lt.append(V2_d)
        
        writer.writerow(store_lt)
