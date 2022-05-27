import os

seed_dict_diamond = {}

for diamond_file in os.listdir("./DiamondOut"):
    if diamond_file.endswith("tsv"):
        name = diamond_file.split("_")[0]
        myfile = open(f'./DiamondOut/{diamond_file}', "r")
        myfile.readline()
        myline = myfile.readline()
        while myline:
            myset = set(myline.strip("{").strip().strip("}").replace("'", "").split(", "))
            if seed_dict_diamond.get(name, None) is None:
                seed_dict_diamond[name] = myset
            else:
                seed_dict_diamond[name] |= myset
            myline = myfile.readline()
        myfile.close()

seed_dict_biosteiner = {}
seed_dict_biosteiner2 = {}
seed_dict_biosteiner9 = {}
for biosteiner_file in os.listdir("./BiosteinerOut"):
    if biosteiner_file.endswith("txt") and ("0.1" in biosteiner_file):
        name = biosteiner_file.split("_")[0]
        myfile = open(f'./BiosteinerOut/{biosteiner_file}', "r")
        myfile.readline()
        myfile.readline()
        myfile.readline()
        myline = myfile.readline()
        while myline:
            myset = set(myline.strip("{").strip().strip("}").replace("'", "").split(", "))
            if seed_dict_biosteiner.get(name, None) is None:
                seed_dict_biosteiner[name] = myset
            else:
                seed_dict_biosteiner[name] |= myset
            myline = myfile.readline()
        myfile.close()
    elif biosteiner_file.endswith("txt") and ("0.2" in biosteiner_file):
        name = biosteiner_file.split("_")[0]
        myfile = open(f'./BiosteinerOut/{biosteiner_file}', "r")
        myfile.readline()
        myfile.readline()
        myfile.readline()
        myline = myfile.readline()
        while myline:
            myset = set(myline.strip("{").strip().strip("}").replace("'", "").split(", "))
            if seed_dict_biosteiner2.get(name, None) is None:
                seed_dict_biosteiner2[name] = myset
            else:
                seed_dict_biosteiner2[name] |= myset
            myline = myfile.readline()
        myfile.close()
    elif biosteiner_file.endswith("txt") and ("0.9" in biosteiner_file):
        name = biosteiner_file.split("_")[0]
        myfile = open(f'./BiosteinerOut/{biosteiner_file}', "r")
        myfile.readline()
        myfile.readline()
        myfile.readline()
        myline = myfile.readline()
        while myline:
            myset = set(myline.strip("{").strip().strip("}").replace("'", "").split(", "))
            if seed_dict_biosteiner9.get(name, None) is None:
                seed_dict_biosteiner9[name] = myset
            else:
                seed_dict_biosteiner9[name] |= myset
            myline = myfile.readline()
        myfile.close()

jaccard_dict = {}
for key, set_biosteiner in seed_dict_biosteiner.items():
    set_diamond = seed_dict_diamond[key]
    jaccard = len(set_diamond.intersection(set_biosteiner)) / len(set_biosteiner)
    jaccard_dict[key] = jaccard

with open("diamond_biosteiner.csv", "w") as f:
    f.write("seed_file,jaccard\n")
    for key, value in jaccard_dict.items():
        f.write(f'{key},{value}\n')
    f.close()

jaccard_dict = {}
for key, set_biosteiner in seed_dict_biosteiner2.items():
    set_diamond = seed_dict_diamond[key]
    jaccard = len(set_diamond.intersection(set_biosteiner)) / len(set_biosteiner)
    jaccard_dict[key] = jaccard

with open("diamond_biosteiner2.csv", "w") as f:
    f.write("seed_file,jaccard\n")
    for key, value in jaccard_dict.items():
        f.write(f'{key},{value}\n')
    f.close()

jaccard_dict = {}
for key, set_biosteiner in seed_dict_biosteiner9.items():
    set_diamond = seed_dict_diamond[key]
    jaccard = len(set_diamond.intersection(set_biosteiner)) / len(set_biosteiner)
    jaccard_dict[key] = jaccard

with open("diamond_biosteiner9.csv", "w") as f:
    f.write("seed_file,jaccard\n")
    for key, value in jaccard_dict.items():
        f.write(f'{key},{value}\n')
    f.close()
