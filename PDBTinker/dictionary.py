import numpy as np

def parameter_builder(prmfile):
    params = open(prmfile,"r")
    params = params.readlines()
    biotypes = []
    for line in params:
        if "biotype" in line:
            biotypes.append(line)

    params_dict={}
    for item in biotypes:
        res_key = item.split("\"")[1].replace("N-Terminal ","N").replace("C-Terminal ","C").replace(" ","_").replace("-","minus").replace("+","plus").replace("\'","prime").replace("3","three").replace("5","five")
        atom_key = item.split("\"")[0].split()[-1].replace("*","\'")
        atom_type = item.split("\"")[2].split()[0]
        if res_key not in params_dict:
            params_dict[res_key] = {}
        if atom_key not in params_dict[res_key]:
            params_dict[res_key][atom_key] = atom_type
    return params_dict

def parameter_loader(filename):
    return np.load(filename+".npy").item()

def parameter_saver(filename,dictionary):
    np.save(filename,dictionary)

def parameter_editor(dictionary):
    keys = dictionary.keys()
    for k in keys:
        print(k)
    key = input("Enter key to change: ")
    if key in keys:
        new_key = input("New key: ")
        dictionary[new_key]=dictionary[key]
        del dictionary[key]
    else:
        print("Key not found.")


key_converter = np.load("dictionaries/key_convert_dict.npy").item()
