
"""
pdbtinker.py
Converts PDBs to Tinker XYZs

Handles the primary functions
"""
import parmed as prm
import sys
import numpy as np
import os
import PDBTinker.dictionaries as dictionaries

dictpath = os.path.dirname(dictionaries.__file__)
amino_acid_list = np.load(dictpath+"/amino_acid_list.npy")
nucleic_acid_list = np.load(dictpath+"/nucleic_acid_list.npy")
water_list = np.load(dictpath+"/water_list.npy")
ion_list = np.load(dictpath+"/ion_list.npy")
non_standard_list = np.load(dictpath+"/non_standard_list.npy")
cap_list = np.load(dictpath+"/cap_list.npy")
param_dict={}
atom_name_dict={}

def interatomic_distance(atom1,atom2):
    x_diff = (atom1.xx - atom2.xx)**2
    y_diff = (atom1.xy - atom2.xy)**2
    z_diff = (atom1.xz - atom2.xz)**2
    return np.sqrt(x_diff + y_diff + z_diff)

def load_pdb(filename):
    system = prm.load_file(filename)
    for atom in system.atoms:
        atom.mass = 0
    disambiguate_residues(system)
    return system

def disambiguate_residues(temp):
    for residue in temp.residues:
        if residue.name == "HIS":
            atomlist=[]
            for atom in residue.atoms:
                atomlist.append(atom.name)
            if "HE2" in atomlist and "HD1" in atomlist:
                residue.name = "HIP"
            elif "HE2" in atomlist and "HD1" not in atomlist:
                residue.name = "HIE"
            elif "HE2" not in atomlist and "HD1" in atomlist:
                residue.name = "HID"

def process_new_parameters(filename):
    global param_dict
    global atom_name_dict
    temp = open(filename,"r")
    params = temp.readlines()
    temp.close()
    for line in params:
        if 'atom' in line[:5]:
    #         print(line)
            key = line.split("\"")[0].split()[1]
            value = line.split("\"")[1]
            value = value.replace("Lysine HN","Lysine HZ").replace("Glycine","GLY").replace("Alanine","ALA").replace("Valine","VAL").replace("Leucine","LEU")
            value = value.replace("Isoleucine","ILE").replace("Serine","SER").replace("Threonine","THR").replace("Cysteine Anion","CYM")
            value = value.replace("Cysteine","CYS").replace("Cystine","CYX").replace("Proline","PRO").replace("Phenylalanine","PHE")
            value = value.replace("Tyrosine Anion","TYD").replace("Tyrosine","TYR").replace("Tryptophan","TRP")
            value = value.replace("Histidine (+)","HIP").replace("Histidine (HD)","HID").replace("Histidine (HE)","HIE")
            value = value.replace("Aspartate","ASP").replace("Aspartic Acid","ASH").replace("Glutamate","GLU").replace("Glutamic Acid","GLH")
            value = value.replace("Asparagine","ASN").replace("Glutamine","GLN").replace("Methionine","MET").replace("Lysine (Neutral)","LYN")
            value = value.replace("Lysine","LYS").replace("Arginine","ARG").replace("Ornithine","ORN")
            value = value.replace("Acetyl Cap","ACE").replace("N-MeAmide Cap","NME").replace("Amide Cap","AMC").replace("N-Terminal PRO","NPRO")
            value = value.replace("N-Terminal","NTERM").replace("C-Terminal","CTERM")
            value = value.replace("Adenine","A").replace("Cytosine","C").replace("Guanine","G").replace("Thymine","T").replace("Uracil","U")
            value = value.replace("Ribose","RN").replace("Deoxyribose","DN").replace("R-Phosphodiester","RPhos").replace("D-Phosphodiester","DPhos")
            value = value.replace("R-5'-Hydroxyl","R5H").replace("R-5'-Phosphate","R5Phos")
            value = value.replace("R-3'-Hydroxyl","R3H").replace("R-3'-Phosphate","R3Phos")
            value = value.replace("D-5'-Hydroxyl","D5H").replace("D-5'-Phosphate","D5Phos")
            value = value.replace("D-3'-Hydroxyl","D3H").replace("D-3'-Phosphate","D3Phos")
            value = value.replace("AMOEBA Water","WAT").replace("Lithium Ion","LI").replace("Sodium Ion","NA").replace("Potassium Ion","K").replace("Rubidium Ion","RB").replace("Cesium Ion","CS").replace("Beryllium Ion","BE").replace("Magnesium Ion","MG").replace("Calcium Ion","CA").replace("Zinc Ion","ZN").replace("Chloride Ion","CL")
            value = value.replace("S-","SG").replace("H2'2","H2''")
            value = value.replace("H5'1","H5'").replace("H5'2","H5''").replace("H2'1","H2'").replace("HO'2","HO''")
            value = value + " " +line.split("\"")[0].split()[3]
            param_dict[tuple(value.split())] = key
    keys = list(param_dict.keys())
    for key in keys:
        if len(key) == 4:
            oldkey = list(key)
            newkeyvar = ''.join([oldkey[0],oldkey[2]])
            oldkey[0] = newkeyvar
            del oldkey[2]
            param_dict[tuple(oldkey)] = param_dict[key]
            del param_dict[key]
    temp = open(filename,"r")
    params = temp.readlines()
    temp.close()
    for line in params:
        if 'atom' in line[:5]:
            key = line.split("\"")[0].split()[1]
            value = line.split("\"")[0].split()[3]
            atom_name_dict[key] = value

def process_n_terminal(residue):
    for atom in residue.atoms:
        if atom.mass == 0:
            if atom.name == "N":
                if ("NTERM","NH3+","N") in param_dict.keys():
                    atom.mass = param_dict[("NTERM","NH3+","N")]
            elif atom.name == "H1" or atom.name == "H2" or atom.name == "H3":
                if ("NTERM","H3N+","H") in param_dict.keys():
                    atom.mass = param_dict[("NTERM","H3N+","H")]

def process_c_terminal(residue):
    for atom in residue.atoms:
        if atom.mass == 0:
            if atom.name == "C":
                atom.mass = param_dict[("CTERM","COO-","C")]
            elif atom.name == "O" or atom.name == "OXT":
                atom.mass = param_dict[("CTERM","COO-","O")]

def process_backbone(residue):
    for atom in residue.atoms:
        if "CY" in residue.name and atom.mass == 0 and atom.name == "CA":
            atom.mass = param_dict[('CYS','CA','CA')]
        if atom.mass == 0:
            if atom.name == "N":
                atom.mass = param_dict[('ALA', 'N', 'N')]
            elif atom.name == "C":
                atom.mass = param_dict[('ALA', 'C', 'C')]
            elif atom.name == "CA":
                atom.mass = param_dict[('ALA', 'CA', 'CA')]
            elif atom.name == "O":
                atom.mass = param_dict[('ALA', 'O', 'O')]
            elif atom.name == "H":
                atom.mass = param_dict[('ALA', 'HN', 'HN')]
            elif atom.name == "HA":
                atom.mass = param_dict[('ALA', 'HA', 'H')]

def process_glycine(residue):
    for atom in residue.atoms:
        if atom.mass == 0:
            if atom.name == "N":
                atom.mass = param_dict[('GLY', 'N', 'N')]
            elif atom.name == "C":
                atom.mass = param_dict[('GLY', 'C', 'C')]
            elif atom.name == "CA":
                atom.mass = param_dict[('GLY', 'CA', 'CA')]
            elif atom.name == "O":
                atom.mass = param_dict[('GLY', 'O', 'O')]
            elif atom.name == "H":
                atom.mass = param_dict[('GLY', 'HN', 'HN')]
            elif atom.name == "HA2" or atom.name == "HA3":
                atom.mass = param_dict[('GLY', 'HA', 'H')]

def process_n_term_proline(residue):
    for atom in residue.atoms:
        if atom.mass == 0:
            if atom.name == "N":
                atom.mass = param_dict[("NPRO","NH2+","N")]
            elif atom.name == "H1" or atom.name == "H2" or atom.name == "H3":
                atom.mass = param_dict[("NPRO","H2N+","HN")]
            elif atom.name == "CA":
                atom.mass = param_dict[('NPRO', 'CA', 'CA')]
            elif atom.name == "C":
                atom.mass = param_dict[('NPRO', 'C', 'C')]
            elif atom.name == "O":
                atom.mass = param_dict[('NPRO', 'O', 'O')]
            elif atom.name == "HA":
                atom.mass = param_dict[('NPRO', 'HA', 'H')]
            elif atom.name == "CD":
                atom.mass = param_dict[('NPRO', 'CD', 'C')]
            elif atom.name == "HD2" or atom.name == "HD3":
                atom.mass = param_dict[('NPRO', 'HD', 'H')]

def get_key(resname,atomname):
    keys = param_dict.keys()
    for key in keys:
        if resname == key[0] and atomname == key[1]:
            return key
    for key in keys:
        if resname == key[0] and atomname[:-1] == key[1]:
            return key
    for key in keys:
        if resname == key[0] and atomname[:-2] == key[1]:
            return key
    for key in keys:
        if resname == "CYX":
            return get_key("CYS",atomname)
    for key in keys:
        if resname == "LYD":
            return get_key("LYS",atomname)
    for key in keys:
        if resname == "TYD":
            return get_key("TYR",atomname)
    return 0

def process_side_chain(residue):
    r_name = residue.name
    for atom in residue.atoms:
        a_name = atom.name
        if atom.mass == 0:
            key = get_key(r_name,a_name)
            if r_name == "THR" and atom.name == "HG1":
                key = get_key("THR","H")
            if key == 0:
                atom.mass = 0
            else:
                atom.mass = param_dict[key]

def process_proline(residue):
    for atom in residue.atoms:
        if atom.mass == 0:
            if atom.name == "N":
                atom.mass = param_dict[("PRO","N","N")]
            elif atom.name == "H1" or atom.name == "H2" or atom.name == "H3":
                atom.mass = param_dict[("PRO","H2N+","HN")]
            elif atom.name == "CA":
                atom.mass = param_dict[('PRO', 'CA', 'CA')]
            elif atom.name == "C":
                atom.mass = param_dict[('PRO', 'C', 'C')]
            elif atom.name == "O":
                atom.mass = param_dict[('PRO', 'O', 'O')]
            elif atom.name == "HA":
                atom.mass = param_dict[('PRO', 'HA', 'H')]
            elif atom.name == "CD":
                atom.mass = param_dict[('PRO', 'CD', 'C')]
            elif atom.name == "HD2" or atom.name == "HD3":
                atom.mass = param_dict[('PRO', 'HD', 'H')]

def process_amino_acid(residue):
    atomlist = []
    for atom in residue.atoms:
        atomlist.append(atom.name)
    if residue.name == "PRO":
        #check if it's n-terminal proline
        if 'H3' in atomlist:
            process_n_term_proline(residue)
#             process_proline(residue)
            #apply n-terminal proline values
        if 'OXT' in atomlist:
            #apply c-terminal backbone
            process_c_terminal(residue)
        process_proline(residue)
    elif residue.name == "GLY":
        #check for n-terminal
        if "H3" in atomlist:
            process_n_terminal(residue)
        #check if it's c-terminal
        if 'OXT' in atomlist:
            process_c_terminal(residue)
        #process backbone
        process_glycine(residue)

    if "H3" in atomlist:
        process_n_terminal(residue)

    #check if it's c-terminal
    if 'OXT' in atomlist:
        process_c_terminal(residue)
    #process backbone
    process_backbone(residue)
    #process side chain atoms
    process_side_chain(residue)

def process_PDB(system):
    global amino_acid_list
    global nucleic_acid_list
    global water_list
    global ion_list
    global cap_list
    for residue in system.residues:
        if residue.name in amino_acid_list:
            process_amino_acid(residue)
        elif residue.name in nucleic_acid_list:
            process_nucleotide(residue)
        elif residue.name in water_list:
            process_water(residue)
        elif residue.name in ion_list:
            process_ion(residue)
        elif residue.name in cap_list:
            process_cap(residue)
    for atom in system.atoms:
        if atom.mass != 0:
            atom.name = atom_name_dict[atom.mass]

def build_xyz(system,filename):
    f = open(filename,"w")
    f.write(str(len(system.atoms))+"\n")
    for atom in system.atoms:
        bondstring = ""
        bondlist = []
        for i in atom.bonds:
            if i.atom1.idx == atom.idx:
                bondlist.append(i.atom2.idx+1)
                # bondstring = bondstring + str(i.atom2.idx+1) + "\t"
            elif i.atom2.idx == atom.idx:
                bondlist.append(i.atom1.idx+1)
                # bondstring = bondstring + str(i.atom1.idx+1) + "\t"
        if atom.name == "SS" and atom.residue.name == "CYX" and len(bondlist) < 2:
            for atom2 in system.atoms:
                if atom2.residue.name == "CYX" and atom2.name == "SS":
                    ia_dist = interatomic_distance(atom,atom2)
                    if ia_dist < 2 and ia_dist > 0:
                        bondlist.append(atom2.idx)
        bondlist.sort()
        bondstring = ''.join(str('{:>6}'.format(x)) for x in bondlist)
        index = atom.idx+1
        name = atom.name
        x = atom.xx
        y = atom.xy
        z = atom.xz
        atomtype = atom.mass
        linestring = "{:>6}  {:<2} {:>12.6f} {:>11.6f} {:>11.6f} {:>5}{}".format(str(index),str(name),x,y,z,str(atomtype),bondstring)
        if atomtype == 0:
            linestring = linestring + "ATOM TYPE NOT FOUND"
        linestring = linestring + "\n"
        f.write(linestring)
    f.close()

def process_nucleotide(residue):
    atomlist = []
    for atom in residue.atoms:
        atomlist.append(atom.name)
    process_nucleic_terminals(residue)
    if "H2''" in atomlist or "H2'2" in atomlist:
        process_deoxy_phosphate(residue)
        process_deoxyribose(residue)
    process_ribose_phosphate(residue)
    process_ribose(residue)
    process_base(residue)

def process_base(residue):
    for atom in residue.atoms:
        a_name = atom.name
        if atom.mass == 0:
            if "A" in residue.name:
                key = get_key("A",a_name)
                if key == 0:
                    atom.mass = 0
                else:
                    atom.mass = param_dict[key]
            if "C" in residue.name:
                key = get_key("C",a_name)
                if key == 0:
                    atom.mass = 0
                else:
                    atom.mass = param_dict[key]
            if "G" in residue.name:
                key = get_key("G",a_name)
                if key == 0:
                    atom.mass = 0
                else:
                    atom.mass = param_dict[key]
            if "T" in residue.name:
                key = get_key("T",a_name)
                if key == 0:
                    atom.mass = 0
                else:
                    atom.mass = param_dict[key]
            if "U" in residue.name:
                key = get_key("U",a_name)
                if key == 0:
                    atom.mass = 0
                else:
                    atom.mass = param_dict[key]

def process_deoxy_phosphate(residue):
    for atom in residue.atoms:
        a_name = atom.name
        if atom.mass == 0 and atom.name =="P":
            key = get_key("DPhos","P")
            if key == 0:
                atom.mass = 0
            else:
                atom.mass = param_dict[key]
        if atom.mass == 0 and "O" in atom.name and "P" in atom.name:
            key = get_key("DPhos","OP")
            if key == 0:
                atom.mass = 0
            else:
                atom.mass = param_dict[key]

def process_ribose_phosphate(residue):
    for atom in residue.atoms:
        a_name = atom.name
        if atom.mass == 0 and atom.name =="P":
            key = get_key("RPhos","P")
            if key == 0:
                atom.mass = 0
            else:
                atom.mass = param_dict[key]
        if atom.mass == 0 and "O" in atom.name and "P" in atom.name:
            key = get_key("RPhos","OP")
            if key == 0:
                atom.mass = 0
            else:
                atom.mass = param_dict[key]

def process_deoxyribose(residue):
    for atom in residue.atoms:
        a_name = atom.name
        if atom.mass == 0:
            if "C" in residue.name or "T" in residue.name:
                key = get_key("DN(CT)",a_name)
                if key == 0:
                    atom.mass = 0
                else:
                    atom.mass = param_dict[key]
            elif "G" in residue.name or "A" in residue.name:
                key = get_key("DN(AG)",a_name)
                if key == 0:
                    atom.mass = 0
                else:
                    atom.mass = param_dict[key]

def process_ribose(residue):
    for atom in residue.atoms:
        a_name = atom.name
        if atom.name == "HO2'":
            a_name = "HO''"
        if atom.mass == 0:
            if "C" in residue.name or "U" in residue.name:
                key = get_key("RN(CU)",a_name)
                if key == 0:
                    atom.mass = 0
                else:
                    atom.mass = param_dict[key]
            elif "G" in residue.name or "A" in residue.name:
                key = get_key("RN(AG)",a_name)
                if key == 0:
                    atom.mass = 0
                else:
                    atom.mass = param_dict[key]

def process_nucleic_terminals(residue):
    for atom in residue.atoms:
        if atom.mass == 0:
            if (atom.name == "HO5'" or atom.name == "H5T") and "D" in residue.name:
                atom.mass = param_dict[('D5H', 'H5T', 'HO')]
            elif atom.name == "HO5'":
                atom.mass = param_dict[('R5H', 'H5T', 'HO')]
            elif (atom.name == "HO3'" or atom.name == "H3T") and "D" in residue.name:
                atom.mass = param_dict[('D3H', 'H3T', 'HO')]
            elif atom.name == "HO3'":
                atom.mass = param_dict[('R3H', 'H3T', 'HO')]
            elif atom.name == "O5'" and "D" in residue.name:
                bondlist = []
                for bond in atom.bonds:
                    if atom.name == bond.atom1.name:
                        bondlist.append(bond.atom2.name)
                    elif atom.name == bond.atom2.name:
                        bondlist.append(bond.atom1.name)
                if "P" not in bondlist:
                    atom.mass = param_dict[('D5H', "O5'T", 'OH')]
            elif atom.name == "O5'":
                bondlist = []
                for bond in atom.bonds:
                    if atom.name == bond.atom1.name:
                        bondlist.append(bond.atom2.name)
                    elif atom.name == bond.atom2.name:
                        bondlist.append(bond.atom1.name)
                if "P" not in bondlist:
                    atom.mass = param_dict[('R5H', "O5'T", 'OH')]
            elif atom.name == "O3'" and "D" in residue.name:
                bondlist = []
                for bond in atom.bonds:
                    if atom.name == bond.atom1.name:
                        bondlist.append(bond.atom2.name)
                    elif atom.name == bond.atom2.name:
                        bondlist.append(bond.atom1.name)
                if "P" not in bondlist:
                    atom.mass = param_dict[('D3H', "O3'T", 'OH')]
            elif atom.name == "O3'":
                bondlist = []
                for bond in atom.bonds:
                    if atom.name == bond.atom1.name:
                        bondlist.append(bond.atom2.name)
                    elif atom.name == bond.atom2.name:
                        bondlist.append(bond.atom1.name)
                if "P" not in bondlist:
                    atom.mass = param_dict[('R3H', "O3'T", 'OH')]

def process_water(residue):
    for atom in residue.atoms:
        if atom.mass == 0:
            if atom.name == "O":
                atom.mass = param_dict[('WAT', 'O', 'O')]
            if atom.name == "H1" or atom.name == "H2":
                atom.mass = param_dict[('WAT', 'H', 'H')]

def process_ion(residue):
    for atom in residue.atoms:
        if atom.mass == 0:
            key = get_ion_key(residue.name)
            if key == 0:
                atom.mass = 0
            else:
                atom.mass = param_dict[key]

def get_ion_key(name):
    keys = param_dict.keys()
    for key in keys:
        if name == key[0]:
            return key
    for key in keys:
        if name[:-1] == key[0]:
            return key
    return 0

def process_cap(residue):
    for atom in residue.atoms:
        if atom.mass == 0:
            key = get_key(residue.name,atom.name)
            if "HH3" in atom.name:
                key = get_key(residue.name,"H3C")
            if atom.name == "H":
                key = get_key(residue.name,"HN")
            if key == 0:
                atom.mass = 0
            else:
                atom.mass = param_dict[key]

if __name__ == "__main__":
    if len(sys.argv) != 4:
        sys.exit("Expected 3 arguments:  \n \n python pdbtinker.py <file.pdb> <parameters.prm> <output.xyz>")
    if len(sys.argv) == 4:
        if sys.argv[1].split(".")[-1] == "pdb" and sys.argv[2].split(".")[1] == "prm" and sys.argv[3].split(".")[1] == "xyz":
            system = load_pdb(sys.argv[1])
            process_new_parameters(sys.argv[2])
            process_PDB(system)
            build_xyz(system,sys.argv[3])
        else:
            sys.exit("Arguments invalid.\n\n\tpython pdbtinker.py <file.pdb> <parameters.prm> <output.xyz>")
