"""
pdbtinker.py
Converts PDBs to Tinker XYZs

Handles the primary functions
"""
def disambiguate_histidines(temp):
    """
    Returns:
        None
    Description:
        Processes through the PDB by residue, searching for ambiguous histidines,
        then determines the protonation state and reassigns the residue name
        appropriate to that state.  This ensures proper histidine processing later.
    """
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

def check_terminals(temp):
    """
    Returns:
        None
    Description:
        Processes PDB system by residue and checks for known N-terminal and
        C-terminal specific atoms in amino acid residues, and known 5'/3' atoms
        in nucleic acids. When detected, it modifies the residue name to include
        the appropriate terminal indicator.
    """
    for residue in temp.residues:
        if residue.name in amino_acid_list:
            atomlist = []
            for atom in residue.atoms:
                atomlist.append(atom.name)
            if 'OXT' in atomlist:
                residue.name = "C" + residue.name
            elif 'H3' in atomlist:
                residue.name = "N" + residue.name
        elif residue.name in nucleic_acid_list:
            atomlist = []
            for atom in residue.atoms:
                atomlist.append(atom.name)
            if "HO5'" in atomlist:
                residue.name = residue.name + "5"
            elif "HO3'" in atomlist:
                residue.name = residue.name + "3"
            if "D" not in residue.name and "R" not in residue.name:
                residue.name = "R" + residue.name

def build_dictionary(temp,temp_param_file):
    """
    Returns:
        Dictionary composed of atom types found in the loaded PDB and with values given by the parameter file.
    """
    ### Build dictionary for all PDB atoms found in loaded file.
    ### This sets the atom type for every atom to 0 initially.
    atom_type_dictionary={}

    for atom in temp.atoms:
        res_name = atom.residue.name
        atom_name = atom.name
        if res_name not in atom_type_dictionary:
            atom_type_dictionary[res_name] = {}
        if atom_name not in atom_type_dictionary[res_name]:
            atom_type_dictionary[res_name][atom_name] = 0


    ### Load the key converter.  This will ensure that the parameter file names for residues are properly translated
    ### This may need to be updated in the future/as the project undergoes more testing.
    key_convert=np.load("dictionaries/key_convert_dict.npy").item()
    params = temp_param_file.readlines()
    temp_param_file.close() ### Always close what you've opened.

    ### Get biotype data from parameter file.  May also need updating later, depending on future parameter sets.
    biotypes = []
    for line in params:
        if "biotype" in line:
            biotypes.append(line)

    ### Build dictionary from parameter file.
    param_set = {}
    for line in biotypes:
        old_key = line.split("\"")[1]
        atom_type = line.split("\"")[2].split()[0]
        atom_name = line.split("\"")[0].split()[-1].replace("*","\'")
        new_key = key_convert[old_key]
        if new_key not in param_set.keys():
            param_set[new_key] = {}
        if atom_name not in param_set[new_key].keys():
            param_set[new_key][atom_name] = atom_type
    ### In case the parameters don't include HIS, I've assumed all HIS residues are HID.  If not, fix your PDB.
    ### HIS is ambiguous, and there is no place for ambiguity here.
    # param_set["HIS"] = param_set["HID"]

    ### Fill the PDB Dictionary set with corresponding atom types from the parameter set.
    ### Atom types that aren't in the parameter file will remain 0 and get flagged on file output.
    ### Don't look at this loop, it's chaotic and inelegant.  But it works so far, at least for amino acids.
    for res in atom_type_dictionary.keys():
        for atom in atom_type_dictionary[res].keys():
            if res in param_set.keys() and atom_type_dictionary[res][atom] == 0:
                if atom in param_set[res].keys():
                    atom_type_dictionary[res][atom] = param_set[res][atom]
                elif ''.join(i for i in atom if not i.isdigit()) in param_set[res].keys():
                    atom_type_dictionary[res][atom] = param_set[res][''.join(i for i in atom if not i.isdigit())]
                elif atom[:-1] in param_set[res].keys():
                    atom_type_dictionary[res][atom] = param_set[res][atom[:-1]]
                elif atom[:-2] in param_set[res].keys():
                    atom_type_dictionary[res][atom] = param_set[res][atom[:-2]]
                elif atom in ["H","H1","H2","H3"] and "HN" in param_set[res].keys():
                    atom_type_dictionary[res][atom] = param_set[res]["HN"]
            if res[1:] in param_set.keys() and atom_type_dictionary[res][atom] == 0:
                if atom in param_set[res[1:]].keys():
                    atom_type_dictionary[res][atom] = param_set[res[1:]][atom]
                elif ''.join(i for i in atom if not i.isdigit()) in param_set[res[1:]].keys():
                    atom_type_dictionary[res][atom] = param_set[res[1:]][''.join(i for i in atom if not i.isdigit())]
                elif atom[:-1] in param_set[res[1:]].keys():
                    atom_type_dictionary[res][atom] = param_set[res[1:]][atom[:-1]]
                elif atom[:-2] in param_set[res[1:]].keys():
                    atom_type_dictionary[res][atom] = param_set[res[1:]][atom[:-2]]
                elif atom in ["H","H1","H2","H3"]:
                    atom_type_dictionary[res][atom] = param_set[res[1:]]["HN"]
    return atom_type_dictionary


def print_tinker_xyz(temp,atom_type_dictionary):
    """
    Returns:
        None
    Description:
        Prints system as Tinker XYZ to the terminal/stdout.
    """
    print(str(len(temp.atoms)))
    for atom in temp.atoms:
        atomstring = str(atom.idx+1) + "\t" + str(atom.name) + "\t" + str(atom.xx) + "\t" + str(atom.xy) + "\t" + str(atom.xz) + "\t" + str(atom_type_dictionary[atom.residue.name][atom.name]) + "\t"
        for i in atom.bonds:
            if i.atom1.idx == atom.idx:
                atomstring = atomstring + str(i.atom2.idx+1) + "\t"
            elif i.atom2.idx == atom.idx:
                atomstring = atomstring + str(i.atom1.idx+1) + "\t"
        if atom_type_dictionary[atom.residue.name][atom.name] == 0:
            atomstring = atomstring + "ATOM TYPE NOT FOUND\t"
        print(atomstring)

def make_tinker_xyz_file(temp, atom_type_dictionary, filename):
    """
    Returns:
        None
    Description:
        Outputs Tinker XYZ file to designated filename.
    """
    f = open(filename,"w+")
    f.write(str(len(temp.atoms))+"\n")
    for atom in temp.atoms:
        atomstring = str(atom.idx+1) + "\t" + str(atom.name) + "\t" + str(atom.xx) + "\t" + str(atom.xy) + "\t" + str(atom.xz) + "\t" + str(atom_type_dictionary[atom.residue.name][atom.name]) + "\t"
        for i in atom.bonds:
            if i.atom1.idx == atom.idx:
                atomstring = atomstring + str(i.atom2.idx+1) + "\t"
            elif i.atom2.idx == atom.idx:
                atomstring = atomstring + str(i.atom1.idx+1) + "\t"
        if atom_type_dictionary[atom.residue.name][atom.name] == 0:
            atomstring = atomstring + "ATOM TYPE NOT FOUND\t"
        f.write(atomstring+"\n")


if __name__ == "__main__":
    # Do something if this file is invoked on its own
    print("Hello")
