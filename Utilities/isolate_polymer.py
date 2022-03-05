import pandas as pd
import numpy as np
filename = 'systemAssembled.data'
filesave = 'systemAIsolated.data'
polymerAtoms = [1,2,3]
def main():
    df, bonds = read_lammpsBond(filename)
    polymer = isolate_polymer(df, polymerAtoms)

    write_lammpsBond(filename, polymer, bonds)   
    
def isolate_polymer(df, AtomType):
    polymer = pd.DataFrame(columns=['AtomCount', 'MoleculeCount', 'AtomType', 'x', 'y', 'z'])
    for index, row in df.iterrows():
        if int(row['AtomType']) in AtomType:
            polymer.loc[polymer.shape[0]] = row
    return polymer
def read_lammpsBond(fname):
    Atoms = pd.DataFrame(columns=['AtomCount', 'MoleculeCount', 'AtomType', 'x', 'y', 'z'])
    Bonds = pd.DataFrame(columns=['BondCount', 'BondType', 'AtomA', 'AtomB'])

    with open(fname) as f:
        line = f.readline()
        while line:
        
            if len(line.split()) > 1 and line.split()[-1] == 'atoms':
                nparticles = int(line.split()[-2])
            if len(line.split()) > 1 and line.split()[-1] == 'bonds':
                nbonds = int(line.split()[-2])
            elif len(line.split()) > 1 and line.split()[0] == 'Atoms':
                f.readline()
                for i in range(nparticles):
                    print(f"\rReading: {round(i/nparticles*100,2)} %", end="")
                    try:
                        tag, mol, type, x, y, z, i, i , i  = map(float, f.readline().split())
                    except:
                        tag, mol, type, x, y, z  = map(float, f.readline().split())
                        
     
                    Atoms.loc[int(tag)-1] = [int(tag), int(mol), int(type), x, y, z]
            elif'Bonds' in line:
                f.readline()
                for i in range(nbonds):
                    tag, type, a, b  = map(float, f.readline().split())
                        
     
                    Bonds.loc[int(tag)-1] = [int(tag), int(type), int(a), int(b)]         
            line = f.readline()
    return Atoms, Bonds





def write_lammpsBond(fname, Atoms, Bonds):
    print("writing")
    with open(fname) as f:
        line = f.readline()
        
        while line:
            if 'xhi' in line:
                dim = float(line.split()[1])
                break
            line = f.readline()

                
               
        
    with open(filesave, 'w') as out_file:
        out_file.write('LAMMPS Description  # full\n\n')

        out_file.write(f'\t{int(max((Atoms["AtomCount"])))} atoms\n')
        out_file.write(f'\t{int(max((Bonds["BondCount"])))} bonds\n')
        out_file.write(f'\t{int(max((Atoms["AtomType"])))} atom types\n')
        out_file.write(f'\t{7} bond types\n')
                    
        out_file.write(f'0.0 {dim} xlo xhi\n')
        out_file.write(f'0.0 {dim} ylo yhi\n')
        out_file.write(f'0.0 {dim} zlo zhi\n\n')


        
        out_file.write(f'Masses\n')
        if int(max((Atoms["AtomType"]))) == 2:
            out_file.write(f"\n1 1.0  # HR\n2 1.0  # PR\n\n")
        elif int(max((Atoms["AtomType"]))) == 3:
            out_file.write(f"\n1 1.0  # HR\n2 1.0  # PR\n3 1.0  # SR\n\n")

        
        
        out_file.write('Atoms  # bond\n\n')

        for index, row in Atoms.iterrows():
            out_file.write(f"{int(row['AtomCount'])}\t{int(row['MoleculeCount'])}\t{int(row['AtomType'])}\t{row['x']}\t{row['y']}\t{row['z']}\n")
        out_file.write('\nBonds\n\n')

        for index, row in Bonds.iterrows():
            out_file.write(f"{int(row['BondCount'])}\t{int(row['BondType'])}\t{int(row['AtomA'])}\t{row['AtomB']}\n")
        

        out_file.close()
        
main()