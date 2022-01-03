import pandas as pd
import numpy as np
# Replaces random water molecules and updates atom and bond type count (Requires system.data)
def main():
    AtomType = [4]

    df = read_lammpsBond("system.data") 
    system = isolate_water(df, AtomType)
    write_lammpsBond("system.data", system)
def isolate_water(df, AtomType):
    
    count = 0 
    for index, row in df.iterrows():
        if int(row['AtomType']) in AtomType:
            if np.random.uniform(0,1)> 0.9:
                df.iloc[index]['AtomType'] = df.iloc[index]['AtomType'] + 1
                count +=1

        if count == 50:
            break
    return df
def read_lammpsBond(fname):
    df = pd.DataFrame(columns=['AtomCount', 'MoleculeCount', 'AtomType', 'x', 'y', 'z'])
    with open(fname) as f:
        line = f.readline()
        while line:
            if len(line.split()) > 1 and line.split()[-1] == 'atoms':
                nparticles = int(line.split()[-2])
            elif len(line.split()) > 1 and line.split()[0] == 'Atoms':
                f.readline()
                for i in range(nparticles):
                    print(f"\rReading: {round(i/nparticles*100,2)} %", end="")
                    try:
                        tag, mol, type, x, y, z, i, i , i  = map(float, f.readline().split())
                    except:
                        tag, mol, type, x, y, z  = map(float, f.readline().split())
                        
     
                    df.loc[int(tag)-1] = [int(tag), int(mol), int(type), x, y, z] 
            line = f.readline()
    return df





def write_lammpsBond(fname, df):
    print()
    with open(fname) as f:
        line = f.readline()
        while line:
            if len(line.split()) > 1 and line.split()[-1] == 'bonds':
                bonds = int(line.split()[-2])
            if 'atom types' in line:
                atomtype = int(line.split()[-3]) + 1
            if 'bond types' in line:
                bondtype = int(line.split()[-3]) + 1
            line = f.readline()

    with open("systembond.data", 'w') as out_file:
        out_file.write('LAMMPS Description  # full\n\n')

        out_file.write(f'\t{int(max((df["AtomCount"])))} atoms\n')
        with open(fname) as f:
            line = f.readline()
            while line:
                if len(line.split()) > 1 and line.split()[-1] == 'bonds':
                    bonds = int(line.split()[-2])
                    out_file.write(f'\t{bonds} bonds\n')
                    #out_file.write(f'\t0 angles\n\t0  dihedrals\n\t0  impropers\n\n')
                if 'atom types' in line:
                    out_file.write(f'\t{atomtype} atom types\n')

                if 'bond types' in line:
                    out_file.write(f'\t{bondtype} bond types\n')
                    
                    while line:
                        line = f.readline()
                        if 'Masses' not in line and len(line)>0:
                            out_file.write(line)
                        else:
                            out_file.write(line)
                            out_file.write(f"\n1 1.0  # HR\n2 1.0  # PR\n3 1.0  # SR\n4 1.0  # W\n5 1.0  # M\n\n")
                            break
                line = f.readline()

        
        
        out_file.write('Atoms  # full\n\n')

        for index, row in df.iterrows():
            out_file.write(f"{int(row['AtomCount'])}\t{int(row['MoleculeCount'])}\t{int(row['AtomType'])}\t{row['x']}\t{row['y']}\t{row['z']}\n")
        out_file.write('\nBonds #iteration\n\n')
        
        with open(fname) as f:
            line = f.readline()
            while line:
                if 'Bonds' in line:
                    f.readline()
                    for i in range(bonds):
                        print(f"\rWriting Bonds: {round(i/bonds*100,2)} %", end="")
                        out_file.write(f.readline())
                line = f.readline()
        
        out_file.close()
        
main()