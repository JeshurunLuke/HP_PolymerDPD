import numpy as np
import random as rd
import pandas as pd
from subprocess import call, Popen, PIPE
from multiprocessing import Process
from itertools import repeat
import time
from multiprocessing import Pool, freeze_support
import os
#Bond length of 1 (If edited modify forcefield.lt as well)
bondlen = [0.5, 0.5]
rc = 1.3
bxSize = 30.0
p = 3
processes = 8

from math import log10, floor
def round_sig(x, sig=2):
    return round(x, sig-int(floor(log10(abs(x))))-1)



def rmv():
    p2 = Popen(['rm','-r',  'output_file.lt'], stdout=PIPE)
    p2 = Popen(['rm','-r',  'tmpbond.txt'], stdout=PIPE)
    p2 = Popen(['rm','-r',  'tmpmon.txt'], stdout=PIPE)
    p2 = Popen(['rm','-r',  'system.lt'], stdout=PIPE)
    p2 = Popen(['rm','-r',  'system.in'], stdout=PIPE)


    
def boxinit(H,P,S, NP):
    mxd = 360
    entry = 19
    len = bxSize
    print(f"Box paramater: {len}")
    with open("system.lt", "w") as f:
        box = len
        f.write('import "output_file.lt\"\n')
        f.write('write_once("Data Boundary") {\n')
        f.write(f'\t0 {box} xlo xhi \n\t0 {box}  ylo yhi\n\t0 {box}  zlo zhi\n ')
        f.write('}\n')
        f.write("polymers = new random([Polymer,\n")
        for i in range(entry):
            rotation = [rd.randrange(0, mxd, 1), rd.randrange(0, mxd, 1), rd.randrange(0, mxd, 3)]
            if i == entry-1:
                f.write(f"                     Polymer.rot({rotation[0]},{rotation[1]},{rotation[2]},1)],\n")
            else:
                f.write(f"                     Polymer.rot({rotation[0]},{rotation[1]},{rotation[2]},1),\n")
        probvec = stringen(np.ones(entry+1)*(1/(entry+1)))
        f.write(f'                     {probvec},\n')
        f.write("                     123456)\n")
        f.write(f"                     [{NP}].move({len/NP}, 0, 0) \n                     [{NP}].move(0, {len/NP} , 0) \n                     [{NP}].move(0, 0, {len/NP})")    
    run()
def stringen(vec):
    string = ""
    for i in range(0,len(vec)):
        if i == 0:
            string = f"[{vec[0]}, "
        elif i == len(vec)-1:
            string += f" {vec[i]}]"
        else:
            string += f" {vec[i]},"
    return string
        

def run():
    p = Popen('moltemplate.sh system.lt -atomstyle full', shell=True)  
    p.wait()
    #equilibrate()
    

def minimize():
    p = Popen(['lmp_mpi', '-in', 'equilibrate2.in'])
    p.wait()
    print("Done Equilibrating")
def polymerize(H,P,S, order):
    polist = []    
    for section in order:
        if section == 'H':
            for i in range(0, H):
                polist.append('H')
        if section == 'P':
            for i in range(0, P):
                polist.append('P')
        if section == 'S':
            for i in range(0, S):
                polist.append('S')

    monlist(polist)
    bondlist(polist)
    filenames = ["polback.txt", "tmpmon.txt", "tmpbond.txt"]

    with open("output_file.lt", "w") as outfile:
        for filename in filenames:
            with open(filename) as infile:
                outfile.write(infile.read())
        outfile.write('}')
        outfile.close()

def monlist(polist):
    rot = 180

    with open('tmpmon.txt', 'w') as f:
        for i,mon in enumerate(polist):
            if i < len(polist)-1:
                if polist[i + 1] == 'H':
                    if i%2 == 0:
                        f.write(f'  mon{i} = new {mon}.move({i*bondlen[0]},0,0)\n')
                    else:
                        f.write(f'  mon{i} = new {mon}.rot({i%2*rot}, {i%2*1},0,0).move({i*bondlen[0]},0,0)\n')
                elif polist[i-1] == 'S':
                    if i%2 == 0:
                        f.write(f'  mon{i} = new {mon}.move({i*bondlen[0]},0,0)\n')
                    else:
                        f.write(f'  mon{i} = new {mon}.rot({i%2*rot}, {i%2*1},0,0).move({i*bondlen[0]},0,0)\n')        
                else:
                    if i%2 == 0:
                        f.write(f'  mon{i} = new {mon}.move({i*bondlen[0]},0,0)\n')
                    else:
                        f.write(f'  mon{i} = new {mon}.rot({i%2*rot}, {i%2*1},0,0).move({i*bondlen[0]},0,0)\n')
            else:
                if polist[i - 1] == 'H':
                    if i%2 == 0:
                        f.write(f'  mon{i} = new {mon}.move({i*bondlen[0]},0,0)\n')
                    else:
                        f.write(f'  mon{i} = new {mon}.rot({i%2*rot}, {i%2*1},0,0).move({i*bondlen[0]},0,0)\n')
                elif polist[i -1] == 'S':
                    if i%2 == 0:
                        f.write(f'  mon{i} = new {mon}.move({i*bondlen[0]},0,0)\n')
                    else:
                        f.write(f'  mon{i} = new {mon}.rot({i%2*rot}, {i%2*1},0,0).move({i*bondlen[0]},0,0)\n')
                else:
                    if i%2 == 0:
                        f.write(f'  mon{i} = new {mon}.move({i*bondlen[0]},0,0)\n')
                    else:
                        f.write(f'  mon{i} = new {mon}.rot({i%2*rot}, {i%2*1},0,0).move({i*bondlen[0]},0,0)\n')
        f.close()

def bondlist(polist):
    with open('tmpbond.txt', 'w') as f:
        f.write('  write("Data Bond List") {\n')
        for i in range(1, len(polist)-1):
            f.write(f'  \t$bond:backbone{i}\t$atom:mon{i}/ca\t$atom:mon{i+1}/ca\n')
        f.write('  }\n')
        f.close()


def read_lammps(fname):
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
                        tag, mol, type, x, y, z , i, i , i = map(float, f.readline().split())
                    except:
                        tag, mol, type, x, y, z  = map(float, f.readline().split())
     
                    df.loc[int(tag)-1] = [int(tag), int(mol), int(type), x, y, z] 
            line = f.readline()
    return df


def solvate(nonSolvated, waterBeads):
    print()

    #mass = 1 #Mass Arbritary?

    ngrid = int(waterBeads**(1/3))
    #ngrid = np.floor(((bxSize ** 3) * 0.6022 / mass) ** (1.0 / 3.0)) #Mass Paramter
    print(f" {ngrid**3}: approximate")
    # put water molecules in the nodes of a regular grid

    S = 2

    x = np.linspace(0, bxSize, int(ngrid) + 1)
    y = np.linspace(0, bxSize, int(ngrid) + 1)
    z = np.linspace(0, bxSize, int(ngrid) + 1)
    process = []
    for i in range(0, S):
        for j in range(0, S):
            for k in range(0, S):
                process.append([np.array_split(x, S)[i], np.array_split(y, S)[j], np.array_split(z, S)[k]])
    
            
    return nonSolvated, process
def search(nonSolvated, rng):
    print()
    solvnt = pd.DataFrame(columns=['AtomCount', 'MoleculeCount', 'AtomType', 'x', 'y', 'z'])

    count = 0 
    for p in rng[0]:
        for q in rng[1]:
            for t in rng[2]:
                flags = []
                print(f"\rSolvating: {round(count/(len(rng[0])*len(rng[1])*len(rng[2]))*100,2)} %", end="")
                for index, row in nonSolvated.iterrows():
                    dist = np.linalg.norm(np.array([row['x'], row['y'], row['z']]) - np.array([p, q, t]))
                    if dist > rc:
                        flags.append(dist > rc)
                    else:
                        break
                if all(flags):
                    count += 1
                    tag, mol, type, charge, x, y, z = count, count, 4, 0, p,q,t
                    solvnt.loc[count-1] =  [int(tag), int(mol), int(type), x, y, z] 
    print()

    return solvnt
def write_lammps(fname, df):
    print()
    with open(fname) as f:
        line = f.readline()
        while line:
            if len(line.split()) > 1 and line.split()[-1] == 'bonds':
                bonds = int(line.split()[-2])
            if 'atom types' in line:
                atomtype = int(line.split()[-3])
            if 'bond types' in line:
                bondtype = int(line.split()[-3])
            line = f.readline()

    with open("systemSolvated.data", 'w') as out_file:
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
                    atomtype = int(line.split()[-3])
                    out_file.write(f'\t{atomtype+1} atom types\n')

                if 'bond types' in line:
                    bondtype = int(line.split()[-3])
                    out_file.write(f'\t{bondtype} bond types\n')
                    
                    while line:
                        line = f.readline()
                        if 'Masses' not in line and len(line)>0:
                            out_file.write(line)
                        else:
                            out_file.write(line)
                            out_file.write(f"\n1 1.0  # HR\n2 1.0  # PR\n3 1.0  # SR\n4 1.0  # W\n\n")
                            break
                line = f.readline()

        
        
        out_file.write('Atoms  # full\n\n')

        for index, row in df.iterrows():
            out_file.write(f"{int(row['AtomCount'])}\t{int(row['MoleculeCount'])}\t{int(row['AtomType'])}\t{row['x']}\t{row['y']}\t{row['z']}\n")
        out_file.write('\nBonds\n\n')
        
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
                    tag, mol, type, charge, x, y, z  = map(float, f.readline().split())
     
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
                atomtype = int(line.split()[-3])
            if 'bond types' in line:
                bondtype = int(line.split()[-3])
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
                    atomtype = int(line.split()[-3])
                    out_file.write(f'\t{atomtype} atom types\n')

                if 'bond types' in line:
                    bondtype = int(line.split()[-3])
                    out_file.write(f'\t{bondtype} bond types\n')
                    
                    while line:
                        line = f.readline()
                        if 'Masses' not in line and len(line)>0:
                            out_file.write(line)
                        else:
                            out_file.write(line)
                            out_file.write(f"\n1 1.0  # HR\n2 1.0  # PR\n3 1.0  # SR\n\n")
                            break
                line = f.readline()

        
        
        out_file.write('Atoms  # full\n\n')

        for index, row in df.iterrows():
            out_file.write(f"{int(row['AtomCount'])}\t{int(row['MoleculeCount'])}\t{int(row['AtomType'])}\t{row['x']}\t{row['y']}\t{row['z']}\n")
        out_file.write('\nBonds\n\n')
        
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
        

if __name__ == "__main__":
    #Number of H and P monomers
    H = 9
    P = 12
    S = 6
    polymerPercentage = 0.05
    order = ['H', 'P', 'S']

    TotalBeads =  (bxSize**3)*p
    polymerBeads = TotalBeads*polymerPercentage

    #Number of polymer (actual amt = 5**3)
    NumofPol = int((polymerBeads/(H+P+S))**(1/3))
    print(f"Actual Percentage: {NumofPol**3/TotalBeads *100}")
    polymerize(H, P, S, order)
    boxinit(H, P,S,  NumofPol)
    df = read_lammpsBond("system.data")
    write_lammpsBond("system.data", df)
    minimize()
    rmv()

    waterBeads = TotalBeads - NumofPol**3

    df = read_lammps("system.data")
    nonSolvated, process = solvate(df,waterBeads)
    print(len(process))
    with Pool(processes) as p:
        J = p.starmap(search, zip( repeat(nonSolvated), process))
    Solvated = pd.DataFrame()
    Solvated = Solvated.append(df)
    for i in J:
        MaxAtom = max(Solvated["AtomCount"])
        MaxMol = max(Solvated["MoleculeCount"])
        i['AtomCount'] += MaxAtom
        i['MoleculeCount'] += MaxMol
        print(len(i))
        Solvated = Solvated.append(i)
                            
    write_lammps("system.data", Solvated)