import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from numpy import linalg as LA


linklength = 1.6
DP = 20
polymeratomtype = [1,2,3]
file = "systemclustered.data"

def main():
    print("Reading")
    df, dim = read_lammpsBond(file)
    print("Isolating")
    polymer = isolate_polymer(df, polymeratomtype)
    print("FOFing")
    groupedP = fof(polymer, dim)
    groupedN = fofN(polymer, dim)
    itergroup(groupedP, groupedN, polymer, dim)

def itergroup(groupedP, groupedN, polymer, dim ):
    print(len(groupedP))
    #save(grouped)
    Naggr = []
    RadiusofGyration = []
    ROG2 = []
    if len(groupedP)>1: 
        for label, groups in enumerate(groupedP):
            trimmed = []
            for i in groupedN:
                for atom in groups:
                    if atom in i:
                        trimmed.append(groups)
                        break
            coords = transform(trimmed, polymer, dim)
            RadiusofGyration.append(ROG(coords))
            ROG2.append(list(gyrationTensor(coords)[1]))
            Naggr.append(len(coords))

            graph(coords, label+1)
    else:
        coords = transform(groupedN, polymer, dim)
        RadiusofGyration.append(ROG(coords))
        ROG2.append(gyrationTensor(coords)[1])

        Naggr.append(len(coords))

        graph(coords, 1)
    Naggr = np.array(Naggr)*(1/DP)
    RadiusofGyration = np.array(RadiusofGyration)
    print(Naggr)
    print(f"Naggr: {Naggr}")

    print(f"Average Naggr: {np.average(Naggr)}")
    print(f"Radius of Gyration: {RadiusofGyration}")
    print(f"Average Radius of Gyration: {np.average(RadiusofGyration)}")
    print(f'Principle Moments: {ROG2}')

def graph(coords,i):
    numpy_array = np.array(coords)
    transpose = numpy_array.T
    fig = plt.figure()
    ax = plt.axes(projection='3d')
    ax.set_title(f'group {i}')
    ax.scatter3D(transpose[0], transpose[1], transpose[2])
    plt.show()
def ROG(coords):
    numpy_array = np.array(coords)
    transpose = numpy_array.T
    xa = np.average(transpose[0])
    ya = np.average(transpose[1])
    za = np.average(transpose[2])
    ROG = 0
    for i in range(len(coords)):
        ROG +=   np.linalg.norm(np.array(coords[i]) - np.array([xa,ya, za]))**2
    return np.sqrt((1/len(coords))*ROG)
def gyrationTensor(coords):
    N =len(coords)
    print(f"atom number: {N}")
    numpy_array = np.array(coords)
    transpose = numpy_array.T   
    xa = np.average(transpose[0])
    ya = np.average(transpose[1])
    za = np.average(transpose[2])
    Rcm = np.array([xa, ya, za])
    
    S = np.zeros([3,3])
    for i in range(3):
        for j in range(3):
            S[i,j] = (1/N)*np.dot(transpose[i]-Rcm[i],transpose[j]-Rcm[j])
    w, v = LA.eig(S)
    return np.sqrt(np.sum(w)), w
    
            
def transform(groupedN, polymer, dim):
    maingroup = []
    length = 0 
    for groupcheck in groupedN:
        if len(groupcheck) > length:
            maingroup = groupcheck
        length = len(groupcheck)
    coords = []
    print("Main Group Found")
    for atom in maingroup:
        atomindex = polymer['AtomCount'].tolist().index(atom)
        coords.append(np.array([polymer['x'][atomindex], polymer['y'][atomindex], polymer['z'][atomindex]]))
    print("Iterating through SubGroups")    
    i = 0
 
    trimmedgroup = [i for i in groupedN if i != maingroup and i != []] 

    for group in trimmedgroup:
        print(f"\rProgress: {round(i/len(groupedN)*100,2)} %", end="")
        displacement = minimize(group,maingroup, polymer, dim)
        print("\nStart Displacement\n")
        for atom in group:
            atomindex = polymer['AtomCount'].tolist().index(atom)
            coords.append(np.array([polymer['x'][atomindex], polymer['y'][atomindex], polymer['z'][atomindex]])-displacement)
        i+=1
    return coords


def minimize(group,maingroup, polymer, dim):
    displacementofgroup = []
    for atom in maingroup:
        for atom2 in group:
            atomindex = polymer['AtomCount'].tolist().index(atom)
            atomindex2 = polymer['AtomCount'].tolist().index(atom2)
            distance, displacement = periodicdifference3d(np.array([polymer['x'][atomindex], polymer['y'][atomindex], polymer['z'][atomindex]]), np.array([polymer['x'][atomindex2], polymer['y'][atomindex2], polymer['z'][atomindex2]]), dim)
            displacementofgroup.append(displacement)

         

    return np.array(most_frequent(displacementofgroup))

def most_frequent(List):
    counter = 0
    num = List[0]
     
    for i in List:
        curr_frequency = List.count(i)
        if(curr_frequency> counter):
            counter = curr_frequency
            num = i
 
    return num

def save(grouped):
    textfile = open("fofresult.txt", "w")
    for element in grouped:
        textfile.write(element + "\n")
    textfile.close()

def read_lammpsBond(fname):
    df = pd.DataFrame(columns=['AtomCount', 'MoleculeCount', 'AtomType', 'x', 'y', 'z'])
    with open(fname) as f:
        line = f.readline()
        
        while line:
            if 'xhi' in line:
                dim = float(line.split()[1])
                break
            line = f.readline()

    with open(fname) as f:
        line = f.readline()
        while line:
            if len(line.split()) > 1 and line.split()[-1] == 'atoms':
                nparticles = int(line.split()[-2])
            elif len(line.split()) > 1 and line.split()[0] == 'Atoms':
                f.readline()
                for i in range(nparticles):
                    print(f"\rReading: {round(i/nparticles*100,2)} %", end="")
                    tag, mol, type, x, y, z, i , i , i  = map(float, f.readline().split())
     
                    df.loc[int(tag)-1] = [int(tag), int(mol), int(type), x, y, z] 
            line = f.readline()
    return df, dim 

def isolate_polymer(df, AtomType):
    polymer = pd.DataFrame(columns=['AtomCount', 'MoleculeCount', 'AtomType', 'x', 'y', 'z'])
    for index, row in df.iterrows():
        if int(row['AtomType']) in AtomType:
            polymer.loc[polymer.shape[0]] = row
    return polymer
def fof(polymer, dim):
    print()
    molecules = polymer['AtomCount'].tolist()
    init = len(molecules)
    grouping = []
    while len(molecules) > 0:

        currentatom = molecules[0]
        group = [currentatom]

        currentatomindex = polymer['AtomCount'].tolist().index(currentatom)
        molecules.remove(molecules[0])
        neigh = []
        for atom in molecules:
            atomindex = polymer['AtomCount'].tolist().index(atom)
            dist, _ = periodicdifference3d(np.array([polymer['x'][currentatomindex], polymer['y'][currentatomindex], polymer['z'][currentatomindex]]),  np.array([polymer['x'][atomindex], polymer['y'][atomindex], polymer['z'][atomindex]]), dim)
            if dist < linklength:
                neigh.append(atom)
                group.append(atom)
        while len(neigh) > 0:
            molecules = [i for i in molecules if i not in neigh] 
            print(f"\rPatching: {round((init-len(molecules))/init*100,2)} %", end="")

            NeighAtom = neigh[0]
            NeighAtomindex = polymer['AtomCount'].tolist().index(NeighAtom)
            for atom in molecules:
                atomindex = polymer['AtomCount'].tolist().index(atom)
                dist, _ = periodicdifference3d(np.array([polymer['x'][NeighAtomindex], polymer['y'][NeighAtomindex], polymer['z'][NeighAtomindex]]),  np.array([polymer['x'][atomindex], polymer['y'][atomindex], polymer['z'][atomindex]]), dim)
                if dist < linklength:
                    neigh.append(atom)
                    group.append(atom)
            neigh.remove(neigh[0])
        grouping.append(group)
    return grouping
def fofN(polymer, dim):
    print()
    molecules = polymer['AtomCount'].tolist()
    init = len(molecules)
    grouping = []
    while len(molecules) > 0:

        currentatom = molecules[0]
        group = [currentatom]

        currentatomindex = polymer['AtomCount'].tolist().index(currentatom)
        molecules.remove(molecules[0])
        neigh = []
        for atom in molecules:
            atomindex = polymer['AtomCount'].tolist().index(atom)
            dist = nonperiodicdifference(np.array([polymer['x'][currentatomindex], polymer['y'][currentatomindex], polymer['z'][currentatomindex]]),  np.array([polymer['x'][atomindex], polymer['y'][atomindex], polymer['z'][atomindex]]), dim)
            if dist < linklength:
                neigh.append(atom)
                group.append(atom)
        while len(neigh) > 0:
            molecules = [i for i in molecules if i not in neigh] 
            print(f"\rPatching: {round((init-len(molecules))/init*100,2)} %", end="")

            NeighAtom = neigh[0]
            NeighAtomindex = polymer['AtomCount'].tolist().index(NeighAtom)
            for atom in molecules:
                atomindex = polymer['AtomCount'].tolist().index(atom)
                dist = nonperiodicdifference(np.array([polymer['x'][NeighAtomindex], polymer['y'][NeighAtomindex], polymer['z'][NeighAtomindex]]),  np.array([polymer['x'][atomindex], polymer['y'][atomindex], polymer['z'][atomindex]]), dim)
                if dist < linklength:
                    neigh.append(atom)
                    group.append(atom)
            neigh.remove(neigh[0])
        grouping.append(group)
    return grouping
def nonperiodicdifference(p1, p2, dim):
    return  np.linalg.norm(p1 - p2)

def periodicdifference3d(p1, p2, dim ):
    variation = [-1, 0, 1]
    seperationVec = []
    displacementvec = []
    for i in variation:
        for j in variation:
            for k in variation:
                displacement = [i*dim,j*dim,k*dim]
                distance = np.linalg.norm(p1 - p2 + np.array(displacement))
                seperationVec.append(distance)
                displacementvec.append(displacement)
    minindex = seperationVec.index(min(seperationVec))


    return min(seperationVec), displacementvec[minindex]


main()