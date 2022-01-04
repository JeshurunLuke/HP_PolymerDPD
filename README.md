# HP_PolymerDPD

**Folders**:

**bond**: 

      Summary: HP polymer where polymerize.py can be used to specify boxSize, H and P monomer count, polymer density and relative bead density
      Modification: 
        polymerize.py: (rc, p, bxSize, processes, H, P, polymerPercentage, bondlen)
                       (water bead distance, bead density, number of parallel process, H count, P count, polymer bead %, polymer bead distance)
      Requires: moltemplate.sh
      Output: systemSolvated.data (polymer and water bonds and locations)
      Subfolder: RunFolder
        Summary: Takes systemSolvated.data to run a DPD simulation ("lmp_mpi -in equilibrate2.in")
        Modification: Can Modify system.in.init as needed to change pair potential or equilibrate2.in to change run parameters
        Output: system.data, dump.soft, log.lammps
       
 
**bond3bead**: 

      **Summary**: HP-S polymer performs similarily to the bond folder but with the addition of the monomer S (reactive polar with current pair potential).
      Modification: 
        polymerize.py: (rc, p, bxSize, processes, H, P, S, polymerPercentage, bondlen)
                       (water bead distance, bead density, number of parallel process, H count, P count, S count, polymer bead %, polymer bead distance)
      Requires: moltemplate.sh
      Output: systemSolvated.data 
      Subfolder: RunFolder
        Summary: Takes systemSolvated.data to run a DPD simulation ("lmp_mpi -in equilibrate2.in")
        Modification: Can Modify system.in.init as needed to change pair potential or equilibrate2.in to change run parameters
        Output: system.data, dump.soft, log.lammps
        
        
 **React**: 
 
      **Summary**: Replaces random water molecules and updates atom and bond type count (Requires system.data)
      Modification:
        replace.py: (AtomType, amtOfreactant)
                    (Water Molecule Type, Amt of Reactant)
      Requires: none
      Output: systembond.data
      Subfolder: Reaction
        Summary: Takes system.data (renamed systembond.data) and iteratively runs a reaction and creates bonds if bond criteria are met
         Modification: Can Modify pair potentials in system.in.settings or runtime in equilibrate2.in
          react.py: (associationdistance, iterations)  
                    (bond distance requirement, run count)
   
**FoFQuantify**

    **Summary**: Reads a data file and groups atoms according to Friends of Friends algorithm 
    Modification:
      fof.py: (DP, linklength,polymeratomtype, file)
              (DP of polymer, linking length of FOF algorithm, atomtype of polymer from atomtype, filename of datafile)
    Requires: none
    Output: plots of the position of all points, Naggr,  Rg of atoms, and the eigenvalues of the gyration tensor

**Utilities**

      **Summary**: Creates datafile of just the polymer from a solvated file
      Modification
        isolate_polymer.py: (filename, filesave, polymerAtoms)
                            (file name of solvated file, filename you want to save to, atomtype of the polymer atoms)
      Requires: none
      Output: data file with the filesave name
      
**Notes**

All equilibrate2.in can be run in parallel (mpirun -np 4) and (bond3bead and bond) requires the subprocess and multiprocessing python modules

**Credits**


Solvating Matrix Method: Idea for grid type method from pysimm

ReadandWrite Functions: Modified versions from pysimm system

Interaction Parameters: Pair Potentials from "https://pubs.rsc.org/en/content/articlehtml/2018/ra/c8ra07023g#imgfig6"

DPD Paramters:  "https://pubs.rsc.org/en/content/articlehtml/2018/ra/c8ra07023g#imgfig6" and https://www.mdpi.com/2073-4360/13/13/2193

Initial Polymer Config: moltemplate
