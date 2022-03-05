import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

def main():
    systemVar = read_log('log.lammps')
    plt.title('Energy Fluctuation')
    plt.plot(systemVar['Step'][-100:], systemVar['TotalEng'][-100:])
    plt.xlabel('Steps')
    plt.ylabel('Energy')

    plt.show()
    

def read_log(fname):

               
    SystemVar = pd.DataFrame(columns=['Step', 'Temp', 'PotEng', 'TotalEng', 'E_pair', 'E_bond'])
      
    with open(fname) as f:
        line = f.readline()
        record = False
        count = 0 
        while line:
            if len(line.split()) > 1 and line.split()[0] == 'Loop':
                record = False
            if record:
                Step, Temp, PotEng, TotEng, Press, Volume, E_pair, E_bond, E_angle, E_dihed  = map(float, line.split())
                SystemVar.loc[count] = [Step, Temp, PotEng, TotEng, E_pair, E_bond]
                count +=1
                    
            if len(line.split()) > 1 and line.split()[0] == 'Step':
                record = True


            line = f.readline()

    return SystemVar


main()