import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from numpy import linalg as LA
def main():
    text = 'sample.txt'
    data = acquire(text)
    
    plt.scatter(data[0], data[1])
    plt.show()
    
    
    plt.scatter(data[0], data[2])
    plt.show()
    plt.title('Micelle Formation Over Time')
    plt.scatter(data[0], data[3])
    plt.xlabel('Steps')
    plt.ylabel('Micelle Count')
    plt.show()
def acquire(text):
    data = []
    with open(text) as f:
        l = f.readline()
        while l:
        
            if len(l.split('\t')) > 1:
                print(l)
                print(l)
                point = [float(x) for x in l.split('\t')]
                data.append(np.array(point))
            l = f.readline()
    return np.transpose(data)
main()