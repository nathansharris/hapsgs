import numpy as np 
from itertools import combinations as comb
import sys

def isbreak(line):
    if 2 in line and 0 in line:
        return True
    else:
        return False

if __name__ == "__main__":
    regions = np.genfromtxt("regions.txt")[1:,1:]
    print(regions)

    full = np.zeros([len(regions),len(regions)])
    previous = "throw"

    for line in sys.stdin:
        line = ''.join(line)
        line = line\
            .replace("0|0","0")\
            .replace("0|1","1")\
            .replace("1|0","1")\
            .replace("1|1","2")
        line=np.array(line.split()).astype(int)
        print(isbreak(line))
