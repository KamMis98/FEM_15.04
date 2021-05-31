import numpy as np
import matplotlib.pyplot as plt
from RysujGeometrie import *

def RysRozw(NODES, ELEMS, WB, u):
    
    RysGeo(NODES,ELEMS,WB)
    
    x = NODES[:,1]
    
    plt.plot(x, u, 'm*')