import numpy as np

def Aij(df_i, df_j, c, f_i, f_j):
    
    # fun_podc = lambda x: -df_i(x)*df_j(x) + c*f_i(x)*f_j(x) 

    return lambda x: -df_i(x)*df_j(x) + c*f_i(x)*f_j(x)  # fun_podc
