import numpy as np

def Aij(df_i, df_j, c, f_i, f_j):
    '''
    Parametry :
    df_i , df_j - pochodne funkcji ksztaltu zwiazanych z i-tym oraz j- tym wezlem
    Zwraca :
    wyrazenie lambda bedace funkcja podcalkowa A_ij
    '''
    
    # fun_podc = lambda x: -df_i(x)*df_j(x) + c*f_i(x)*f_j(x) 

    return lambda x: -df_i(x)*df_j(x) + c*f_i(x)*f_j(x)  # fun_podc
