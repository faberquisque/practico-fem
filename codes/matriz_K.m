function K=matriz_K(a,b)
K=pi*(b-a)*gauss_int(@(XI)matriz_BTDB(XI,a,b),4);
