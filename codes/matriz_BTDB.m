function BTDB=matriz_BTDB(XI,a,b)
BTDB=matriz_B(XI,a,b)'*matriz_D*matriz_B(XI,a,b)*radio(XI,a,b);
end
