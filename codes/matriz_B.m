function B=matriz_B(XI,a,b)
B=[-1/(b-a),1/(b-a);0.5*(1-XI)/radio(XI,a,b),0.5*(1+XI)/radio(XI,a,b)];
end