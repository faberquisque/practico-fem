function u=vector_U(a,b,p_in,p_out,n)
K=matriz_K_total(a,b,n);
f=vector_ft(a,b,p_in,p_out,n);
u=K\f;
end
