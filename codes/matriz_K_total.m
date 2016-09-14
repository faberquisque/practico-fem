function K=matriz_K_total(a,b,n)
m=linspace(a,b,n+1);
K=zeros(n+1,n+1);
for i=1:n
    Ka=matriz_K(m(i),m(i+1));
    K([i,i+1],[i,i+1])=K([i,i+1],[i,i+1])+Ka;
end
end

    
