function S=vector_Sig(E)
n=size(E,2);
S=zeros(2,n);
for i=1:n
S(:,i)=matriz_D*E(:,i);
end