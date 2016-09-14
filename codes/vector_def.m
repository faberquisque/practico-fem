function E=vector_def(u,a,b)
n=length(u)-1;
m=linspace(a,b,n+1);
E=zeros(2,n);
for i=1:n
E(:,i)=matriz_B(0.5,m(i),m(i+1))*u([i,i+1]);
end
