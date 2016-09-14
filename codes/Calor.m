function [T_g,x]=Calor(n)
a=0.1;
T=[0;0;0;10];
T_inf=[0;20;0;0];
q=[0;-100;0;0];
h=[0;10;0;0];
k=14.9;
borde=[true;true;true;false];
n_nodos=(n+1)^2;        % n_elementos=2*n^2
K_g=zeros(n_nodos);
T_g=zeros(n_nodos,1);
f_g=zeros(n_nodos,1);
flag=true(n_nodos,1);
Q_vol=10000;

for i=1:n
    for j=1:n
        n_i=i+(n+1)*(j-1);
        n_j=(i+1)+(n+1)*(j-1);
        n_k=(i+1)+(n+1)*j;
        n_l=i+(n+1)*j;
        x_i=[(i-1)*a/n;(j-1)*a/n;0]';
        x_j=[i*a/n;(j-1)*a/n;0]';
        x_k=[i*a/n;j*a/n;0]';
        x_l=[(i-1)*a/n;j*a/n;0]';
        
        %elemento 1
        L=[modulo(x_i-x_j);modulo(x_j-x_k);modulo(x_k-x_i)];
        A=[x_i;x_j;x_k];
        A=circshift(A,1)-circshift(A,2);
        b=-A(:,2);
        c=A(:,1);
        Area=modulo(cross(x_j-x_i,x_k-x_i))/2;
        K_a=1/4/Area*k*((b*b')+(c*c'));
        tipo_2=false(3,1);
        h1=zeros(3,1);
        q1=h1;
        tinf1=h1;
        if j==1
            if ~borde(1)
                flag(n_i)=borde(1);
                T_g(n_i)=T(1);
                flag(n_j)=borde(1);
                T_g(n_j)=T(1);
            end
            tipo_2(1)=borde(1);
            h1(1)=h(1);
            q1(1)=q(1);
            tinf1(1)=T_inf(1);
        end
        if i==n
            if ~borde(2)
                flag(n_j)=borde(2);
                T_g(n_j)=T(2);
                flag(n_k)=borde(2);
                T_g(n_k)=T(2);
            end
            tipo_2(2)=borde(2);
            h1(2)=h(2);
            q1(2)=q(2);
            tinf1(2)=T_inf(2);
        end
        aux=h1.*L/6.*(tipo_2);
        K_aux=[0,aux(1),aux(3);0,0,aux(2);0,0,0];
        K_a=K_a+K_aux+K_aux'+2*diag(aux)+2*diag(circshift(aux,1));
        f_Q=Q_vol*Area/3*ones(3,1);
        f_c=-0.5*L.*(q1-h1.*tinf1).*(tipo_2);
        f_Q=f_Q+f_c+circshift(f_c,1);
        nodo=[n_i;n_j;n_k];
        K_g(nodo,nodo)=K_g(nodo,nodo)+K_a;
        f_g(nodo)=f_g(nodo)+f_Q;
        %elemento 2
        L=[modulo(x_i-x_l);modulo(x_l-x_k);modulo(x_k-x_i)];
        A=[x_i;x_l;x_k];
        A=circshift(A,1)-circshift(A,2);
        b=-A(:,2);
        c=A(:,1);
        Area=modulo(cross(x_l-x_i,x_k-x_i))/2;
        K_a=1/4/Area*k*((b*b')+(c*c'));
        tipo_2=false(3,1);
        h1=zeros(3,1);
        q1=h1;
        tinf1=h1;
        if i==1
            if ~borde(4)
                flag(n_i)=borde(4);
                T_g(n_i)=T(4);
                flag(n_l)=borde(4);
                T_g(n_l)=T(4);
            end
            tipo_2(1)=borde(4);
            h1(1)=h(4);
            q1(1)=q(4);
            tinf1(1)=T_inf(1);
        end
        if j==n
            if ~borde(3)
                flag(n_l)=(borde(3));
                T_g(n_l)=T(3);
                flag(n_k)=(borde(3));
                T_g(n_k)=T(3);
            end
            tipo_2(2)=(borde(3));
            h1(2)=h(3);
            q1(2)=q(3);
            tinf1(2)=T_inf(3);
        end
        aux=h1.*L/6.*(tipo_2);
        K_aux=[0,aux(1),aux(3);0,0,aux(2);0,0,0];
        K_a=K_a+K_aux+K_aux'+2*diag(aux)+2*diag(circshift(aux,1));
        f_Q=Q_vol*Area/3*ones(3,1);
        f_c=-L/2.*(q1-h1.*tinf1).*(tipo_2);
        f_Q=f_Q+f_c+circshift(f_c,1);
        
        nodo=[n_i;n_l;n_k];
        
        K_g(nodo,nodo)=K_g(nodo,nodo)+K_a;
        f_g(nodo)=f_g(nodo)+f_Q;
    end
end
%
T_g(flag)=K_g(flag,flag)\(f_g(flag)-K_g(flag,~flag)*T_g(~flag));
T_g=rot90(reshape(T_g,n+1,n+1));
x=linspace(0,a,n+1);