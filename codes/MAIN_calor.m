%% Parametros iniciales
a=0.1;
q=-100;
h=10;
T0=10;
Tinf=20;
q_vol=10000;
k=14.9;
max=6;
%constante A de la solucion analitica
A=(q_vol*a*(1+h/k*a/2)-q+h*(Tinf-T0))/(h*a+k);
E=zeros(1,max);
N=zeros(1,max);
for i=1:max
    %parametro n
    n=2^i;
    %calculo numerico del campo de T
    [T2,x2]=Calor(n);
    %calculo analitico del campo de T
    T2t=-q_vol/k/2.*x2.*x2+A.*x2+T0;
    %Calculo de los residuos
    E1=abs(T2t-T2(n/2+1,:))./T2t;
    E2=(T2t(n+1)*ones(n+1,1)-T2(:,n+1))/T2t(n+1);
    %Graficos
    figure(1)
    plot(x2,E1)
    hold on
    figure(5)
    semilogy(x2,abs(T2t-T2(1,:))./T2t)
    hold on
    figure(2)
    plot(x2,E2)
    hold on
    E(i)=E2(end);
    N(i)=n;
end
%parametros de los graficos
dir='C:\Users\Gaston\Documents\IB\IB\...
     6 - Sexto Cuatrimestre\FEM en Solidos\Informe_2\';
figure(1)
xlabel('x [m] (y=50 mm)')
ylabel('Residuo')
figure4latex([dir,'calor1'])
figure(2)
xlabel('y [m] (x=100 mm)')
ylabel('Residuo')
figure4latex([dir,'calor2'])
figure(3)
H=a./N;
loglog(N,E)
xlabel('Parametro n')
ylabel('Residuo en (x=100 mm;y=100 mm)')
figure4latex([dir,'calor3'])
%orden de convergencia del metodo
n=polyfit(log(H),log(E),1);
figure(4)
[X,Y]=meshgrid(x2,x2);
mesh(X,Y,T2)
axis square
xlabel('x [m]')
ylabel('y [m]')
zlabel('T [°C]')
figure4latex([dir,'calor4'])
figure(5)
xlabel('x [m] (y=0 mm)')
ylabel('Residuo')
figure4latex([dir,'calor5'])