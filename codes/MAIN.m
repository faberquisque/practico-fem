%parametros iniciales
a=30e-3;
b1=35.2e-3;
b2=35e-3;
c=40e-3;
%calculo solucion analitica de b y P_z
[P_z,b_t]=analitica(a,b1,b2,c);
tol=1e-10;
%semilla para iniciar la iteracion
p_0=8.236020397600508e+07;
%numero de elementos
n=10;
error=1;
%Iteracion hasta encontrar P y b por el metodo numerico
while abs(error)>tol
    u1=vector_U(a,b1,0,p_0,n);
    u2=vector_U(b2,c,p_0,0,n);
    error=(b1+u1(end)-(b2+u2(1)))/b_t;
    p_0=p_0*(1+error);
end
%Calculo del campo de tensiones por el metodo numerico propio
S1=vector_Sig(vector_def(u1,a,b1));
S2=vector_Sig(vector_def(u2,b2,c));
%Calculo del b y de los residuos de b y P_z del metodo numerico
b=(b1+u1(end)+b2+u2(1))/2;
error_b=(b-b_t)/b_t;
error_P=(p_0-P_z)/P_z;
% Calculo del campo de desplazamiento y tensiones 
% por el metodo analitico
m1=linspace(a,b1,n+1);
m2=linspace(b2,c,n+1);
u1_t=vector_U_analit(m1,a,b1,0,P_z)';
u2_t=vector_U_analit(m2,b2,c,P_z,0)';
S1_t=vector_Sig(vector_def(u1_t,a,b1));
S2_t=vector_Sig(vector_def(u2_t,b2,c));
%Datos extraidos del Abaqus de desplazamiento y tensiones
u1_ab=[-8.99807e-005,-8.95299e-005,-8.91129e-005,-8.87272e-005,...
    -8.83722e-005,-8.80454e-005,-8.77466e-005,-8.74734e-005,...
    -8.72257e-005,-8.70012e-005,-8.67993e-005];
u2_ab=[1.13150e-004,1.12493e-004,1.11868e-004,1.11273e-004,...
    1.10707e-004,1.10168e-004,1.09656e-004,1.09169e-004,...
    1.08707e-004,1.08269e-004,1.07853e-004];
S1_ab=[-5.25624e+06,-1.51319e+07,-2.45207e+07,-3.3455e+07,...
    -4.19626e+07,-5.00723e+07,-5.78048e+07,-6.51838e+07,...
    -7.22244e+07,-7.89355e+07;
    -5.94612e+08,-5.84741e+08,-5.75354e+08,-5.66423e+08,...
    -5.57919e+08,-5.49814e+08,-5.42093e+08,-5.34724e+08,...
    -5.27707e+08,-5.20983e+08];

S2_ab=[-7.74315e+07,-6.78079e+07,-5.85631e+07,-4.96883e+07,...
    -4.11667e+07,-3.29807e+07,-2.5115e+07,-1.75525e+07,...
    -1.02786e+07,-3.27876e+06;
    6.16652e+08,6.07016e+08,5.97792e+08,5.88928e+08,...
    5.80417e+08,5.72237e+08,5.64376e+08,5.56818e+08,...
    5.49548e+08,5.42554e+08];
%Calculo de los residuos a graficar
n1=linspace(a,b1,n);
n2=linspace(b2,c,n);
E1=(S1-S1_t)./S1_t;
E1_ab=(S1_ab-S1_t)./S1_t;
E2=(S2-S2_t)./S2_t;
E2_ab=(S2_ab-S2_t)./S2_t;
leyenda=['Tubo In-propio';
         'Tubo In-abaqus';
         'Tubo Ex-propio';
         'Tubo Ex-abaqus'];
dir='C:\Users\Gaston\Documents\IB\IB\...
     6 - Sexto Cuatrimestre\FEM en Solidos\Informe_2\';
%Figuras
figure(1)
semilogy(n1,E1(1,:),n1,E1_ab(1,:),n2,abs(E2(1,:)),n2,abs(E2_ab(1,:)))
title('Residuo de la Tension Radial')
xlabel('r [m]')
ylabel('Residuo S_r')
legend(leyenda)
legend('Location','North')
xlim([a,c])
figure4latex([dir,'figure1'])
figure(2)
semilogy(n1,E1(2,:),n1,E1_ab(2,:),n2,abs(E2(2,:)),n2,abs(E2_ab(2,:)))
title('Residuo de la Tension Tangencial')
xlabel('r [m]')
ylabel('Residuo S_t')
legend(leyenda)
legend('Location','East')
xlim([a,c])
figure4latex([dir,'figure2'])
figure(3)
e1=abs((u1-u1_t)./u1_t);
e2=abs((u2-u2_t)./u2_t);
e1_ab=abs((u1_ab'-u1_t)./u1_t);
e2_ab=abs((u2_ab'-u2_t)./u2_t);
semilogy(m1,e1,m1,e1_ab,m2,e2,m2,e2_ab)
title('Residuo del Desplazamiento')
xlabel('r [m]')
ylabel('Residuo u')
legend(leyenda)
legend('Location','East')
xlim([a,c])
figure4latex([dir,'figure3'])
figure(4)
e=zeros(10,1);
for i=1:50
    r=linspace(a,b1,i+1);
    u=vector_U(a,b1,0,P_z,i);
    ut=vector_U_analit(r',a,b1,0,P_z);
    e(i)=abs((u(end)-ut(end))/ut(end));
end
loglog(1:50,e)
title('Residuo del desplazamiento en funcion del numero de elementos')
xlabel('numero de elementos')
ylabel('Residuo u1 (en r=b1)')
figure4latex([dir,'figure4'])
%Calculo del orden de convergencia
N=polyfit(log((b1-a)./(1:50))',log(e),1);
%Calculo del residuo de b obtenido en abaqus
error_b_ab=((b1+u1_ab(end)+b2+u2_ab(1))/2-b_t)/b_t;
