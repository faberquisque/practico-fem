function D=matriz_D
nu=0.3;
E=200e9;
%deformacion plana
%D=E/(1+nu)/(1-2*nu)*[1-nu,nu,0;nu,1-nu,0;0,0,(1-2*nu)/2];
%tension plana
D=E/(1-nu^2)*[1,nu;nu,1];
end
