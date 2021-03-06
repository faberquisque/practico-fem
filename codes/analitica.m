function [P,b]=analitica(a,b1,b2,c)
% a=0.03;
% b1=0.0352;
% b2=0.035;
% c=0.04;
E=200e9;
nu=0.3;
C1=(1-nu)*b2^3/(c^2-b2^2)+(1+nu)*c^2*b2/(c^2-b2^2);
C2=(1-nu)*b1^3/(b1^2-a^2)+(1+nu)*a^2*b1/(b1^2-a^2);
b=(C1*b1+C2*b2)/(C1+C2);
P=E*(b-b2)/C1;