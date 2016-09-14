function u=vector_U_analit(r,a,b,p_in,p_out)
nu=0.3;
E=200e9;
C1=(1-nu)/E*(a*a*p_in-b*b*p_out)/(b*b-a*a);
C2=(1+nu)/E*a*a*b*b*(p_in-p_out)/(b*b-a*a);
u=C1.*r+C2./r;