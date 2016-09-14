function I=gauss_int(fun,orden_int)
switch orden_int
    case 2
        w=[1,1];
        xi=0.5773502692*[-1,1];
    case 3
        w=[5/9,8/9,5/9];
        xi=0.7745966692*[-1,0,1];
    case 4
        w=[0.3478548451,0.6521451549,0.6521451549,0.3478548451];
        xi=[-0.8611363116,-0.3399810436,0.3399810436,0.8611363116];     
end
I=0;
for i=1:orden_int
    I=I+w(i)*fun(xi(i));
end
end
