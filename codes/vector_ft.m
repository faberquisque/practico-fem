function f=vector_ft(a,b,p_in,p_out,n)
f=2*pi*[p_in*a;zeros(n-1,1);-p_out*b];
end