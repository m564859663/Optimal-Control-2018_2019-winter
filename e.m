function e_1 = e(k)
M=5; rad=1.5;
e_1(1)=M*sin(2*k*rad)+5*M*sin(k*rad)+2*M*sin(5*k*rad)+M*cos(k*rad);
e_1(2)=3*M*sin(2*k*rad)+5*M*cos(3*k*rad)+2*M*sin(5*k*rad);
e_1=e_1';
end

