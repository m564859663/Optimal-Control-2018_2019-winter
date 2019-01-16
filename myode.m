function dx = myode(x,A,B,Li,Q,R)
dx = zeros(5,1);
dx(1:4) = (A-B*Li)*x(1:4);
dx(5) = x(1:4)'*(Q+Li'*R*Li)*x(1:4);
end