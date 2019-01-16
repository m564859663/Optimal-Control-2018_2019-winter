%%
clear all,close all,clc
A=[-0.0389 0.0271 0.0188 -0.4555;
    0.0482 -1.010 0.0019 -4.0208;
    0.1024 0.3681 -0.707 1.4200;
    0.0000 0.0000 1.0000 0.0000];
B=[ 0.4422 0.1291;
    3.5446 -7.5922;
    -6.0214 4.4900;
    0       0];
n = 4; %system state num
Q=eye(4); R=eye(2);
L=[5 -2 -1 -2;
    5 -5 2 -1];%(A-BK1)Hurwitz
P=zeros(n);
dT = 0.01; %time period T
X0 = [0;0.1;0;0.1;0]; %system inital state
time_end = 10;
threshold = 1e-6; %epsilon
num1 = 0; %count iteration time

if all(real(eig(A-B*L))<0)
    disp('The initial L is Horwitz!')
else
    disp('The initial L is not Horwitz!')
    return
end

%%
X=X0;
while(1)    
    tspan = [num1 * dT, dT + num1 * dT];
    [t,x] = ode45(@(t,x) myode(x,A,B,L,Q,R),tspan,X);
    figure(1);
    plot(t,x(:,1),'r',t,x(:,2),'b',t,x(:,3),'g',t,x(:,4),'k'); hold on;
    X = x(end,:)';
    
    P_1=P; L_1=L;
    P=lyap((A-B*L)', (Q+L'*R*L));
    L = R \ (B') * P;
    
    if (all(eig(P)>0) == 0) % Matrix_P is not positive definite
        P=P_1;
        break;
    end
   
    if(norm((L - L_1),1) < threshold )
        break;
    end
    
    figure(2)
    plot(num1,P(1,1),'ro',num1,P(2,3),'bx',num1,P(2,4),'gs',num1,P(4,4),'k*'),hold on;
    num1 = num1 + 1;
    
end

%%
string1 = sprintf('Total iterate %d times!', num1);disp(string1);
string2 = sprintf('time is %f', tspan(end));disp(string2);
norm=norm((L - L_1),1)
P
L
J=X0(1:n)' * P * X0(1:n)
K = care(A,B,Q,R) % Riccati
L_r=R\B'*K
if tspan(end)<time_end
    tspan = [tspan(end), time_end];
    [t,x] = ode45(@(t,x) myode(x,A,B,L,Q,R),tspan,X);
    figure(1);
    plot(t,x(:,1),'r',t,x(:,2),'b',t,x(:,3),'g',t,x(:,4),'k'); hold on;
    legend('x1', 'x2', 'x3', 'x4')
    xlabel('time(s)')


    figure(2);
    legend('P(1,1)', 'P(2,3)', 'P(2,4)', 'P(4,4)')
    xlabel('time(number)')
    set(gca,'XTick',0:1:num1)
end