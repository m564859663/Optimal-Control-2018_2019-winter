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
Q = eye(4);
R = eye(2);
L=[5 -2 -1 -2;
    5 -5 2 -1];%(A-BK1)Hurwitz
dT = 0.01; %time period T
X0 = [0;10;0;10;0]; %system inital state
N=20; % linear equations numbers
time_end = 10; % end time
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
L_1=zeros(size(L));
P=zeros(n);
P(1)=1000;
iterate_num=0;
syms k11 k12 k13 k14 k22 k23 k24 k33 k34 k44; %10 variations
while(1)
    P_a=[k11 k12 k13 k14;
         k12 k22 k23 k24;
         k13 k23 k33 k34;
         k14 k24 k34 k44];
    for i = 1:N
        tspan = [(i - 1) * dT + num1 * dT * N, i * dT + num1 * dT * N];
        [t,x] = ode45(@(t,x) myode(x,A,B,L,Q,R),tspan,X);
        X = x(end,:)';
        figure(1);
        plot(t,x(:,1),'r',t,x(:,2),'b',t,x(:,3),'g',t,x(:,4),'k'); hold on;

        X_i_t = x(1, 1:4)*P_a*x(1, 1:4)';
        X_i_t_T = x(end, 1:4)*P_a*x(end, 1:4)';
        Y = x(end,end) - x(1,end);
        equ= X_i_t - X_i_t_T -Y;
        [Ai, Bi]=equationsToMatrix(equ,[k11 k12 k13 k14 k22 k23 k24 k33 k34 k44]);
        Ai_n(i,:)=double(Ai);
        Bi_n(i,:)=double(Bi);
    end
    P_1=P;
    P_solve = pinv(Ai_n)*Bi_n;
%     P_solve = Ai_n\Bi_n;

    P=[P_solve(1) P_solve(2) P_solve(3) P_solve(4);
      P_solve(2) P_solve(5) P_solve(6) P_solve(7);
      P_solve(3) P_solve(6) P_solve(8) P_solve(9);
      P_solve(4) P_solve(7) P_solve(9) P_solve(10)]    
  
    if t(end)==10
        P=P_1;
        break;
    end
    
    if (all(eig(P)>0) == 0) || P(1)>P_1(1) % Matrix_P is not positive definite
        iterate_num=iterate_num+1;
        if iterate_num == 4
            P=P_1
            break;
        end        
        Ai_n=[zeros(N, 10); Ai_n];
        Bi_n=[zeros(N,1); Bi_n];
        P=P_1
        num1=num1+1;
        continue;
    end
    iterate_num=0;
    clear Ai_n Bi_n;

    L_1=L; 
    L = R \ (B') * P;
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
norm_num=norm((L - L_1),1)
P
L
J=X0(1:n)' * P * X0(1:n)
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