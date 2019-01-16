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
L0=[5 -2 -1 -2;
    5 -5 2 -1];%(A-BK1)Hurwitz
L=L0;
dT = 0.01; %time period T
X0 = [0;50;0;50]; %system inital state
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
X=X0; L_1=zeros(size(L));
time_num=60; time_end = 10;
iterate_num=0;
P=zeros(n); P(1)=1000;

syms k11 k12 k13 k14 k22 k23 k24 k33 k34 k44; %10 variations
syms l11 l12 l13 l14 l21 l22 l23 l24; %8 variations
while(1)
    P_a=[k11 k12 k13 k14;
         k12 k22 k23 k24;
         k13 k23 k33 k34;
         k14 k24 k34 k44];
     L_i1=[l11 l12 l13 l14;
          l21 l22 l23 l24];
    for i = 1:N        
        tspan = [(i - 1) * dT + num1 * dT * N, i * dT + num1 * dT * N];
        t=linspace(tspan(1),tspan(2),time_num);t=t';
        del_t=(tspan(2)-tspan(1))/time_num;
        x=zeros(time_num, 4);
        x(1,:)=X';
        sum1=0; sum2=0;
        for j=1:time_num-1
            x_k=x(j,:)';
            x_k_1=x_k + del_t*(A*x_k - B*L*x_k + B*e(t(j)));
            x(j+1,:)=x_k_1';
            sum1=sum1 + x_k' * (Q+L'*R*L) * x_k * del_t;
            sum2=sum2+e(t(j))' * R * L_i1 *x_k * del_t;
        end
        
        X = x(end,:)';
        figure(1);
        plot(t,x(:,1),'r',t,x(:,2),'b',t,x(:,3),'g',t,x(:,4),'k'); hold on;

        X_i_t = x(1, :)*P_a*x(1, :)';
        X_i_t_T = x(end, :)*P_a*x(end, :)';
        
        equ= X_i_t_T - X_i_t + sum1 - 2 * sum2;
        [Ai, Bi]=equationsToMatrix(equ,[k11 k12 k13 k14 k22 k23 k24 k33 k34 k44 l11 l12 l13 l14 l21 l22 l23 l24]);
        Ai_n(i,:)=double(Ai);
        Bi_n(i,:)=double(Bi);
    end
    P_1=P;L_1=L;
    P_solve = Ai_n\Bi_n;
%     P_solve = pinv(Ai_n)*(Bi_n);
    P=[P_solve(1) P_solve(2) P_solve(3) P_solve(4);
      P_solve(2) P_solve(5) P_solve(6) P_solve(7)
      P_solve(3) P_solve(6) P_solve(8) P_solve(9)
      P_solve(4) P_solve(7) P_solve(9) P_solve(10)]
    L_i1=[P_solve(11) P_solve(12) P_solve(13) P_solve(14);
            P_solve(15) P_solve(16) P_solve(17) P_solve(18)];

    if t(end)==10
        break;
    end
    
    if (all(eig(P)>0) == 0) || P(1)>P_1(1) % Matrix_P is not positive definite
        iterate_num=iterate_num+1;
        if iterate_num == 4
            P=P_1
            break;
        end        
        Ai_n=[zeros(N, 18); Ai_n];
        Bi_n=[zeros(N, 1); Bi_n];
        P=P_1
        num1=num1+1;
        continue;
    end
    iterate_num=0;
    P_history(:,:,num1+1)=P;   
    clear Ai_n Bi_n;

    L = L_i1;
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
    x(:,5)=0;X(5)=0;
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