function lsq_local_mRNA

clear all
clc
% initial search value
params = csvread('./parameters_global.csv');
k0 = params(1,1:8);
lb = [10^(1) 10^2 10^(-3) 10^(-6) 10^(-4) 10^(-4) 10^(-3) 10^(-3)];
ub = [10^4 10^5 10^(4) 10^(-1) 10^(-1) 10^4 10^1 10^3];
A = xlsread('./data/il10(mRNA).xlsx');
t = A(1:7,1);
yexp = A(1:7,2);
%options = optimoptions('fmincon')
options = optimoptions('lsqcurvefit','Algorithm','Levenberg-Marquardt','FiniteDifferenceType','central');
[khat, resnorm] = lsqcurvefit(@ObjFunc4Fmincon, k0, t, yexp, lb, ub, options);
fprintf('\tk1 = %.4f\n',khat(1))
fprintf('\tk2 = %.4f\n',khat(2))
fprintf('\tk3 = %.4f\n',khat(3))
fprintf('\tk4 = %.4f\n',khat(4))
fprintf('\tk5 = %.4f\n',khat(5))
fprintf('\tk6 = %.4f\n',khat(6))
fprintf('\tk7 = %.4f\n',khat(7))
fprintf('\tk8 = %.4f\n',khat(8))
fprintf('  The sum of the squares is: %.1e\n\n',resnorm)
k_fmincon = khat;
csvwrite('./parameter_local(mRNA).csv',khat)
% simulation
y0 = [0 1];
tspan = [0 : 5 : 30];
[tt yy] = ode45(@(t,Y) KineticEqs(t,Y,khat),tspan,y0);  
plot(t,yexp(:,1),'ro')
hold on
plot(tt,yy(:,2),'bo-')
xlabel('t/h')
ylabel('fold-change')

end

function f = ObjFunc4Fmincon(k,t0)

y0 = [0 1];
[t, Y] = ode45(@(t,Y) KineticEqs(t,Y,k),t0,y0);
f = Y(:,2);

end

function dYdt = KineticEqs(t,Y,k)

R = 3.5 * 10^(-4);
V = zeros(2,1);
V(1) = 1;
V(2) = V(1)*k(1) - 2*Y(1);
dYdt = [ ((k(2)*V(2))/(k(3)+V(2))-k(4)*Y(1))
((k(7)+(k(6)*Y(1))/(k(5)+Y(1)))-k(8)*Y(2))];

end