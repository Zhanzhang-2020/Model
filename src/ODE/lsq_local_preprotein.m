function lsq_local_preprotein

clear all
clc
% initial search value
params = csvread('./parameters_global.csv');
k0 = params(1,9:10);
lb = [10^(-3) 10^(-3)];
ub = [10^(3) 10^(3)];
A = xlsread('/Users/zhangzhan/Desktop/matlab/model/data/il10(flow).xlsx');
t = A(1:11,1);
yexp = A(1:11,2);
%options = optimoptions('fmincon')
options = optimoptions('lsqcurvefit','Algorithm','Levenberg-Marquardt','FiniteDifferenceType','central');
[khat, resnorm] = lsqcurvefit(@ObjFunc4Fmincon, k0, t, yexp, lb, ub, options);
fprintf('\tk1 = %.4f\n',khat(1))
fprintf('\tk2 = %.4f\n',khat(2))
fprintf('  The sum of the squares is: %.1e\n\n',resnorm)
k_fmincon = khat;
csvwrite('./parameter_local(preprotein).csv',khat)
% simulation
y0 = [0 1 0];
tspan = [0 : 5 : 50];
[tt yy] = ode45(@(t,Y) KineticEqs(t,Y,khat),tspan,y0);  
plot(t,yexp(:,1),'ro')
hold on
plot(tt,yy(:,3),'bo-')
xlabel('t/h')
ylabel('A.U.')

end

function f = ObjFunc4Fmincon(k,t0)

y0 = [0 1 0];
[t, Y] = ode45(@(t,Y) KineticEqs(t,Y,k),t0,y0);
f = Y(:,3);

end

function dYdt = KineticEqs(t,Y,k)

params_1 = csvread('./parameter_local(mRNA).csv');
c = params_1(1,:);
R = 3.5 * 10^(-4);
V = zeros(2,1);
V(1) = 1;
V(2) = V(1)*c(1) - 2*Y(1);
dYdt = [ ((c(2)*V(2))/(c(3)+V(2))-c(4)*Y(1))
((c(7)+(c(6)*Y(1))/(c(5)+Y(1)))-c(8)*Y(2))
(k(1)*Y(2)-k(2)*Y(3))];

end