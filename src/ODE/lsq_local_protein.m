function lsq_local_protein

clear all
clc
% initial search value
params = csvread('./parameters_global.csv');
k0 = params(1,11:12);
lb = [10^(-4), 10^(-4)];
ub = [10^4, 10^4];
A = xlsread('./data/il10(light time).xlsx');
t = A(1:11,1);
yexp = A(1:11,3);
%options = optimoptions('fmincon')
options = optimoptions('lsqcurvefit','Algorithm','Levenberg-Marquardt','FiniteDifferenceType','central');
[khat, resnorm] = lsqcurvefit(@ObjFunc4Fmincon, k0, t, yexp, lb, ub, options);
fprintf('\tk1 = %.4f\n',khat(1))
fprintf('\tk2 = %.4f\n',khat(2))
fprintf('  The sum of the squares is: %.1e\n\n',resnorm)
k_fmincon = khat;
csvwrite('./parameter_local(protein).csv',khat)
% simulation
y0 = [0 1 0 0];
tspan = [0 : 5 : 50];
[tt yy] = ode45(@(t,Y) KineticEqs(t,Y,khat),tspan,y0);  
plot(t,yexp(:,1),'ro')
hold on
plot(tt,yy(:,4),'bo-')
xlabel('t/h')
ylabel('pg/ml')

end


function f = ObjFunc4Fmincon(k,t0)

y0 = [0 1 0 0];
[t, Y] = ode45(@(t,Y) KineticEqs(t,Y,k),t0,y0);
f = Y(:,4);

end

function dYdt = KineticEqs(t,Y,k)

params(1,1:8) = csvread('./parameter_local(mRNA).csv');
params(1,9:10) = csvread('./parameter_local(preprotein).csv');
c = zeros(10);
c(1:8) = params(1,1:8);
c(9:10) = params(1,9:10);

R = 3.5 * 10^(-4);
V = zeros(2,1);
V(1) = 5 * 10^4;
V(2) = V(1)*c(1) - 2*Y(1);
dYdt = [ ((c(2)*V(2))/(c(3)+V(2))-c(4)*Y(1))
((c(7)+(c(6)*Y(1))/(c(5)+Y(1)))-c(8)*Y(2))
(c(9)*Y(2)-c(10)*Y(3))
(k(1)*Y(3)-k(2)*Y(4))];

end