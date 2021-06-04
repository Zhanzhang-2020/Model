function PS0global
tic
clear all
clc

lb = [1e-4 1e-4 1e-4 1e-5 1e-4 1e-4 1e-4 1e-4 1e-4 1e-4 1e-4];
ub = [1e4 1e4 1e4 1e5 1e4 1e4 1e4 1e4 1e4 1e4 1e4];
nvars = 11;
%x0 = zeros(20,nvars); % set 20 individuals as row vectors

options = optimoptions('particleswarm','SwarmSize',150,'FunctionTolerance',1e-8,'MaxStallIterations',500,'PlotFcn',@pswplotbestf);
%options.InitialSwarmMatrix = x0;
rng default 
options.HybridFcn = @fmincon; % particleswarm calls the hybrid function after particleswarm finishes its iterations
[k,fval,exitflag,output] = particleswarm(@ObjFunc4Fmincon,nvars,lb,ub,options);
fprintf('The number of iteration was : %d\n', output.iterations);
fprintf('The number of function evaluations was : %d\n', output.funccount);
fprintf('The best fitness value was : %d\n', fval);
csvwrite('./parameters_global.csv',k)
toc
t = toc;
fprintf('TIME Consumed = %.4f\n',t)
disp(output.hybridflag)

end

function f = ObjFunc4Fmincon(k)
y0 = [0 1 0 0];
A = xlsread('./data/il10(light time).xlsx');
t = A(1:11,1);
yexp = A(1:11,2);
tspan = 0 : 5 : 50;
[t, Y] = ode45(@(t,Y) KineticEqs(t,Y,k),tspan,y0,[]);
f = sum((Y(4)-yexp).^2);
end

function dYdt = KineticEqs(t,Y,k)

R = 3.5 * 10^(-4);
V = zeros(2,1);
V(1) = 5 * 10^4;
V(2) = V(1)*k(1) - 2*Y(1);
dYdt = [ ((k(2)*V(2))/(k(3)+V(2))-k(4)*Y(1))
((k(7)+(k(6)*Y(1))/(k(5)+Y(1)))-k(8)*Y(2))
(k(9)*Y(2)-k(10)*Y(3))
(V(1)*R*k(10)*Y(3)-k(11)*Y(4))];

end