function normal_plot

clear all
clc 
% initial search value
k(1:8) = csvread('./parameter_local(mRNA).csv');
k(9:10) = csvread('./parameter_local(preprotein).csv');
k(11:12) = csvread('./parameter_local(protein).csv');
y0 = [0 1 0 0];

tspan = [0 : 5 : 50];
[tt yy] = ode45(@KineticEqs,tspan,y0,[],k); 

%figure(1);
figure('Renderer', 'painters', 'Position', [10 10 9000 6000])
subplot(2,2,1,'align');
plot(tt,yy(:,1),'ro-');
xlabel('t/h')
ylabel('The concentration of GAVPO(bound)')
legend({'simulation'},'Location','southeast')
title('The dynamic change of GAVPO(bound)')
set(gca,'FontSize',15,'FontName','Arial');

subplot(2,2,2,'align');
A_1 = xlsread('/Users/zhangzhan/Desktop/matlab/model/data/il10(mRNA).xlsx');
t_1 = A_1(1:7,1);
yexp_1 = A_1(1:7,2);
plot(tt,yy(:,2),'ro-');
hold on
plot(t_1,yexp_1(:,1),'b')
xlabel('t/h')
ylabel('Normalized mRNA (fold-change)')
legend({'simulation','Experimental Data'},'Location','southeast')
title('The dynamic change of mRNA')
set(gca,'FontSize',15,'FontName','Arial');

subplot(2,2,3,'align');
A_2 = xlsread('/Users/zhangzhan/Desktop/matlab/model/data/il10(flow).xlsx');
t_2 = A_2(1:11,1);
yexp_2 = A_2(1:11,2);
plot(tt,yy(:,3),'ro-');
hold on
plot(t_2,yexp_2(:,1),'b')
xlabel('t/h')
ylabel('The concentration of Intracellular protein (a.f.u)')
legend({'simulation','Experimental Data'},'Location','southeast')
title('The dynamic change of Intracellular protein')
set(gca,'FontSize',15,'FontName','Arial');

subplot(2,2,4,'align');
A_3 = xlsread('/Users/zhangzhan/Desktop/matlab/model/data/il10(light time).xlsx');
t_3 = A_3(1:11,1);
yexp_3 = A_3(1:11, 3);
plot(tt,yy(:,4),'ro-');
hold on
plot(t_3,yexp_3(:,1),'b')
xlabel('t/h')
ylabel('The concentration of Secreted protein (pg/ml)')
legend({'simulation','Experimental Data'},'Location','southeast')
title('The dynamic change of Secreted protein')
set(gca,'FontSize',15,'FontName','Arial');

saveas(gcf,'./figure/test.png');

end

function dYdt = KineticEqs(t,Y,k)

R = 3.5 * 10^(-4);
V = zeros(2,1);
V(1) = 5 * 10^4;
V(2) = V(1)*k(1) - 2*Y(1);
dYdt = [ ((k(2)*V(2))/(k(3)+V(2))-k(4)*Y(1))
((k(7)+(k(6)*Y(1))/(k(5)+Y(1)))-k(8)*Y(2))
(k(9)*Y(2)-k(10)*Y(3))
(k(11)*Y(3)-k(12)*Y(4))];

end