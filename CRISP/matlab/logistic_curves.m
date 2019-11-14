t = 0:50;
L = 1;

f1 = L./(1+exp(-0.5.*(t-4)))
f2 = L./(1+exp(-0.4.*(t-16)))
f3 = L./(1+exp(-0.3.*(t-24)))
%f4 = L./(1+exp(-0.2.*(t-48)))
%f5 = L./(1+exp(-0.1.*(t-48)))

figure
hold on
plot(t,f1)
plot(t,f2)
plot(t,f3)
plot([4 4], [0 1], 'k')
%plot(t,f4)
%plot(t,f5)
xlabel('time (hr)')
ylabel('Prob(com fails)')
legend('param p = 0.5  & x_o = 4', 'p = 0.4 & x_o = 16', 'p = 0.3 & x_o = 24')%, 'p = 2 & x_o = 48', 'p = 1 & x_o = 48')
set(gca, 'fontsize',13)
hold off