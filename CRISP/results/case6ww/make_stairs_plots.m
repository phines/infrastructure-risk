data = xlsread('test2_restoration_case6ww_1.csv')
time = data(:,1)
ls = data(:,3)
nlo = data(:,4)
figure
subplot(2,1,1)
[Rx,Ry] = stairs(time, ls)
plot(Rx,Ry)
hold on
ylabel('Load served (%)')
axis([min(time), max(time)+1, min(ls)-0.1, 1.01])
hold off
subplot(2,1,2)
[Nx,Ny] = stairs(time, nlo)
plot(Nx,Ny)
hold on
ylabel('No. lines out')
xlabel('time (min)')
axis([min(time), max(time)+1, -0.1, max(nlo)+1])
hold off
