% make step plot of load served and lines out over time
data = xlsread('../results/RestorationFixed_case39_9.csv');
time = data(:,1);
linesout = data(:,4);
loadserved = data(:,3);

[Lin,t1] = stairs(linesout,time);
[Load,t2] = stairs(loadserved,time);
figure
subplot(2,1,1)
plot(t2,Load)
hold on
ylabel('Load served (%)')
hold off
set(gca, 'fontsize',13)
subplot(2,1,2)
plot(t1,Lin)
hold on
ylabel('No. lines out')
xlabel('time (min)')
hold off
set(gca, 'fontsize',13)