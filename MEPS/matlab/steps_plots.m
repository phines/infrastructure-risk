% make step plot of load served and lines out over time
% data = xlsread('../results/case6ww/test2_restoration_case6ww_1.csv');
data = xlsread('../results/RestorationFixed_case39_9.csv');
time = [-100; data(:,1)];
linesout = [0; data(:,4)];
loadserved = [1; data(:,3)];

[Lin,t1] = stairs(linesout,time);
[Load,t2] = stairs(100*loadserved,time);
t1 = t1/60;
t2 = t2/60;
figure
subplot(2,1,1)
plot(t2,Load)
hold on
ylabel('Load served (%)')
set(gca, 'fontsize',13)
set(gca,'xticklabels',[])
p = get(gca,'position'); p(2) = 0.5;
set(gca, 'position', p)
ylim(100*[0 1.1])
xlim([-0.3 23])
hold off
subplot(2,1,2)
plot(t1,Lin)
hold on
ylabel('No. lines out')
xlabel('time (hours)')
set(gca, 'fontsize',13)
ylim([-1 40])
xlim([-0.3 23])
hold off
