
n = 1000;
x = exp(randn(n,1));

[x,pdf] = emperical_pdf(x,40);

figure(1); clf;
plot(x,pdf)
set(gca,'xscale','log');
set(gca,'yscale','log');

area = trapz(x,pdf)
% 
% figure(2); clf;
% [yvalues,xvalues] = histc(x,20)
% p = histcounts(yvalues,xvalues,'Normalization','pdf')
% plot(xvalues,p);
