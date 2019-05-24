function [Costs, bin_ls, bin_rt] = make_contour( load_shed, restoration_time, ENS )

LS = load_shed(load_shed~=0);
minLS = min(LS);
maxLS = max(LS);
RT = restoration_time(load_shed~=0);
minRT = min(RT);
maxRT = max(RT);
costs = ENS(load_shed~=0);
plot3(LS,RT,costs,'*')
bin_ls = zeros(length(LS));
LS1 = LS;
bin_rt = zeros(length(RT));
RT1 = RT;
l=0;r=0;
for i = 1:length(LS)
    if ~isempty(LS1)
        bin_ls(i) = min(LS1);
        LS1 = LS1(LS1~=min(LS1));
    elseif l==0
        l=1;
        bin_ls = bin_ls(1:i-1);
    end
    if ~isempty(RT1)
        bin_rt(i) = min(RT1);
        RT1 = RT1(RT1~=min(RT1));
    elseif r==0
        r=1;
        bin_rt = bin_rt(1:i-1);
    end
end
Costs = zeros(length(bin_ls),length(bin_rt));
for n = 1:length(bin_ls)
    l = (abs(LS - bin_ls(n))<=0.001);
    for m = 1:length(bin_rt)
        r = (abs(RT - bin_rt(m))<=0.0001);
        if isempty(costs(l & r))
        else
        Costs(n,m) = sum(costs(l & r))./length(costs);
        end
    end
end
figure
contour(bin_ls,bin_rt,Costs')

