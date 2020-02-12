
p1 = 1;
p2 = 1;
factor = 1.2; % calibrate
n = 1e6;
results = zeros(n,1);

for i = 1:n
    rest_time = wblrnd(p1,p2);
    % repeatedly increase rest_time, if rest_time is big
    while 1
        p_inc = 1-1/log(rest_time); % magic; calibrate
        if rand<p_inc
            rest_time = rest_time*factor;
        else
            break;
        end
    end
    % record the results
    results(i) = rest_time;
end

%% figure
figure(1); clf;
hist(results,20);

figure(2); clf;
x = sort(results);
p = (n:-1:1)/n;
loglog(x,p);
