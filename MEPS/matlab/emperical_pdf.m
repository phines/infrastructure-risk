function [x,density] = emperical_pdf(data,points_per_point)

sorted_data = sort(data);
n = length(data);
x = [];
density = [];

start = 1;
while 1
    finish = start + points_per_point;
    if finish>n
        break;
    end
    data_subset = sorted_data(start:finish);
    dx = data_subset(end) - data_subset(1); % probably need to fix dx for bias
    mid_point = mean([data_subset(end),data_subset(1)]);
    pdf_at_midpoint = points_per_point/n/dx;
    x = [x;mid_point]; %#ok<AGROW>
    density = [density;pdf_at_midpoint]; %#ok<AGROW>
    start = start + 1;
end

end