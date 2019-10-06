function square_plots(data,arrangement,colors)

if nargin<3
    colors = [0 1 1; 1 1 0; 0 1 0; 0 0 1];
end
n_sets = size(data,1);


%if any(any(data<0))
    %error('Does not work for <0 data');
%end
if prod(arrangement)~=n_sets
    error('Bad dimensions');
end

% normalize the data
mx = max(max(data));
data = data./mx;
% choose the centers
all_centers = zeros(n_sets,2);
n_rows = arrangement(1);
n_cols = arrangement(2);
for row = 1:n_rows
    index = (1:n_cols) + (row-1)*n_cols;
    all_centers(index,1) = 0:2:((n_cols-1)*2);
    all_centers(index,2) = (row-1)*2;
end

for row = 1:n_sets
    center_xy = all_centers(row,:);
    % find the patch dimensions
    % lower left
    if data(row,1)<=0
    else
    area = data(row,1);
    side = sqrt(area);
    lower_left_xy = center_xy - side;
    pos = [lower_left_xy side side];
    rectangle('Position',pos,'FaceColor',colors(1,:),'EdgeColor','none');
    end
    % lower right
    if data(row,2)<=0
    else
    area = data(row,2);
    side = sqrt(area);
    lower_left_xy = center_xy - [0 side];
    pos = [lower_left_xy side side];
    rectangle('Position',pos,'FaceColor',colors(2,:),'EdgeColor','none');
    end
    % upper left
    if data(row,3)<=0
    else
    area = data(row,3);
    side = sqrt(area);
    lower_left_xy = center_xy - [side 0];
    pos = [lower_left_xy side side];
    rectangle('Position',pos,'FaceColor',colors(3,:),'EdgeColor','none');
    end
    % upper right
    if data(row,4)<=0
    else
    area = data(row,4);
    side = sqrt(area);
    lower_left_xy = center_xy;
    pos = [lower_left_xy side side];
    rectangle('Position',pos,'FaceColor',colors(4,:),'EdgeColor','none');
    end
end

axis off;

end