function square_3x3_plots(d,arrangement,colors)

if nargin<3
    colornames_wbc('xkcd','Dark Blue')
    colors = [colornames_wbc('xkcd','Dark Blue'); colornames_wbc('xkcd','Cerulean Blue');colornames_wbc('xkcd','Duck Egg Blue');...
            colornames_wbc('xkcd','Aubergine');colornames_wbc('xkcd','Darkish Purple');colornames_wbc('xkcd','Dusty Lavender');...
            colornames_wbc('xkcd','Claret');colornames_wbc('xkcd','Cherry');colornames_wbc('xkcd','Blush Pink')]
    %colors = [0 1 1; 1 1 0; 0 1 0; 0 0 1];
end
n_sets = size(d,1);


%if any(any(data<0))
    %error('Does not work for <0 data');
%end
if prod(arrangement)~=n_sets
    error('Bad dimensions');
end

% normalize the data
mx = max(max(d));
data = d./mx;
% choose the centers
all_centers = zeros(n_sets,2);
n_rows = arrangement(1);
n_cols = arrangement(2);
for row = 1:n_rows
    index = (1:n_cols) + (row-1)*n_cols;
    all_centers(index,1) = 0:5:((n_cols-1)*5);
    all_centers(index,2) = (row-1)*5;
end
for row = 1:n_sets
    center_xy = all_centers(row,:);
    % find the patch dimensions
    % upper left
    xy1 = center_xy - [1.5 -2.1];
    %pos = [lower_left_xy 0.1 0.1];
    %rectangle('Position',pos,'FaceColor','none','EdgeColor',colors(7,:));
    text(xy1(1),xy1(2),num2str(sprintf('%0.3e',sum(d(row,[3;2;1])))))
    if data(row,3)<=0
    else
    area = data(row,3);
    side = sqrt(area);
    lower_left_xy = center_xy - [1.5 -1];
    pos = [lower_left_xy side side];
    rectangle('Position',pos,'FaceColor',colors(1,:),'EdgeColor','none');
    end
    % middle left
    if data(row,2)<=0
    else
    area = data(row,2);
    side = sqrt(area);
    lower_left_xy = center_xy - [1.5 0];
    pos = [lower_left_xy side side];
    rectangle('Position',pos,'FaceColor',colors(2,:),'EdgeColor','none');
    end
    % lower left
    if data(row,1)<=0
    else
    area = data(row,1);
    side = sqrt(area);
    lower_left_xy = center_xy - [1.5 1];
    pos = [lower_left_xy side side];
    rectangle('Position',pos,'FaceColor',colors(3,:),'EdgeColor','none');
    end
    % upper middle
    xy1 = center_xy - [0.5 -2.1];
    %pos = [lower_left_xy 0.1 0.1];
    %rectangle('Position',pos,'FaceColor','none','EdgeColor',colors(7,:));
    text(xy1(1),xy1(2),sprintf('%0.3e',sum(d(row,[6;5;4]))))
    if data(row,6)<=0
    else
    area = data(row,6);
    side = sqrt(area);
    lower_left_xy = center_xy - [0.5 -1];
    pos = [lower_left_xy side side];
    rectangle('Position',pos,'FaceColor',colors(4,:),'EdgeColor','none');
    end
        % middle middle
    if data(row,5)<=0
    else
    area = data(row,5);
    side = sqrt(area);
    lower_left_xy = center_xy - [0.5 0];
    pos = [lower_left_xy side side];
    rectangle('Position',pos,'FaceColor',colors(5,:),'EdgeColor','none');
    end
        % lower middle
    if data(row,4)<=0
    else
    area = data(row,4);
    side = sqrt(area);
    lower_left_xy = center_xy - [0.5 1];
    pos = [lower_left_xy side side];
    rectangle('Position',pos,'FaceColor',colors(6,:),'EdgeColor','none');
    end
        % upper right
    xy1 = center_xy - [-0.5 -2.1];
    %pos = [lower_left_xy 0.1 0.1];
    %rectangle('Position',pos,'FaceColor','none','EdgeColor',colors(7,:));
    text(xy1(1),xy1(2),sprintf('%0.3e',sum(d(row,[9;8;7]))))
    xy = center_xy - [-1.6 -1.5];
    %pos = [lower_left_xy 0.1 0.1];
    %rectangle('Position',pos,'FaceColor','none','EdgeColor',colors(7,:));
    text(xy(1),xy(2),sprintf('%0.3e',sum(d(row,[9;6;3]))))
    if data(row,9)<=0
    else
    area = data(row,9);
    side = sqrt(area);
    lower_left_xy = center_xy - [-0.5 -1];
    pos = [lower_left_xy side side];
    rectangle('Position',pos,'FaceColor',colors(7,:),'EdgeColor','none');
    end
        % middle right
     xy = center_xy - [-1.6 -.5];
    %pos = [lower_left_xy 0.1 0.1];
    %rectangle('Position',pos,'FaceColor','none','EdgeColor',colors(7,:));
    text(xy(1),xy(2),sprintf('%0.3e',sum(d(row,[8;5;2]))))
    if data(row,8)<=0
    else
    area = data(row,8);
    side = sqrt(area);
    lower_left_xy = center_xy - [-0.5 0];
    pos = [lower_left_xy side side];
    rectangle('Position',pos,'FaceColor',colors(8,:),'EdgeColor','none');
    end
        %lower right
    xy = center_xy - [-1.6 .5];
    %pos = [lower_left_xy 0.1 0.1];
    %rectangle('Position',pos,'FaceColor','none','EdgeColor',colors(7,:));
    text(xy(1),xy(2),sprintf('%0.3e',sum(d(row,[8;5;2]))))
    if data(row,7)<=0
    else
    area = data(row,7);
    side = sqrt(area);
    lower_left_xy = center_xy - [-0.5 1];
    pos = [lower_left_xy side side];
    rectangle('Position',pos,'FaceColor',colors(9,:),'EdgeColor','none');
    end
end

axis off;

end