function data = sort_data2(size1,data1,size2,data2,costs)

pick_data1 = (data1 >= size1(1) & data1 <= size1(2));
pick_data2 = (data2 >= size2(1) & data2 <= size2(2));
n_data = costs(pick_data1 & pick_data2);
data = sort(n_data);
end