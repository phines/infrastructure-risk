function sorted_data = sort_data(size,nlines,data)

pick_data = false(length(nlines),1);
for k = 1:length(size)
    pick_data = max((nlines == size(k)),pick_data);
end
n_data = data(pick_data);
sorted_data = sort(n_data);
end