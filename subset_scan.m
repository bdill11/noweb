function idx = subset_scan(list1,list2)
% function idx = subset_scan(list1,list2)
%
% Returns index list idx where
% if idx(i) != 0 then list1(i) == list2(idx(i)), and
% if idx(i) == 0 then list1(i) ~= list2(j) for any j.
idx = zeros(size(list1));
[list1,px1] = sort(list1);
[list2,px2] = sort(list2);
i1 = 1;
i2 = 1;
while i1 <= length(list1) && i2 <= length(list2)
    if list1(i1) < list2(i2)
        i1 = i1 + 1;
    elseif list1(i1) > list2(i2)
        i2 = i2 + 1;
    else
        idx(px1(i1)) = px2(i2);
        i1 = i1 + 1;
        i2 = i2 + 1;
    end
end
