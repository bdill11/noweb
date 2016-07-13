function [edges,t_idx1,t_idx2] = internal_edges(p,t)
% function [edges,t_idx1,t_idx2] = internal_edges(p,t)
%
% Return list of internal edges as edges, m x 2 array.
% edges(i,:) is a list of indexes into rows of p
% edges(i,:) is common to triangle t(t_idx1(i),:) and t(t_idx2(i),:)
nt = size(t,1);
temp = [t(:,[1 2]), (1:nt)'; t(:,[2 3]), (1:nt)'; t(:,[3 1]), (1:nt)'];
temp(:,[1,2]) = sort(temp(:,[1,2]),2);  % sort each row
temp = sortrows(temp);
[temp2,idx1] = unique(temp(:,[1,2]),'rows','first');
[temp2,idx2] = unique(temp(:,[1,2]),'rows','last');
idx = find(idx1 ~= idx2);
t_idx1 = temp(idx1(idx),3);
t_idx2 = temp(idx2(idx),3);
edges = temp(idx1(idx),[1,2]);
end

