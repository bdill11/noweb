function [bedges,bnodes,t_idx1,t_idx2] = get_internal_boundary2d(t,t_list)
% function [bedges,bnodes,t_idx1,t_idx2] = get_internal_boundary2d(t,t_list)
%
% Returns the boundary between the triangles in t_list and its complement.
% t_list is a list of row indexes i into t(i,j)
% bedges is an m x 2 array listing the edges in the boundary.
% bnodes is a  p x 1 array listing the nodes in the boundary.
% t_idx1(i) is the triangle containing bedges(i) in t_list.
% t_idx2(i) is the triangle containing bedges(i) in the complement of t_list.
t_list = t_list';
% compute complement of t_list
tf = zeros(size(t,1),1);
tf(t_list) = 1;
ct_list = find(tf == 0);
% compute boundary of each part ...
[bedges1,bnodes1,t_idx1a] = boundary2d(t( t_list,:));
[bedges2,bnodes2,t_idx2a] = boundary2d(t(ct_list,:));
% ... and find the common part
temp = sortrows([bedges1, t_list(t_idx1a); 
                 bedges2,ct_list(t_idx2a)+size(t,1)]);
[temp2,idx1] = unique(temp(:,1:2),'rows','first');
[temp2,idx2] = unique(temp(:,1:2),'rows','last');
difflist = find(idx1 ~= idx2);
bedges = temp(idx1(difflist),1:2);
t_idx1 = temp(idx1(difflist),3)';
t_idx2 = temp(idx2(difflist),3)' - size(t,1);
bnodes = unique(sort(bedges(:)));
