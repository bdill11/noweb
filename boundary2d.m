function [bedges,bnodes,t_index] = boundary2d(t)
% function [bedges,bnodes,t_index] = boundary2d(t)
%
% Construct boundary edge list from triangle list t
% t is ntriangles x 3, bedges = nedges x 2
% Edge k joins points bd(k,1) and bd(k,2).
%
% Simply check when edges only appear once in the triangle list.
%
% Also returns the triangle index for each boundary edge
t = sort(t,2); % sort each row of t
bd1 = sortrows([t(:,1),t(:,2),(1:size(t,1))';
                t(:,2),t(:,3),(1:size(t,1))';
                t(:,1),t(:,3),(1:size(t,1))']);
[bd2,idx1] = unique(bd1(:,1:2),'rows','first');
[bd2,idx2] = unique(bd1(:,1:2),'rows','last');
eqlist = find(idx1 == idx2);
bedges = bd1(idx1(eqlist),1:2);
t_index = bd1(idx1(eqlist),3);
bnodes = unique(sort(bedges(:)));
