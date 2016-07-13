function plot_boundary2d(p,t,bb)
% function plot_boundary2d(p,t,bb)
%
% Plots the boundary of a mesh in 2D given by p and t.
% The i'th point is p(i,:), and triangle j is given by
% points with indexes t(j,1), t(j,2) & t(j,3).
%
% The bounding box is given in bb = [xmin, xmax, ymin,ymax].
%
% See distmesh.m etc.
[bedges,bnodes] = boundary2d(t);
bdry_tri = [bedges(:,1),bedges(:,2),bedges(:,2)];
triplot(bdry_tri,p(:,1),p(:,2));
axis(bb)
