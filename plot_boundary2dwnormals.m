function plot_boundary2dwnormals(p,t,bb)
% function plot_boundary2dwnormals(p,t,bb)
%
% Plots the boundary of a mesh in 2D given by p and t.
% The i'th point is p(i,:), and triangle j is given by
% points with indexes t(j,1), t(j,2) & t(j,3).
% The normals are also plotted from the center of each edge.
% The length of the normal vectors is the length of the edge.
%
% The bounding box is given in bb = [xmin, xmax, ymin,ymax].
%
% See distmesh.m etc.
[bedges,bnodes,normals] = boundary2d_2(p,t);
bdry_tri = [bedges(:,1),bedges(:,2),bedges(:,2)];
midpts = 0.5*(p(bedges(:,1),:)+p(bedges(:,2),:));
len_edges = sqrt(sum((p(bedges(:,1),:)-p(bedges(:,2),:)).^2,2));
arrowpts = midpts + diag(sparse(len_edges))*normals;
hold on
triplot(bdry_tri,p(:,1),p(:,2));
arrow(midpts,arrowpts);
axis(bb);
hold off
