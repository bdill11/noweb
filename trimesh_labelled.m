function trimesh_labelled(p,t)
% function trimesh_labelled(p,t)
%
% Produces plot of triangulation, labelled by vertex number
trimesh(t,p(:,1),p(:,2));
for i = 1:size(p,1)
    text(p(i,1),p(i,2),num2str(i));
end
% put triangle numbers
for i = 1:size(t,1)
    text(sum(p(t(i,:),1))/3,sum(p(t(i,:),2))/3,['[',num2str(i),']'])
end
