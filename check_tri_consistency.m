function tf = check_tri_consistency(p,t)
% function tf = check_tri_consistency(p,t)
%
% Check consistency of orientation of a triangulation
% to avoid overlapping triangles
[elist,tidx1,tidx2] = internal_edges(p,t);
for i = 1:size(elist,1)
    [T1,b1] = gen_transform2d(p(t(tidx1(i),:),:));
    [T2,b2] = gen_transform2d(p(t(tidx2(i),:),:));
end
