function pvlist = get_pvlist(fht,np)
% function pvlist = get_pvlist(fht,np)
%
% Get list of variable indexes for the points.
% That is, pvlist(i) is the variable number for the point
% in the triangulation with index i.
% As usual, np is the number of points.
pvlist = zeros(np,1);
for i = 1:np
    pvlist(i) = fht(get_feature_ref(i,np));
end
