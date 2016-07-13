function [vlist,slist] = get_var_triangle(tri,fht,elt,np)
% function [vlist,slist] = get_var_triangle(tri,fht,elt,np)
%
% Get the list of variables (vlist) and the list of sign changes (slist)
% for a given triangle tri using the feature ahstable (fht) for
% the given element type (see elt data structure).
% Note that np is the number of points.
%
% tri is a 1 x 3 array of indexes into the p array of points in the
% triangulation
tri2 = [0,tri];
flist = elt.flist;
flist = tri2(flist+1); % use point indexes
vlist = [];
slist = [];
for i = 1:size(flist,1) % for each feature
f = flist(i,:);
f = f(find(f ~= 0)); % strip zeros from f
[fn,px] = sort(f); % normalize f: fn(px) == f
ref = get_feature_ref(fn,np);
if ~ isKey(fht,ref)
    error('flexPDE:missing value','get_var_triangle: Missing feature',fn,ref)
return
else
    fvlist = fht(ref);
    [pxvars,fslist] = elt.pxfeature(px);
    fvlist = fvlist(pxvars);
end
% concatenate the list of variable indexes & signs
vlist = [vlist,fvlist];
slist = [slist,fslist];
end
