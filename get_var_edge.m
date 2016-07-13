function [vlist] = get_var_edge(edge,fht,np)
% function [vlist] = get_var_edge(edge,fht,np)
%
% Returns list of variable indexes for given edge (including end-point
% variables). The feature hash table (fht) is used to look up variable
% lists. Also np is the number of points in the triangulation.
vlist = [];
ref = get_feature_ref(sort(edge),np);
if isKey(fht,ref)
    vlist = [vlist, fht(ref)];
end
ref = get_feature_ref(edge(1),np);
if isKey(fht,ref)
    vlist = [vlist, fht(ref)];
end
ref = get_feature_ref(edge(2),np);
if isKey(fht,ref)
    vlist = [vlist, fht(ref)];
end
