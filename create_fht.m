function [fht,v2tnum,v2fnum,v2fidx] = create_fht(p,t,elt)
% function [fht,v2tnum,v2fnum,v2fidx] = create_fht(p,t,elt)
%
% Create feature hash table (fht) which shows what variables
% are associated with which geometric features.
% The geometric features inserted into fht must be
% in normalized form (that is, a sorted vector of point indexes).
%
% This routine also returns a variable-to-triangle-number array (v2tnum)
% a variable-to-feature-number array (v2fnum), and a variable-to-feature-index
% array (v2fidx).
% The variable (or basis function) with global index k is the v2fidx(k)'th
% basis function associated with the v2fnum(k)'th geometric feature of 
% triangle v2tnum(k).
fht = containers.Map('KeyType','int64','ValueType','any');
nvars = elt.nvars;
flist = elt.flist; % list of features with associated variables
v2tnum = [];
v2fnum = [];
v2fidx = [];
% flist is assumed normalized except for trailing zeros
np = size(p,1);
counter = 0;
for i = 1:size(t,1) % for each triangle ...
    triangle = [0, t(i,:)];
    tflist = triangle(flist+1);
    for j = 1:size(flist,1)
        % for each feature ...
        f = tflist(j,:);
        f = f(find(f ~= 0));
        f = sort(f);
        % Is this feature already in fht? If not add its variables.
        ref = get_feature_ref(f,np);
        if ~ isKey(fht,ref)
            fht(ref) = [(counter+1):(counter + nvars(j))];
            counter = counter + nvars(j);
            v2tnum = [v2tnum, i*ones(1,nvars(j))];
            v2fnum = [v2fnum, j*ones(1,nvars(j))];
            v2fidx = [v2fidx, 1:nvars(j)];
        end
    end % for each feature
    size_fht = size(fht);
end % for each triangle
end % function create_fht
