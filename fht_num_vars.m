function nvars = fht_num_vars(fht)
% function nvars = fht_num_vars(fht)
%
% Returns the total number of variables for a discretization
% based on the feature hash table (fht), which stores
% variable index lists for each geometric feature.
nvars = 0;
vals = values(fht);
for i = 1:length(vals)
    nvars = nvars + length(vals{i});
end
end
