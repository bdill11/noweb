function ref = get_feature_ref(f,np)
% function ref = get_feature_ref(f,np)
%
% Return a unique integer for the given feature, for
% use in the feature hashtable.
% np is the number of points in the triangulation.
ref = sum(int64(f) .* int64(np).^int64(0:length(f)-1));
end
