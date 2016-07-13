function compmethod = int2d_comp(pr,tr,intmethod)
% function compmethod = int2d_comp(pr,tr,intmethod)
% 
% Returns function handle for composite integration
% method that uses intmethod() as the basic method.
% replicated across all triangles of the triangulation
% (pr,tr) of the reference element.
compmethod = @int2d_comp_func(pr,tr,intmethod);
end function

function [p_int,w_int] = int2d_comp_function(pr,tr,intmethod)
% function [p_int,w_int] = int2d_comp_function(pr,tr,intmethod)
%
% Internal function for int2d_comp().
% This is where the work gets done.
[p_base,w_base] = intmethod(); % base method points & weights
base_len = size(p_base,1);
p_int = zeros(size(tr,1)*base_len,2);
w_int = zeros(size(tr,1)*base_len,1);
for j = 1:size(tr,1)
    % For sub-triangle j...
    % Create affine transformation
    i1 = tr(i,1);  i2 = tr(i,2);  i3 = tr(i,3);
    T = [pr(i2,:)'-pr(i1,:)', pr(i3,:)'-pr(i1,:)'];
    b0 = pr(i1,:)';
    % transform weights and points and add to list
    detT = abs(det(T));
    p_int(((j-1)*base_len+1):(j*base_len),:) = p_base*T'+b0';
    w_int(((j-1)*base_len+1):(j*base_len))   = w_int*detT;
end % for
end % function
