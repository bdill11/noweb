function [p,t] = ref_triangle_submesh(n)
% function [p,t] = ref_triangle_submesh(n)
%
% Creates a standard mesh on the reference triangle
% (vertices at (0,0), (1,0) and (0,1)).
% n+1 is the number of grid points on each edge
% Generate points
p = zeros(n*(n+1)/2,2);
k = 1;
for i = 0:n
    for j = 0:n-i
        p(k,:) = [i,j];
        k = k+1;
    end
end
p = p / n;
% Create triangulation
t = zeros(n*n,3);
k = 1;
idx = 1;
for i = 0:n-1
    for j = 0:n-i-1
        if j > 0
            t(k,:) = [idx, idx+n-i, idx+n-i+1];
            k = k+1;
        end
        t(k,:) = [idx, idx+n-i+1, idx+1];
        k = k+1;
        idx = idx+1;
    end % for j
    idx = idx+1;
end % for i
