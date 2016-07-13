function match = match_edge_triangle(edge,triangle)
% function match = match_edge_triangle(edge,triangle)
%
% Returns index vector (2 elements) match so that
% edge(i) = triangle(match(i)), i = 1, 2.
for i = 1:2
    for j = 1:3
        if edge(i) == triangle(j)
            match(i) = j;
        end
    end % for j
end % for i
end % function
