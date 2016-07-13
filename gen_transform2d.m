function [T,b] = gen_transform2d(p)
% function [T,b] = gen_transform2d(p)
%
% Returns T, b so that x \mapsto T*x+b maps
% the reference triangle co{(0,0),(1,0),(0,1)} to
% the given triangle, mapping (0,0) to p(1,:)',
% (1,0) to p(2,:)' and (0,1) to p(3,:)'
b = p(1,:)';
T = [p(2,:)' - p(1,:)', p(3,:)' - p(1,:)'];
