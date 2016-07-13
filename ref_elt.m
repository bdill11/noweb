function [p,t] = ref_elt()
% function [p,t] = ref_elt()
% 
% Returns the triangulation of the reference element.
% This is a trivial task, but provides a convenient
% benchmark for other tasks.
p = [0 0; 1 0; 0 1];
t = [1 2 3];
