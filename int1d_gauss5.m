function [p,w] = int1d_gauss5()
% function [p,w] = int1d_gauss5()
%
% Returns points and weights for a 5-point Gauss rule
% in 1-D for the interval [0,1]
p = [1/2; (1+sqrt(5-2*sqrt(10/7))/3)/2; (1-sqrt(5-2*sqrt(10/7))/3)/2; ...
     (1+sqrt(5+2*sqrt(10/7))/3)/2; (1-sqrt(5+2*sqrt(10/7))/3)/2];
w = 0.5*[128/225; (322+13*sqrt(70))/900; (322+13*sqrt(70))/900; ...
     (322-13*sqrt(70))/900; (322-13*sqrt(70))/900];
