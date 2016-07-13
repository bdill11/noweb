function [p,w] = int2d_radon7()
% function [p,w] = int2d_radon7()
%
% Returns the points (p) and weights (w) of J. Radon's 7-point
% integration formula for the triangle with vertices (0,0), (1,0), (0,1).
% This formula is exact for polynomials up to degree 5.
% Points are the rows of p.
%
% Reference: J. Radon, Zur mechanischen Kubatur. (German)
% Monatsh. Math. 52, (1948), pp. 286-300.
p = [1/3,              1/3;
    (6+sqrt(15))/21,   (9-2*sqrt(15))/21;
    (9-2*sqrt(15))/21, (6+sqrt(15))/21;
    (6+sqrt(15))/21,   (6+sqrt(15))/21;
    (6-sqrt(15))/21,   (9+2*sqrt(15))/21;
    (9+2*sqrt(15))/21, (6-sqrt(15))/21;
    (6-sqrt(15))/21,   (6-sqrt(15))/21];
w = [9/80;
    (155+sqrt(15))/2400;
    (155+sqrt(15))/2400;
    (155+sqrt(15))/2400;
    (155-sqrt(15))/2400;
    (155-sqrt(15))/2400;
    (155-sqrt(15))/2400];
end
