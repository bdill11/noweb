function elt = cub2d_elt()
% function elt = cub2d_elt()
%
% Returns cubic 2-D (10-point) element data structure
nvars = [1;1;1;2;2;2;1];
flist = [1 0 0;
         2 0 0;
         3 0 0;
         1 2 0;
         1 3 0;
         2 3 0;
         1 2 3];
vnodes = [0   0;
          1   0;
          0   1;
          0   1/3;
          0   2/3;
          1/3 0;
          2/3 0;
          1/3 2/3;
          2/3 1/3;
          1/3 1/3];
elt = struct('Aphihat',@cub2d_Aphihat, ...
    'nvars',nvars,'flist',flist, ...
    'pxfeature',@cub2d_pxfeature,'vnodes',vnodes, ...
    'trans_Aphihat',@trans2d_Aphilist);
end % function

function Aphihat = cub2d_Aphihat(xhat,order)
% function Aphihat = cub2d_Aphihat(xhat,order)
%
% Returns basis function values, gradients and Hessian
% entries for Lagrangian cubic basis functions on the
% reference triangle with vertices (0,0), (1,0), and (0,1).
% Each row contains the value, 1st derivatives, and 2nd derivatives 
% of the corresponding basis function on the reference element.
x = xhat(1);
y = xhat(2);
Aphihat0 = [(9/2)*(1-x-y)*(2/3-x-y)*(1/3-x-y);
    (9/2)*x*(x-1/3)*(x-2/3);
    (9/2)*y*(y-1/3)*(y-2/3);
    (27/2)*x*(2/3-x-y)*(1-x-y);
    (27/2)*x*(x-1/3)*(1-x-y);
    (27/2)*y*(2/3-x-y)*(1-x-y);
    (27/2)*y*(y-1/3)*(1-x-y);
    (27/2)*x*y*(x-1/3);
    (27/2)*x*y*(y-1/3);
    27*x*y*(1-x-y)];
if order >= 1
    Aphihat1 = [ ...
      18*x + 18*y - 27*x*y - (27*x^2)/2 - (27*y^2)/2 - 11/2, 18*x + 18*y - 27*x*y - (27*x^2)/2 - (27*y^2)/2 - 11/2;
      (27*x^2)/2 - 9*x + 1, 0;
      0, (27*y^2)/2 - 9*y + 1;
      (81*x^2)/2 + 54*x*y - 45*x + (27*y^2)/2 - (45*y)/2 + 9, (9*x*(6*x + 6*y - 5))/2;
      36*x + (9*y)/2 - 27*x*y - (81*x^2)/2 - 9/2, -(27*x*(x - 1/3))/2;
      (9*y*(6*x + 6*y - 5))/2, (27*x^2)/2 + 54*x*y - (45*x)/2 + (81*y^2)/2 - 45*y + 9;
      -(27*y*(y - 1/3))/2, (9*x)/2 + 36*y - 27*x*y - (81*y^2)/2 - 9/2;
      (9*y*(6*x - 1))/2, (27*x*(x - 1/3))/2;
      (27*y*(y - 1/3))/2, (9*x*(6*y - 1))/2;
      -27*y*(2*x + y - 1), -27*x*(x + 2*y - 1)];
end
if order >= 2
    Aphihat2 = [ ...
      18 - 27*y - 27*x, 18 - 27*y - 27*x, 18 - 27*y - 27*x;
      27*x - 9, 0, 0;
      0, 0, 27*y - 9;
      81*x + 54*y - 45, 54*x + 27*y - 45/2, 27*x;
      36 - 27*y - 81*x, 9/2 - 27*x, 0;
      27*y, 27*x + 54*y - 45/2, 54*x + 81*y - 45;
      0, 9/2 - 27*y, 36 - 81*y - 27*x;
      27*y, 27*x - 9/2, 0;
      0, 27*y - 9/2, 27*x;
      -54*y, 27 - 54*y - 54*x, -54*x];
end
if order == 0
    Aphihat = Aphihat0;
elseif order == 1
    Aphihat = [Aphihat0, Aphihat1];
elseif order == 2
    Aphihat = [Aphihat0, Aphihat1, Aphihat2];
end % if
end % function

function [px_vars,signs] = cub2d_pxfeature(px)
% function [px_vars,signs] = cub2d_pxfeature(px)
%
% Returns the permutation of the variables (px_vars),
% and the sign changes (signs) resulting from a permutation (px)
% applied to a feature of the appropriate dimension (== length(px)).
% This is for the quadratic 2-D triangle elements.
dim = sum(px ~= 0)-1;
switch dim
    case 0 % points
        px_vars = [1]; signs = [1];
    case 1 % edges
        px_vars = px; signs = [1 1];
    case 2 % triangles
        px_vars = [1]; signs = [1];
    otherwise % not a valid feature
        px_vars = []; signs = [];
end % switch
end % function

