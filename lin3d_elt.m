function elt = lin3d_elt()
% function elt = lin3d_elt()
%
% Returns the linear 3-D (4-point) element data structure.
nvars = [1;1;1;1];
flist = [1 0 0 0;
         2 0 0 0;
         3 0 0 0;
         4 0 0 0];
vnodes = [0 0 0;
          1 0 0;
          0 1 0;
          0 0 1];
elt = struct('Aphihat',@lin3d_Aphihat, ...
    'nvars',nvars,'flist',flist, ...
    'pxfeature',@lin3d_pxfeature,'vnodes',vnodes, ...
    'trans_Aphihat',@trans3d_Aphilist);
end

function Aphilist = lin3d_Aphihat(xhat,order)
% Aphilist = lin3d_Aphihat(xhat,order)
%
% Returns array of basis function values, their gradient and Hessian entries
x = xhat(1);  y = xhat(2);  z = xhat(3);
% Basis function values
Aphilist0 = [1-x-y-z; 
             x; 
             y;
             z];
if order >= 1
    % Basis gradient values (along rows)
    Aphilist1 = [-1 -1 -1;
                  1  0  0;
                  0  1  0;
                  0  0  1];
end
if order >= 2
    % Basis Hessian values (along rows: dx1^2, dx1.dx2, dx1.dx3, dx2^2, dx2.dx3, dx3^2)
    Aphilist2 = [0 0 0 0 0 0;
                 0 0 0 0 0 0;
                 0 0 0 0 0 0];
end
if order == 0
    Aphilist = Aphilist0;
elseif order == 1
    Aphilist = [Aphilist0,Aphilist1];
elseif order == 2
    Aphilist = [Aphilist0,Aphilist1,Aphilist2];
end % if
end % function

function [px_vars,signs] = lin3d_pxfeature(px)
% function [px_vars,signs] = lin3d_pxfeature(px)
%
% Returns the permutation of the variables (px_vars),
% and the sign changes (signs) resulting from a permutation (px)
% applied to a feature of the appropriate dimension (== length(px)).
% This is for the linear (or affine) 3-D tetrahedral elements.
dimp1 = sum(px ~= 0);  % dimp1 == dimension plus 1
switch dimp1
    case 1 % points
        px_vars = [1]; signs = [1];
    otherwise % not a valid feature
        px_vars = []; signs = [];
end % switch
end % function

