function elt = abf2d_elt()
% function elt = abf2d_elt()
%
% ABF vector element for velocity field.
% This is the ABF scalar element "x2".
% The rationale for using this set of basis functions is
% given in Arnold, Brezzi and Fortin, 
%
elt = eltx2_elt(abfs2d_elt());
end % function

function elt = abfs2d_elt()
% function elt = abfs2d_elt()
%
% Returns the ABF scalar 2-D (3-point) element data structure.
% The basis functions for this element are the same as for
% the piecewise linear element, plus the "bubble" function
% \phi_4(x,y) = 27xy(1-x-y). 
nvars = [1;1;1;1]
flist = [1 0 0
         2 0 0
         3 0 0
         1 2 3];
elt = struct('Aphihat',@abfs2d_Aphihat, ...
    'nvars',nvars,'flist',flist, ...
    'pxfeature',@abfs2d_pxfeature,'vnodes',abfs2d_vnodes(), ...
    'trans_Aphihat',@trans2d_Aphilist);
end

function Aphilist = abfs2d_Aphihat(xhat,order)
% Aphilist = abfs2d_Aphihat(xhat,order)
%
% ABF scalar element: linear basis functions plus "bubble" function
% \lambda_1\lambda_2\lambda_3 in barycentric coordinates.
x = xhat(1);  y = xhat(2);
% Basis function values
Aphilist0 = [1-x-y; 
             x; 
             y
             27*x*y*(1-x-y)];
if order >= 1
    % Basis gradient values (along rows)
    Aphilist1 = [-1 -1;
                  1  0;
                  0  1;
                 27*y*(1-y-2*x), 27*x*(1-x-2*y)];
end
if order >= 2
    % Basis hessian values (along rows: dx1^2, dx1.dx2, dx2^2)
    Aphilist2 = [0 0 0;
                 0 0 0;
                 0 0 0;
                 -54*y, 27 - 54*y - 54*x, -54*x];
end
if order == 0
    Aphilist = Aphilist0;
elseif order == 1
    Aphilist = [Aphilist0,Aphilist1];
elseif order == 2
    Aphilist = [Aphilist0,Aphilist1,Aphilist2];
end % if
end % function

function [px_vars,signs] = abfs2d_pxfeature(px)
% function [px_vars,signs] = abfs2d_pxfeature(px)
%
% Returns the permutation of the variables (px_vars),
% and the sign changes (signs) resulting from a permutation (px)
% applied to a feature of the appropriate dimension (== length(px)).
% This is for the linear (or affine) 2-D triangle elements.
dimp1 = sum(px ~= 0);  % dimp1 == dimension plus 1
switch dimp1
    case 1 % points
        px_vars = [1]; signs = [1];
    case 3 % triangles
        px_vars = [1]; signs = [1];
    otherwise % not a valid feature
        px_vars = []; signs = [];
end % switch
end % function

function vnodes = abfs2d_vnodes()
% function vnodes = abfs2d_vnodes()
% 
% Returns the positions of the variable nodes
% with respect to the reference element.
% Same format as p (2 x m)
vnodes = [0   0;
          1   0;
          0   1;
          1/3 1/3];
end % function
