function elt = const2d_elt()
% function elt = const2d_elt()
%
% Returns constant 2-D triangle element
nvars  = [1];
flist  = [1 2 3];
vnodes = [1/3, 1/3];
elt = struct('Aphihat',@const2d_Aphihat, ...
    'nvars',nvars,'flist',flist, ...
    'pxfeature',@const2d_pxfeature,'vnodes',vnodes, ...
    'trans_Aphihat',@trans2d_Aphilist);
end % function

function Aphilist = const2d_Aphihat(xhat,order)
% function Aphilist = const2d_Aphihat(xhat,order)
%
% Returns array of basis function values, their gradient and Hessian entries
% for constant basis functions on a 2-D reference triangle at xhat.
% Basis function values
Aphilist0 = [1];
if order >= 1
    % Basis gradient values (along rows)
    Aphilist1 = [0 0];
end
if order >= 2
    % Basis hessian values (along rows: d1^2, d1.d2, d2^2)
    Aphilist2 = [0 0 0];
end
if order == 0
    Aphilist = Aphilist0;
elseif order == 1
    Aphilist = [Aphilist0,Aphilist1];
elseif order == 2
    Aphilist = [Aphilist0,Aphilist1,Aphilist2];
end % if
end % function

function [px_vars,signs] = const2d_pxfeature(px)
% function [px_vars,signs] = const2d_pxfeature(px)
%
% Returns the permutation of the variables (px_vars),
% and the sign changes (signs) resulting from a permutation (px)
% applied to a feature of the appropriate dimension (== length(px)).
% This is for the quadratic 2-D triangle elements.
dim = sum(px ~= 0)-1;
switch dim
    case 2 % triangles
        px_vars = [1]; signs = [1];
    otherwise % not a valid feature
        px_vars = []; signs = [];
end % switch
end % function
