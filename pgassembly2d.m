function [A,b] = pgassembly2d(A,b,pde,p,t,elt1,fht1,elt2,fht2,intmethod)
% function [A,b] = pgassembly2d(A,b,pde,p,t,elt1,fht1,elt2,fht2,intmethod)
% 
% Petrov-Galerkin matrix assembly.
% Adds the assembled matrix and vector representing the
% given PDE (pde) to the A matrix & b vector.
% This uses a given elements (elt1, elt2) with the triangulation given by (p,t).
% The feature hash tables (fht1 for elt1, fht2 for elt2) are used to obtain 
% variable indexes for given features. These are obtained by create_fht().
%
% elt1 represents the test functions, while elt2 represents the basis
% functions.
%
% The two elements can be quite independent, but the triangulation must be
% the same for the two sets of variables.
%
% A must be nv1 x nv2 and b must be nv1 x 1 where nv1 is the total
% number of variables for elt1 and nv2 is the total number of variables
% for elt2 (as returned by fht_num_vars()).
% Reference triangle has vertices (0,0), (1,0), (0,1).
[p_int,w_int] = intmethod(); % points and weights for reference triangle
% np is the total number of points in the triangulation
np = size(p,1);
% compute total numbers of variables
nv1 = fht_num_vars(fht1);
nv2 = fht_num_vars(fht2);
% nv_elt is the number of variables in one element
nv_elt1 = sum(elt1.nvars);
nv_elt2 = sum(elt2.nvars);
% order is the order of derivatives used in the assembly;
% we need 0 <= order <= 2
order = pde.order;
intval1 = zeros(nv_elt1,nv_elt2);
intval2 = zeros(nv_elt1,1);
% Save Aphihat() values for all the integration points
% on the reference element
Aphihatvals1 = cell(length(w_int),1);
Aphihatvals2 = cell(length(w_int),1);
for k = 1:length(w_int)
    Aphihatvals1{k} = elt1.Aphihat(p_int(k,:),order);
    Aphihatvals2{k} = elt2.Aphihat(p_int(k,:),order);
end
for i = 1:size(t,1) % for all triangles ...
    % obtain variable list and signs for this triangle
    [vlist1,slist1] = get_var_triangle(t(i,:),fht1,elt1,np);
    [vlist2,slist2] = get_var_triangle(t(i,:),fht2,elt2,np);
    % set up affine transformation xhat :-> x = T.xhat + b
    i1 = t(i,1); i2 = t(i,2); i3 = t(i,3);
    T = [p(i2,:)'-p(i1,:)', p(i3,:)'-p(i1,:)'];
    b0 = p(i1,:)';
    % form weighted sum of integrand at integration points
    intval1 = 0;
    intval2 = 0;
    for k = 1:length(w_int)
        Aphival1 = elt1.trans_Aphihat(T,Aphihatvals1{k},order);
        Aphival2 = elt2.trans_Aphihat(T,Aphihatvals2{k},order);
        Dmat   = pde.coeffs(T*p_int(k,:)'+b0);
        rhsvec = pde.rhs(   T*p_int(k,:)'+b0);
        integrand_val1 = Aphival1*Dmat*Aphival2';
        integrand_val2 = Aphival1*rhsvec;
        intval1 = intval1 + w_int(k)*integrand_val1;
        intval2 = intval2 + w_int(k)*integrand_val2;
    end
    detT = abs(det(T));
    intval1 = intval1*detT; % scale by Jacobian
    intval2 = intval2*detT;
    intval1 = diag(slist1)*intval1*diag(slist2); % change signs if needed
    intval2 = slist1'.*intval2;
    A(vlist1,vlist2) = A(vlist1,vlist2) + intval1; % add to matrix & vec
    b(vlist1)        = b(vlist1)        + intval2;
end % for i
