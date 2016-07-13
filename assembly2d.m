function [A,b] = assembly2d(A,b,pde,elt,p,t,fht,intmethod)
% function [A,b] = assembly2d(A,b,pde,elt,p,t,fht,intmethod)
% Adds the assembled matrix and vector representing the
% given PDE (pde) to the A matrix & b vector.
% This uses a given element (elt) with the triangulation given by (p,t).
% The feature hash table (fht) is used to obtain variable indexes
% for given features. This is obtained by create_fht().
%
% A must be nv x nv and b must be nv x 1 where nv is the total
% number of variables (as returned by fht_num_vars()).
% Reference triangle has vertices (0,0), (1,0), (0,1).
[p_int,w_int] = intmethod(); % points and weights for reference triangle
% np is the total number of points in the triangulation
np = size(p,1);
% compute nv = total number of variables
nv = fht_num_vars(fht);
% nv_elt is the number of variables in one element
nv_elt = sum(elt.nvars);
% order is the order of derivatives used in the assembly;
% we need 0 <= order <= 2
order = pde.order;
intval1 = zeros(nv_elt,nv_elt);
% intval2 = zero(nv_elt,1);
intval2 = zeros(nv_elt,size(b,2));
% Save Aphihat() values for all the integration points
% on the reference element
Aphihatvals = cell(length(w_int),1);
for k = 1:length(w_int)
    Aphihatvals{k} = elt.Aphihat(p_int(k,:),order);
end
for i = 1:size(t,1) % for all triangles ...
    % obtain variable list and signs for this triangle
    [vlist,slist] = get_var_triangle(t(i,:),fht,elt,np);
    % set up affine transformation xhat :-> x = T.xhat + b0
    i1 = t(i,1);  i2 = t(i,2);  i3 = t(i,3);
    T = [p(i2,:)'-p(i1,:)', p(i3,:)'-p(i1,:)'];
    b0 = p(i1,:)';
    % form weighted sum of integrand at integration points
    intval1 = 0;
    intval2 = 0;
    for k = 1:length(w_int)
        Aphival = elt.trans_Aphihat(T,Aphihatvals{k},order);
        Dmat    = pde.coeffs(T*p_int(k,:)'+b0);
        rhsvec  = pde.rhs(T*p_int(k,:)'+b0);
        integrand_val1 = Aphival*Dmat*Aphival';
        integrand_val2 = Aphival*rhsvec;
        intval1 = intval1 + w_int(k)*integrand_val1;
        intval2 = intval2 + w_int(k)*integrand_val2;
    end
    detT = abs(det(T));
    intval1 = intval1*detT; % scale by Jacobian
    intval2 = intval2*detT;
    intval1 = diag(slist)*intval1*diag(slist); % change signs if needed
    % intval2 = slist'.*intval2;
    intval2 = bsxfun(@times,intval2,slist');
    A(vlist,vlist) = A(vlist,vlist) + intval1; % add to matrix & vec
    % b(vlist) = b(vlist) + intval2;
    b(vlist,:) = b(vlist,:) + intval2;
end
