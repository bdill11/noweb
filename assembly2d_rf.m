function [A,b] = assembly2d_rf(A,b,pde,p,t,elt1,fht1,intmethod)
% function [A,b] = assembly2d_rf(A,b,pde,p,t,elt1,fht1,intmethod)
% Adds the assembled matrix and vector representing the
% given PDE (pde) to the A matrix & b vector.
% This uses a given element (elt1) with the triangulation given by (p,t).
% The feature hash table (fht1) is used to obtain variable indexes
% for given features. This is obtained by create_fht().
%
% A must be nv1 x nv1 and b must be nv1 x 1 where nv1 is the total
% number of variables (as returned by fht_num_vars(fht1)).
% Reference triangle has vertices (0,0), (1,0), (0,1).
[p_int,w_int] = intmethod(); % points and weights for reference triangle
% np is the total number of points in the triangulation
np = size(p,1);
% order is the order of derivatives used in the assembly;
% we need 0 <= order <= 2
order = pde.order;
% nv_elt1 is the number of variables in one element
nv_elt1 = sum(elt1.nvars);
update_mat = zeros(nv_elt1,nv_elt1);
update_vec = zeros(nv_elt1,1);
% Save Aphihat() values for all the integration points
% on the reference element
Aphihatvals = cell(length(w_int),1);
for k = 1:length(w_int)
    Aphihatvals{k} = elt1.Aphihat(p_int(k,:),order);
end
for i = 1:size(t,1) % for all triangles ...
    % obtain variable list and signs for this triangle
    [vlist,slist] = get_var_triangle(t(i,:),fht1,elt1,np);
    % set up affine transformation xhat :-> x = T.xhat + b0
    i1 = t(i,1);  i2 = t(i,2);  i3 = t(i,3);
    T = [p(i2,:)'-p(i1,:)', p(i3,:)'-p(i1,:)'];
    b0 = p(i1,:)';

    % form weighted sum of integrand at integration points
    update_mat = 0;
    update_vec = 0;
    for k = 1:length(w_int)
        Aphival = elt1.trans_Aphihat(T,Aphihatvals{k},order);
        Dmat    = pde.coeffs(T*p_int(k,:)'+b0);
        update_mat = update_mat + w_int(k)*(Aphival*Dmat*Aphival');
        rhsvec  = pde.rhs(T*p_int(k,:)'+b0);
        update_vec = update_vec + w_int(k)*(Aphival*rhsvec);
    end
    detT = abs(det(T));
    update_mat = update_mat*detT;                    % scale by Jacobian
    update_mat = diag(slist)*update_mat*diag(slist); % change signs if needed
    A(vlist,vlist) = A(vlist,vlist) + update_mat;    % add to matrix
    update_vec = update_vec*detT;                    % scale by Jacobian
    update_vec = slist'.*update_vec;                 % change signs if needed
    b(vlist)   = b(vlist)           + update_vec;    % add to vector
end
