function [pv,tv,vals] = get_submesh_vals(p,t,fht,elt,vars,p_ref,t_ref,order)
% function [pv,tv,vals] = get_submesh_vals(p,t,fht,elt,vars,p_ref,t_ref,order)
%
% Return triangulation (pv,tv) and values (vals) for 
% the given variable values (vars).
% Each triangle in the master triangulation (p,t)
% is subdivided according to the triangulation given for
% the reference element (p_ref,t_ref).
%
% The relationship between the elements of the master
% triangulation (p,t) is given by fht and elt.
%
% Values and derivatives up to the given order are returned in vals.
% Get the values for the basis functions on p_ref
Aphihat = cell(size(p_ref,1),1);
for i = 1:size(p_ref,1)
    Aphihat{i} = elt.Aphihat(p_ref(i,:),order);
end
np = size(p,1);
np_ref = size(p_ref,1);
nt_ref = size(t_ref,1);
pv = zeros(size(t,1)*np_ref,2);
tv = zeros(size(t,1)*nt_ref,3);
vals = zeros(size(t,1)*np_ref,size(Aphihat{1},2));
for i = 1:size(t,1)
    T = [p(t(i,2),:)'-p(t(i,1),:)', p(t(i,3),:)'-p(t(i,1),:)'];
    b = p(t(i,1),:)';
    pv((i-1)*np_ref+(1:np_ref),:) = p_ref*T'+ones(np_ref,1)*b';
    tv((i-1)*nt_ref+(1:nt_ref),:) = t_ref+(i-1)*np_ref;
    % get variable indexes
    [vlist,slist] = get_var_triangle(t(i,:),fht,elt,np);
    % basis function values
    for j = 1:np_ref
        Aphival = elt.trans_Aphihat(T,Aphihat{j},order);
        vals((i-1)*np_ref+j,:) = (vars(vlist)'.*slist)*Aphival;
    end % for j
end % for i
end % function
