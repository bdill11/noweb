function [p_ref,t_ref,master_ref,idx_ref,rfht,master2rr_pt] = ...
create_refinement(p,t,refine_list,refref)
% function [p_ref,t_ref,master_ref,idx_ref,rfht,master2rr_pt] = ...
% create_refinement(p,t,refine_list,refref)
%
% Create refined mesh using reference refinement refref.
% (p,t) represent the unrefined "master" triangulation,
% where refine_list is the list of triangle indexes that are
% to be refined: the triangles to be refined are
% t(refine_list(i),:), i = 1, ..., length(refine_list).
% The returned items are:
% (p_ref,t_ref) are the triangulation of the refined mesh
% master_ref(i) is the row index into t of the master triangle
% of triangle t_ref(i,:)
% idx_ref(i) is the row index into refref.t identifies
% which sub-triangle of triangle master_ref(i) is t_ref(i,:).
% Note that idx_ref(i) == 0 iff t_ref(i,:) is an unrefined triangle.
%
% master2rr_pt(i,:) is the ordered list of indexes of points in
% master triangle i that are in the refined mesh.
% master2rr_pt(i,j) is the index into p_ref of master triangle i
% and point refref.p(j,:).

np = size(p,1);
rfht = containers.Map('KeyType','int64','ValueType','any');
rr_npts = refref.npts;
rr_flist = refref.flist;
p_ref = p; % include all points of unrefined triangles
t_ref = [];
master_ref = [];
nt_ref = 0;
np_ref = size(p,1);
ismaster = zeros(size(t,1),1);
ismaster(refine_list) = 1;
master2rr_pt = zeros(size(t,1),size(refref.p,1));
p_rr_trx = zeros(size(refref.p,1),1);
idx_ref = [];
for i = 1:size(t,1)
    % master_triangle_i = t(i,:)
    % master_triangle_points = [p(t(i,1),:); p(t(i,2),:); p(t(i,3),:)] 
    p_rr_trx(:) = 0;
    if ismaster(i)
        % vertices of master triangle
        p1 = p(t(i,1),:); p2 = p(t(i,2),:); p3 = p(t(i,3),:);
        % transformation from ref element to master triangle
        T = [p2'-p1', p3'-p1'];
        b = p1';
        % transform points in reference refinement
        % Tp_rr = zeros(size(refref.p));
        % for j = 1:size(refref.p,1)
        % Tp_rr(j,:) = refref.p(j,:)*T' + b';
        % end
        Tp_rr = bsxfun(@plus,refref.p*T',b');
        p_rr_idx = 1; % index into refref.p
        % for each geometric feature of the master triangle ...
        for k = 1:size(rr_flist,1) 
            npts_k = refref.npts(k);
            if length(find(rr_flist(k,:))) ~= 1 % not a vertex
                % if geometric feature not in hash table, add the points
                f = rr_flist(k,:);
                f = f(find(f));
                [f,px] = sort(t(i,f));
                fref = get_feature_ref(f,np);
                pt_px = refref.pxfeature(px);
                if ~ isKey(rfht,fref)
                    pt_list = np_ref + (1:npts_k);
                    rfht(fref) = pt_list;
                    p_ref = [p_ref; Tp_rr(p_rr_idx-1+pt_px,:)];
                    np_ref = np_ref + npts_k;
                else % isKey(rfht,fref)
                    pt_list = rfht(fref);
                end % if
                p_rr_trx(p_rr_idx-1+(1:npts_k)) = pt_list(pt_px);
            else   % is a vertex of a master triangle, so use the corresponding 
                % vertex of the master triangle
                p_rr_trx(p_rr_idx) = t(i,rr_flist(k,1));
            end % if
            p_rr_idx = p_rr_idx + npts_k;
        end % for k
        % p_rr_trx
        % refined_points = p_ref(p_rr_trx,:)
        % now add the refined triangles to the mesh
        % new_t_ref = p_rr_trx(refref.t)
        t_ref = [t_ref; p_rr_trx(refref.t)];
        idx_ref = [idx_ref; (1:size(refref.t,1))'];
        nt_ref = nt_ref + size(refref.t,1);
        master_ref = [master_ref, i*ones(1,size(refref.t,1))];
        master2rr_pt(i,:) = p_rr_trx'; % turn p_rr_trx into row vector
    else % ~ ismaster(i)
        % no extra points, just one extra triangle
        nt_ref = nt_ref + 1;
        t_ref = [t_ref; t(i,:)];
        master_ref(nt_ref) = i;
        idx_ref(nt_ref) = 0; % not a refined triangle
        master2rr_pt(i,1:3) = t(i,:);
    end % if ismaster(i)
end % for i
