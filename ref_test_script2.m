% Test script for refinement
p_refref = [0 0; 1 0; 0 1; 1/2 1/2; 0 1/2; 1/2 0];
t_refref = [1 5 6; 2 4 6; 3 4 5; 4 5 6];
npts = [1; 1; 1; 1; 1; 1];
flist = [1 0 0; 2 0 0; 3 0 0; 2 3 0; 1 3 0; 1 2 0];
Brefref = [1 0 0 0 1/2 1/2;
0 1 0 1/2 0 1/2;
0 0 1 1/2 1/2 0];
px_refref1 = @(px)[1];
refref2d = struct('p',p_refref,'t',t_refref,'npts',npts,'flist',flist, ...
'Brefref',Brefref,'pxfeature',px_refref1);
p_testref = [0 0; 1.5 0; 2.5 0; 3.5 0; 3.5 1; ...
2.5 0.7; 1.5 1; 0 1; 1.5 2; 2.5 1.5];
t_testref = [1 2 7; 2 7 6; 2 3 6; 4 3 6; 4 5 6; ...
10 5 6; 6 10 7; 1 7 8; 8 7 9; 7 9 10];
% trimesh(t_testref,p_testref(:,1),p_testref(:,2))
trimesh_labelled(p_testref,t_testref)
axis([-.3 4 -.3 3])
refine_list = [1,2,3,4,5,7,8]
[p_ref,t_ref,master_ref,idx_refref,rfht,master2rr_pt] = ...
create_refinement(p_testref,t_testref,refine_list,refref2d)
% trimesh(t_ref,p_ref(:,1),p_ref(:,2))
trimesh_labelled(p_ref,t_ref)
'Press any key to continue'
pause
fprintf('=============================\n\n')

p_refref2 = [0 0; 1 0; 0 1; 2/3 1/3; 1/3 2/3; 0 1/3; 0 2/3; 1/3 0; 2/3 0; 1/3 1/3];
t_refref2 = [1 8 6; 6 8 10; 8 9 10; 9 4 10; 9 2 4; 7 6 10; 10 5 7; 10 4 5; 7 5 3];
npts2 = [ 1; 1; 1; 2; 2; 2; 1];
flist2 = [1 0 0; 2 0 0; 3 0 0; 2 3 0; 1 3 0; 1 2 0; 1 2 3];
px_refref2 = @(px)ifte(length(px)==2,px,[1])
refref2d_2 = struct('p',p_refref2,'t',t_refref2,'npts',npts2,'flist',flist2, ...
'Brefref',[],'pxfeature',px_refref2);
[p_ref2,t_ref2,master_ref2,idx_refref2,rfht2,master2rr_pt2] = ...
create_refinement(p_testref,t_testref,refine_list,refref2d_2)
% trimesh(t_ref2,p_ref2(:,1),p_ref2(:,2))
trimesh_labelled(p_ref2,t_ref2)
'Press any key to continue'
pause
fprintf('=============================\n\n')

npts3 = [ 1; 1; 1; 1; 1; 1; 3];
flist3 = [1 0 0; 2 0 0; 3 0 0; 2 3 0; 1 3 0; 1 2 0; 1 2 3];
p_refref3 = [0 0; 1 0; 0 1; 1/2 1/2; 0 1/2; 1/2 0; 1/6 1/6; 2/3 1/6; 1/6 2/3];
t_refref3 = [1 7 6; 6 7 8; 6 8 2; 7 5 1; 7 5 9; 7 8 9; 8 9 4; 2 4 8; 5 9 3; 9 4 3];
size_p_refref3 = size(p_refref3)
size_t_refref3 = size(t_refref3)
px_refref3 = @(px)ifte(length(px)==3,px,[1])
refref2d_3 = struct('p',p_refref3,'t',t_refref3,'npts',npts3,'flist',flist3, ...
'Brefref',[],'pxfeature',px_refref3)
[p_ref3,t_ref3,master_ref3,idx_refref3,rfht3,master2rr_pt3] = ...
create_refinement(p_testref,t_testref,refine_list,refref2d_3)
% trimesh(t_ref3,p_ref3(:,1),p_ref3(:,2))
trimesh_labelled(p_ref3,t_ref3)
'Press any key to continue'
pause
fprintf('=============================\n\n')

% Now we take the next step in the process: get the edges and points in the
% refined triangulation associated with the internal boundary between
% refined and unrefined triangles.
rr2 = refref2d_2

% Get edges in the reference refinement for each edge of the reference
% element.
reference_edges = [1 2; 2 3; 1 3]
elists = cell(size(reference_edges,1),1)
rrtlist = cell(size(reference_edges,1),1)
ref_edge_table = zeros(12,1);
for i = 1:3
    [elists{i},rrtlist{i}] = rr_get_redges(reference_edges(i,:),rr2);
    ref_edge_table(get_feature_ref(reference_edges(i,:),3)) = i;
end
fprintf('=============================\n\n')


% Get an element type and create the Brefref matrix for this element type
% using the reference refinement.
elt = lin2d_elt()
ls_rr = struct('coeffs',@(x)[1],'rhs',@(x)elt.Aphihat(x,0)','order',0)
[fht_rr2,v2t_rr2,v2f_rr2,v2fidx_rr2] = create_fht(rr2.p,rr2.t,elt)
[p_re,t_re] = ref_elt()
[fht_re, v2t_re, v2f_re, v2fidx_re]  = create_fht(p_re,t_re,elt)
nvars_rr2 = fht_num_vars(fht_rr2)
A_ls = zeros(nvars_rr2,nvars_rr2);
b_ls = zeros(nvars_rr2,sum(elt.nvars));
size_A_ls = size(A_ls)
size_b_ls = size(b_ls)
[A_ls,b_ls] = assembly2d(A_ls,b_ls,ls_rr,elt,rr2.p,rr2.t,fht_rr2,@int2d_radon7);
Brefref = A_ls \ b_ls
fprintf('=============================\n\n')

% Assemble matrix for Poisson equation on refined triangulation
% but without patching the "hanging nodes"
pde = struct('coeffs',@(x)[0 0 0; 0 1 0; 0 0 1],'rhs',@(x)[cos(x(1)+x(2));0;0],'order',1)
[fht,v2t,v2f,v2fidx] = create_fht(p_ref2,t_ref2,elt);
nvars = fht_num_vars(fht)
A = sparse(nvars,nvars)
b = zeros(nvars,1);
[A,b] = assembly2d(A,b,pde,elt,p_ref2,t_ref2,fht,@int2d_radon7);
nnz_A = nnz(A)
fprintf('=============================\n\n')

% Next, get boundary in the master triangulation between refined &
% unrefined triangles
np_ref = size(p_ref2,1);
[be,bn,t_idx1,t_idx2] = get_internal_boundary2d(t_testref,refine_list)
for i = 1:size(be,1)
    internal_boundary_edge = be(i,:)
    master_edge = subset_scan(be(i,:),t_testref(t_idx1(i),:))
    trx = master2rr_pt2(t_idx1(i),:)
    [master_edge,px] = sort(master_edge)
    edges_refelt = elists{ref_edge_table(get_feature_ref(master_edge,3))}
    edges = trx(edges_refelt)
    pts_refelt = unique(sort([edges_refelt(:,1);edges_refelt(:,2)]))
    pts = trx(pts_refelt)
    tri_refelt = rrtlist{ref_edge_table(get_feature_ref(master_edge,3))}
    tris = trx(rr2.t(tri_refelt',:))
    mvlist = [];
    if isKey(fht,get_feature_ref(sort(be(i,:)),np_ref))
        mvlist = [mvlist,fht(get_feature_ref(sort(be(i,:)),np_ref))];
    end
    mvlist = [mvlist, fht(be(i,1)), fht(be(i,2))]
    revlist = [];
    if isKey(fht_re,get_feature_ref(master_edge,size(p_re,1)))
        revlist = [revlist, fht_re(get_feature_ref(master_edge,size(p_re,1)))];
    end
    revlist = [revlist, fht_re(master_edge(1)), fht_re(master_edge(2))]
    rrvlist = [];
    for j = 1:size(edges_refelt)
        if isKey(fht_rr2,get_feature_ref(edges_refelt(j,:),size(p_re,1)))
            rrvlist = [rrvlist, fht_rr2(get_feature_ref(edges_refelt(j,:),size(p_re,1)))];
        end
    end
    for j = 1:size(pts_refelt)
        if isKey(fht_rr2,pts_refelt(j))
            rrvlist = [rrvlist, fht_rr2(pts_refelt(j))];
        end
    end
    rrvlist
    vlist = [];
    for j = 1:size(edges)
        if isKey(fht,get_feature_ref(edges(j,:),np_ref))
            vlist = [vlist, fht(get_feature_ref(edges(j,:),np_ref))];
        end
    end
    for j = 1:length(pts)
        if isKey(fht,pts(j))
            vlist = [vlist, fht(pts(j))];
        end
    end
    vlist
    fprintf('-------\n')
end
fprintf('=============================\n\n')
