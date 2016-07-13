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
trimesh(t_testref,p_testref(:,1),p_testref(:,2))
axis([-.3 4 -.3 3])
master = [1,2,3,4,5,7,8]
[p_ref,t_ref,master_ref,idx_refref,rfht] = ...
create_refinement(p_testref,t_testref,master,refref2d)
trimesh(t_ref,p_ref(:,1),p_ref(:,2))
'Press any key to continue'
pause

p_refref2 = [0 0; 1 0; 0 1; 2/3 1/3; 1/3 2/3; 0 1/3; 0 2/3; 1/3 0; 2/3 0; 1/3 1/3];
t_refref2 = [1 8 6; 6 8 10; 8 9 10; 9 4 10; 9 2 4; 7 6 10; 10 5 7; 10 4 5; 7 5 3];
npts2 = [ 1; 1; 1; 2; 2; 2; 1];
flist2 = [1 0 0; 2 0 0; 3 0 0; 2 3 0; 1 3 0; 1 2 0; 1 2 3];
px_refref2 = @(px)ifte(length(px)==2,px,[1])
refref2d_2 = struct('p',p_refref2,'t',t_refref2,'npts',npts2,'flist',flist2, ...
'Brefref',[],'pxfeature',px_refref2);
[p_ref2,t_ref2,master_ref2,idx_refref2,rfht2] = ...
create_refinement(p_testref,t_testref,master,refref2d_2)
trimesh(t_ref2,p_ref2(:,1),p_ref2(:,2))
'Press any key to continue'
pause

npts3 = [ 1; 1; 1; 1; 1; 1; 3];
flist3 = [1 0 0; 2 0 0; 3 0 0; 2 3 0; 1 3 0; 1 2 0; 1 2 3];
p_refref3 = [0 0; 1 0; 0 1; 1/2 1/2; 0 1/2; 1/2 0; 1/6 1/6; 2/3 1/6; 1/6 2/3];
t_refref3 = [1 7 6; 6 7 8; 6 8 2; 7 5 1; 7 5 9; 7 8 9; 8 9 4; 2 4 8; 5 9 3; 9 4 3];
size_p_refref3 = size(p_refref3)
size_t_refref3 = size(t_refref3)
px_refref3 = @(px)ifte(length(px)==3,px,[1])
refref2d_3 = struct('p',p_refref3,'t',t_refref3,'npts',npts3,'flist',flist3, ...
'Brefref',[],'pxfeature',px_refref3)
[p_ref3,t_ref3,master_ref3,idx_refref3,rfht3] = ...
create_refinement(p_testref,t_testref,master,refref2d_3)
trimesh(t_ref3,p_ref3(:,1),p_ref3(:,2))
