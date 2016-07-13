%fd = @(p)ddiff(drectangle(p,-1,1,-1,1),dcircle(p,0,0,0.5))
fd = @(p)drectangle(p,-1,1,-1,1);
%fh = @(p)min(4*sqrt(sum(p.^2,2))-1,2)
fh = @huniform;
[p,t]=distmesh2d(fd,fh,0.5,[-1,-1;1,1],[-1,-1;-1,1;1,-1;1,1]);
np = size(p,1)

lin2d = lin2d_elt()
fht = create_fht(p,t,lin2d)
nv = fht_num_vars(fht)
f = @(x)(10*x(1)^2*exp(x(2)))
pde = struct('coeffs',@(x)diag([0,1,1]),'rhs',@(x)[f(x);0;0],'order',1)
% Initialize A and b 
A = sparse(nv,nv); 
b = zeros(nv,1); 
% Assemble matrix and vector  
[A,b] = assembly2d(A,b,pde,lin2d,p,t,fht,@int2d_radon7); 
[bedges,bnodes,t_index] = boundary2d(t);
g = @(x)(cos(x(1))*x(2))
pde2 = struct('coeffs',@(x)[1],'rhs',@(x)g(x),'order',0)
[Ab,bb,bvlist] = assembly2dbdry(pde2,lin2d,p,t,bedges,t_index,fht,@int1d_gauss5);
g1 = Ab(bvlist,bvlist) \ bb(bvlist);
% find non-boundary variables (cbvlist) 
v_array = ones(nv,1); v_array(bvlist) = 0; 
cbvlist = find(v_array ~= 0);

u_int = A(cbvlist,cbvlist) \ (b(cbvlist) - A(cbvlist,bvlist)*g1);
u = zeros(nv,1); 
u(cbvlist) = u_int; 
u( bvlist) = g1; 

pvlist = get_pvlist(fht,np); 
figure(2)
trimesh(t,p(:,1),p(:,2),u(pvlist))

% % Now using piecewise quadratic elements
% quad2d = quad2d_elt()
% fht2 = create_fht(p,t,quad2d)
% nv2 = fht_num_vars(fht2)
% % Initialize A2 and b2 
% A2 = sparse(nv2,nv2); 
% b2 = zeros(nv2,1); 
% % Assemble matrix and vector  
% [A2,b2] = assembly2d(A2,b2,pde,quad2d,p,t,fht2,@int2d_radon7); 
% [Ab2,bb2,bvlist2] = ...
%     assembly2dbdry(pde2,quad2d,p,t,bedges,t_index,fht2,@int1d_gauss5);
% % find non-boundary variables (cbvlist2) 
% v_array = ones(nv2,1); v_array(bvlist2) = 0; 
% cbvlist2 = find(v_array ~= 0);
% 
% g2 = Ab2(bvlist2,bvlist2) \ bb2(bvlist2);
% u_int = A2(cbvlist2,cbvlist2) \ (b2(cbvlist2) - A2(cbvlist2,bvlist2)*g2);
% u2 = zeros(nv2,1); 
% u2(cbvlist2) = u_int; 
% u2( bvlist2) = g2; 
% [pr,tr] = ref_triangle_submesh(4);
% [pv,tv,vals] = get_submesh_vals(p,t,fht2,quad2d,u2,pr,tr,0);
% figure(3)
% trimesh(tv,pv(:,1),pv(:,2),vals)
% f = @(x)10;
% %g = @(x)(0.5*cos(x(1)));
% g = @(x)0;
% h = @(x)exp(x(2));
% w = @(x)[1; -2];
% alpha = @(x)1;
% beta  = @(x)10;
% pde  = struct('coeffs',@(x)[0,w(x)'; [0;0],alpha(x)*eye(2)], ...
%          'rhs',@(x)[f(x);0;0], 'order',1)
% pdeb = struct('coeffs',beta, 'rhs',h, 'order',0)
% 
% in_gamma1_bnodes = (sqrt(p(bnodes,1).^2+p(bnodes,2).^2) > 3/4);
% in_gamma1 = zeros(size(p,1),1);
% in_gamma1(bnodes) = in_gamma1_bnodes;
% in_gamma1_bedges = in_gamma1(bedges(:,1)) & in_gamma1(bedges(:,2));
% bedges1  = bedges (find( in_gamma1_bedges),:);
% t_index1 = t_index(find( in_gamma1_bedges));
% bedges2  = bedges (find(~in_gamma1_bedges),:);
% t_index2 = t_index(find(~in_gamma1_bedges));
% figure(4)
% triplot(t,p(:,1),p(:,2),'k') % plot triangulation
% % plot \Gamma_1 in blue and \Gamma_2 in red
% triplot([bedges1(:,1),bedges1(:,2),bedges1(:,2)],p(:,1),p(:,2),'b')
% triplot([bedges2(:,1),bedges2(:,2),bedges2(:,2)],p(:,1),p(:,2),'r')
% intmethod = @int2d_radon7;
% [A,b,bvlist2] = assembly2dbdry(pdeb,lin2d,p,t, ...
%       bedges2,t_index2,fht,@int1d_gauss5);
% [A,b] = assembly2d(A,b,pde,lin2d,p,t,fht,intmethod);
% dir_bc_pde = struct('coeffs',@(x)[1],'rhs',@(x)g(x),'order',0)
% Ab = sparse(nv,nv);
% bb =  zeros(nv,1);
% [Ab,bb,dir_bc_vlist] = ... 
%       assembly2dbdry(dir_bc_pde,lin2d,p,t, ...
%       bedges1,t_index1,fht,@int1d_gauss5);
% u1 = Ab(dir_bc_vlist,dir_bc_vlist) \ bb(dir_bc_vlist);
% % Get complement to dir_bc_vlist
% varray = zeros(nv,1);
% varray(dir_bc_vlist) = 1;
% cvlist = find(varray == 0);
% % Now solve linear system
% u2 = A(cvlist,cvlist) \ (b(cvlist) - A(cvlist,dir_bc_vlist)*u1);
% u = zeros(nv,1);
% u(dir_bc_vlist) = u1;
% u(cvlist)       = u2;
% figure(5)
% trimesh(t,p(:,1),p(:,2),u(pvlist))
