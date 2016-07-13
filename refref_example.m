function rr = refref_example()
% function rr = refref_example()
%
% Returns a reference refinement (rr)
% which is just a simple division of the
% reference (0,0, (1,0), (0,1) triangle
% into four congruent parts.
p_refref = [0 0; 1 0; 0 1; 1/2 1/2; 0 1/2; 1/2 0];
t_refref = [1 5 6; 2 4 6; 3 4 5; 4 5 6];
npts = [1; 1; 1; 1; 1; 1]; 
% number of points in each geometric feature
flist = [1 0 0; 2 0 0; 3 0 0; 2 3 0; 1 3 0; 1 2 0];
Brefref = [1 0 0 0 1/2 1/2;
    0 1 0 1/2 0 1/2;
    0 0 1 1/2 1/2 0];
rr = struct('p',p_refref,'t',t_refref,'npts',npts,'flist',flist, ...
    'Brefref',Brefref,'pxfeature',@(px)1);
end % function
