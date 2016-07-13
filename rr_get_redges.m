function [elist,rrtlist] = rr_get_redges(edge,rr) 
% function [elist,rrtlist] = rr_get_redges(edge,rr) 
% 
% Return list of refined edges in edge in the 
% reference refinement.  Note that edge % is a pair of indexes into rr.p for the 
% reference triangle. 0; 
% First find list of points in refined triangulation 
% in edge from master triangulation. 
ptlist = []; 
npts  = rr.npts; 
cnpts = cumsum(rr.npts); 
for k = 1:size(rr.flist,1)
    f = rr.flist(k,:);
    f = f(find(f));
    idx = subset_scan(f,edge);
    if min(idx) > 0 % then f is a subset of edge
        ptlist = [ptlist, (cnpts(k)-npts(k)+1):cnpts(k)];
    end 
end 
% ptlist 
% Now find edges in reference refinement that are 
% in edge in the master triangulation 
elist = []; 
rrtlist = []; 
for i = 1:size(rr.t)
    % triangle = rr.t(i,:)
    idx = subset_scan(rr.t(i,:),ptlist);
    if sum(idx > 0) == 2
        % 2 of the 3 vertices are in ptlist
        elist = [elist; rr.t(i,find(idx>0))];
        rrtlist = [rrtlist, i];
    end 
end 
elist = sort(elist,2); 
