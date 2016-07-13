function ptlist = rr_get_edgepts(edge,rr)
% function ptlist = rr_get_edgepts(edge,rr)
%
% Return list of points in edge in the
% reference refinement.  Note that edge
% is a pair of indexes into rr.p for the
% reference triangle.
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