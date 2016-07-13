function eltx2 = eltx2_elt(elt)
% function eltx2 = eltx2_elt(elt)
%
% Returns element data structure with the same basis
% functions as defined by elt, but with 2 components for
% each component of elt.
nvars = elt.nvars;
nvarsx2 = 2*nvars;
flist = elt.flist;
flistx2 = flist;
trans_Aphihat = @(T,Aphilist,order)(transx2(elt.trans_Aphihat,T,Aphilist,order));
pxfeature = @(px)(pxfeaturex2(elt.pxfeature,px));
vnodes = elt.vnodes();
vnodesx2 = zeros(2*size(vnodes,1),size(vnodes,2));
vnodesx2(2*(1:size(vnodes,1))-1,:) = vnodes;
vnodesx2(2*(1:size(vnodes,1)) ,:) = vnodes;
eltx2 = struct('Aphihat',elt.Aphihat, ...
    'nvars',nvarsx2,'flist',flistx2, ...
    'pxfeature',pxfeature,'vnodes',vnodesx2, ...
    'trans_Aphihat',trans_Aphihat);
end % function

function [px_varsx2,signsx2] = pxfeaturex2(base_pxfeature,px)
% function [px_varsx2,signsx2] = pxfeaturex2(base_pxfeature,px)
%
% Uses base_pxfeature to create pxfeature() function for the "x2" element
[px_vars,signs] = base_pxfeature(px);
px_varsx2 = zeros(1,2*length(px_vars));
signsx2   = zeros(1,2*length(px_vars));
px_varsx2(2*(1:length(px_vars))-1) = 2*px_vars-1;
signsx2(2*(1:length(px_vars))-1)   = signs;
px_varsx2(2*(1:length(px_vars)) )  = 2*px_vars;
signsx2(2*(1:length(px_vars)) )    = signs;
end % function

function Aphilistx2 = transx2(base_trans,T,Aphilist,order)
% function Aphilistx2 = transx2(base_trans,T,Aphilist,order)
%
% Uses base_trans to create trans_Aphilist function for the "x2" element
Aphilist = base_trans(T,Aphilist,order);
nb = size(Aphilist,1);
Aphilistx2 = zeros(2*nb,size(Aphilist,2));
Aphilistx2(2*(1:nb)-1,1) = Aphilist(:,1);
Aphilistx2(2*(1:nb) ,2) = Aphilist(:,1);
if order >= 1
    Aphilistx2(2*(1:nb)-1,3:4) = Aphilist(:,2:3);
    Aphilistx2(2*(1:nb) ,5:6) = Aphilist(:,2:3);
end
if order >= 2
    Aphilistx2(2*(1:nb)-1,7:9) = Aphilist(:,4:6);
    Aphilistx2(2*(1:nb) ,10:12) = Aphilist(:,4:6);
end
end % function

