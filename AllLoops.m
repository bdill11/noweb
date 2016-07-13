function [LoopMatrix]= AllLoops(B,MSTedges,root,distlist)


%function AllLoops(nonMSTedges)
nonMSTedges=NonMSTedges(B,MSTedges);
%LoopMatrix=zeros(length(nonMSTedges),size(B,2));
LoopMatrix=[];
for i=1:length(nonMSTedges)
    loop=MSTgenloop(B,MSTedges,nonMSTedges(i),root,distlist); 
    LoopMatrix=[LoopMatrix ; loop];
  
end

end
