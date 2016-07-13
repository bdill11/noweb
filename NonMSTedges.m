function [edges]=NonMSTedges(B,MSTedges)

%function NonMSTedges gives the edges not in the minimum spanning tree of a
%graph.
%input: Matrix B=node-incidence matrix , MSTedges found from function MST


edges= setdiff(1:size(B,2),MSTedges);


end
