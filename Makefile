files1 = usage.m \
	 assembly2d.m \
	 pgassembly2d.m \
	 assembly2d-rf.m \
	 assembly2d-nl.m \
	 assembly2dbdry.m \
	 get-feature-ref.m \
	 create-fht.m \
	 get-var-triangle.m \
	 get-var-edge.m \
	 fht-num-vars.m \
	 boundary2d.m \
	 match-edge-triangle.m \
	 ref-elt.m \
	 gen-transform2d.m \
	 lin2d-elt.m \
	 trans2d-Aphilist.m \
	 quad2d-elt.m \
	 cub2d-elt.m \
	 const2d-elt.m \
	 eltx2-elt.m \
	 abf2d-elt.m \
	 hct2d-elt.m \
	 lin3d-elt.m \
	 trans3d-Aphilist.m \
	 int2d-centroid1.m \
	 int2d-radon7.m \
	 int2d-gatermann12.m \
	 int2d-dunavant33.m \
	 int1d-gauss5.m \
	 int3d-centroid1.m \
	 int2d-comp.m \
	 create-refinement.m \
	 ref-test-script.m \
	 get-internal-boundary2d.m \
	 rr-get-redges.m \
	 subset-scan.m \
	 get-pvlist.m \
	 plot-boundary2d.m \
	 plot-boundary2dwnormals.m \
	 trimesh-labelled.m \
	 ref-triangle-submesh.m \
	 get-submesh-vals.m \
	 ifte.m \
	 iftev.m \
	 filenamehack.bash \
	 lyx2code.bash \
	 check-derivs.m \
	 test-geom.m \
	 dummy.txt 
files2 = pde-code.tex filelist
files = $(files1) $(files2)
source = pde-code.nw
all: $(files)
$(files): $(source)
	notangle -R$@ $(source) > $@
pde-code.tex: $(source)
	noweave -delay -index $(source) > $@
