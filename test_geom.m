p = [0 0; 1 0; 0 1; 1 2.5; 1.5 2.5; ...
      1 1; 2 0; 2 1; 2.5 2.5; 3 3; 3 2; 3 1];
t = [1 2 6; 8 11 9; 11 10 9; 7 8 12; 3 1 6; 6 4 3; ...
      2 7 6; 8 11 12; 7 6 8; 4 5 6; 5 8 6; 9 5 8];
[bedges,bnodes,t_index] = boundary2d(t)
bedges_exact  = [1 2; 1 3; 2 7; 3 4; 4 5; 5 9; 7 12; 9 10; 10 11; 11 12]
t_index_exact = [1  ; 5  ; 7  ; 6  ; 10 ; 12 ; 4   ; 3   ; 3    ; 8    ]
bnodes_exact  = [1; 2; 3; 4; 5; 7; 9; 10; 11; 12] % only 6 & 8 are not
quad2d = quad2d_elt()
fht_quad2d = create_fht(p,t,quad2d)
np = size(p,1)
