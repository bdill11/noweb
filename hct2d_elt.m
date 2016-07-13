function elt = hct2d_elt()
% function elt = hct2d_elt()
%
% Hsieh-Clough-Tocher element in two dimensions.
nvars = [3;3;3;1;1;1];
flist = [1 0 0;  1 0 0;  1 0 0;
         2 0 0;  2 0 0;  2 0 0;
         3 0 0;  3 0 0;  3 0 0;
         1 2 0;
         1 3 0;
         2 3 0];
vnodes = [0 0;  0 0;  0 0;
          1 0;  1 0;  1 0;
          0 1;  0 1;  0 1;
          1/2 0; 0 1/2; 1/2 1/2];
elt = struct('Aphihat',@hct2d_Aphihat, ...
    'nvars',nvars,'flist',flist, ...
    'pxfeature',@hct2d_pxfeature,'vnodes',vnodes, ...
    'trans_Aphihat',@hct2d_trans_Aphilist);
end
end % function
 
