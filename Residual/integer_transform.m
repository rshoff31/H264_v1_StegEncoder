function [W]= integer_transform(X)

% X is a 4x4 block of data
% W is the trasnsformed coefficients

a = 1/4;
b = 1/10;
c = sqrt(1/40);

% E is the scaling factor matrix, to make it similar to DCT ?
% Refer to MPEG4-AVC slides (simplified from the H.264 white paper)

E = [a c a c
     c b c b
     a c a c
     c b c b];

 % C is the core transform matrix
C =  [1 1 1 1
      2 1 -1 -2
      1 -1 -1 1
      1 -2 2 -1];
 
  W = (C*X*C'); 
% W = (C*X*C').*E;  

end