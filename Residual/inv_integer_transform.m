function [Y] = inv_integer_transform(W)

a = 1/4;
b = 1/10;
c = sqrt(1/40);

% E is the scaling factor matrix
% Refer to MPEG4-AVC slides (simplified from the H.264 white paper)

E = [a c a c
     c b c b
     a c a c
     c b c b];

 % Ci is the inverse core transform matrix
Ci =  [1 1 1 1
      1 1/2 -1/2 -1
      1 -1 -1 1
      1/2 -1 1 -1/2];

 Y = Ci'*W*Ci;
%  Y = Ci'*(W.*E)*Ci;
 
end