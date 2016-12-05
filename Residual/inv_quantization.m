function [Wi]= inv_quantization(Z,QP)

% q is qbits
q = 15 + floor(QP/6);

% The scaling factor matrix V depend on the QP and the position of the
% coefficient.
%   delta lambda miu
SM = [10 16 13
      11 18 14
      13 20 16
      14 23 18
      16 25 20
      18 29 23];
 
 x = rem(QP,6);
 
 % find delta, lambda and miu values
 d = SM(x+1,1);
 l = SM(x+1,2);
 m = SM(x+1,3);

 V = [d m d m
      m l m l
      d m d m
      m l m l];
  
 % find the inverse quantized coefficients
  Wi = Z.*V;
  Wi = bitshift(Wi,q-15,'int64');
 
end