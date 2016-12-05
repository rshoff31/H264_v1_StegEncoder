function [Z] = quantization(W,QP)

% q is qbits
q = 15 + floor(QP/6);

% M is the multiplying factor which is found from QP value
% MF is the multiplying factor matrix
% rem(QP,6) alpha   beta    gamma
%           (a)     (b)      (g)
% 0         13107   5243    8066
% 1         11916   4660    7490
% 2         10082   4194    6554
% 3         9362    3647    5825
% 4         8192    3355    5243
% 5         7282    2893    4559

MF =[13107 5243 8066
     11916 4660 7490
     10082 4194 6554
     9362  3647 5825
     8192  3355 5243
     7282  2893 4559];
 
x = rem(QP,6);
 
a = MF(x+1,1);
b = MF(x+1,2);
g = MF(x+1,3);

M = [a g a g
     g b g b
     a g a g
     g b g b];

% scaling and quantization 
Z = round(W.*(M/2^q));

% Z = W.*M;
% Z = round(bitshift(Z,-q));



 
end