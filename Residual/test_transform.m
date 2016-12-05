
clear all;
clc;

load table;
global Table_coeff0 Table_coeff1 Table_coeff2 Table_coeff3
global Table_run Table_zeros


% X = [5 11 8 10
%      9 8 4 12
%      1 10 11 4
%      19 6 15 7];
%  
%  X =[
%     54    42   179   173
%     25   170    39     9
%    210   228   243   206
%     45   132   138   191];
 
X = [255 255 255 255
     255 255 255 255
     255 255 255 255
     255 255 255 255];
 
 X = [255 0 255 0
     0 255 0 255
     255 0 255 0
     0 255 0 255];
 
 X = [0 255 0 255
     255 0 255 0
     0 255 0 255
     255 0 255 0];
 
 QP = 0;
 
 W = integer_transform(X);
 
 Z = quantization(W,QP);
 
 [bits] = enc_cavlc(Z, 0, 0);
 
 [Zi,i] = dec_cavlc(bits,0,0);
 
 diff = Z - Zi;
 
 Wi = inv_quantization(Zi,QP);

 Y = inv_integer_transform(Wi);

 %  post scaling - very important 
 Xi = round(Y/64);
 
 