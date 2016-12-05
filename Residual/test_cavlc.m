clear all;
clc;

load table;
global Table_coeff0 Table_coeff1 Table_coeff2 Table_coeff3
global Table_run Table_zeros

data = [0 3 -1 0
        0 -1 1 0
        1 0 0 0
        0 0 0 0];
% data = [-2 4 0 -1
%         3 0 0 0
%         -3 0 0 0
%         0 0 0 0];
% data = [0 0 1 0
%         0 0 0 0
%         1 0 0 1
%         -1 0 0 0];

% data = [200 -1 -3 1
%         -7 -4 -5 2
%         -1 -6 2 4
%          1 -1 20 1]; 

nL = 0;
nU = 0;

[bits] = enc_cavlc(data, nL, nU);

[data_rec] = dec_cavlc(bits,nL,nU);

x = data - data_rec