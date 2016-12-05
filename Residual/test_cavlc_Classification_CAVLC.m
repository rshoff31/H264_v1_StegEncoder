clear all;
clc;
load table;
global Table_coeff0 Table_coeff1 Table_coeff2 Table_coeff3
global Table_run Table_zeros
load 00AES_Bits.mat               % AES_Bits is 128,000 random bits
global AES_count AES_Bits             % AES_Bits(count)= 1 or 0 randomly

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

% data = [128 -1 -3 1
%         -7 -4 -5 2
%         -1 -6 2 4
%          1 -1 20 1]; 

nL = 0;
nU = 0;
AES_count=0;
[bits] = enc_cavlc_Classification_CAVLC(data, nL, nU);
AES_count=0;
[data_rec] = dec_cavlc_Classification_CAVLC(bits,nL,nU);

diff = data - data_rec