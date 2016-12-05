function [mb_type,mode_diff,k]=dec_mb_header(k,bits)

if bits(k)=='0'
   mb_type = 0;
   disp('Prediciton block size 16x16')
   k = k + 1;
else
    mb_type = 0;
    disp('Prediciton block size 4x4')
    k = k + 1;
end

% find the mode difference
% [mode_diff,k]= dec_golomb(k,bits,1);

% disp('mode_diff')
% disp(mode_diff)