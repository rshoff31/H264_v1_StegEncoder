function [bits, mode_prev]= mb_header(mb_type,mode,mode_prev)

bits = '';

bits = [bits dec2bin(mb_type)];

bits = [bits enc_golomb(mode - mode_prev, 1)];

mode_prev = mode;