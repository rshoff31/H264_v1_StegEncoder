function [bits] = header(h,w,QP,Frame_start,Frame_end)

bits = '';
bits = [bits dec2bin(h,8) dec2bin(w,8) dec2bin(QP,8),dec2bin(Frame_start,8) dec2bin(Frame_end,8)];