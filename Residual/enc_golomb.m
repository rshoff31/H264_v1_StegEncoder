function [bits] = enc_golomb(symbol, sign)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This script/function was created by 
% Abdullah Al Muhit
% contact - almuhit@gmail.com
% website - https://sites.google.com/site/almuhit/
% Please use it at your own risk. Also, Please cite the following paper:
% A A Muhit, M R Pickering, M R Frater and J F Arnold, “Video Coding using Elastic Motion Model and Larger Blocks,” IEEE Trans. Circ. And Syst. for Video Technology, vol. 20, no. 5, pp. 661-672, 2010. [Impact factor – 3.18] [PDF]
% A A Muhit, M R Pickering, M R Frater and J F Arnold, “Video Coding using Geometry Partitioning and an Elastic Motion Model,” accepted for publication in Journal of Visual Communication and Image Representation. [Impact factor – 1.33] [PDF]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Encodes exponential golomb codes for a SINGLE SYMBOL
% if singed_symbols=1, singed mapping is used
% otherwise unsigned mapping is used

% symbols = -3;
% signed_symbols = 1;

bits = '';

% If signed_symbol flag is 1
if (sign)
    if (symbol ==0)
%         symbol = symbol;
    elseif (symbol>0)
        symbol = 2*symbol -1;
    else 
        symbol = (-2)*symbol;
    end
% if unsigned integers are used    
else
%     symbol = symbol;
end

% Here code_num = symbol
% M is prefix, info is suffix
M = floor(log2(symbol + 1));
info = dec2bin(symbol + 1 - 2^M,M);

for j=1:M
    bits = [bits '0'];
end
bits = [bits '1'];
bits = [bits info];
    
