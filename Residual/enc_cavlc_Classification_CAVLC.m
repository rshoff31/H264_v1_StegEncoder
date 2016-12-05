function [bits] = enc_cavlc_Classification_CAVLC(data, nL, nU)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This script/function was created by
% Abdullah Al Muhit
% contact - almuhit@gmail.com
% website - https://sites.google.com/site/almuhit/
% Please use it at your own risk. Also, Please cite the following paper:
% A A Muhit, M R Pickering, M R Frater and J F Arnold, “Video Coding using Elastic Motion Model and Larger Blocks,” IEEE Trans. Circ. And Syst. for Video Technology, vol. 20, no. 5, pp. 661-672, 2010. [Impact factor – 3.18] [PDF]
% A A Muhit, M R Pickering, M R Frater and J F Arnold, “Video Coding using Geometry Partitioning and an Elastic Motion Model,” accepted for publication in Journal of Visual Communication and Image Representation. [Impact factor – 1.33] [PDF]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% CAVLC Encoder.
% takes in 4x4 block of residual data and produces output bits

% load table.mat;
global Table_coeff0 Table_coeff1 Table_coeff2 Table_coeff3
global Table_run Table_zeros

global QP AES_count AES_Bits
bits = '';

% Convert 4x4 matrix data into a 1x16 data of zig-zag scan
[row, col] = size(data);

% check the correct size of the block
if((row~=4)||(col~=4))
    disp('Residual block size mismatch - exit from CAVLC')
    return;
end

scan = [1,1;1,2;2,1;3,1;2,2;1,3;1,4;2,3;3,2;4,1;4,2;3,3;2,4;3,4;4,3;4,4];
for i=1:16
    m=scan(i,1);
    n=scan(i,2);
    l(i)=data(m,n); % l contains the reordered data
end

i_last = 16;
% find the last non-zero co-eff in reverse order
while ((i_last>0)&(l(i_last)==0))
    i_last = i_last - 1;
end

i_total = 0; % Total non zero coefficients
i_total_zero = 0; % Total zeros
i_trailing = 0;
sign = ''; % find sign for trailing ones
idx = 1;

%% find level, trailing ones(sign), run and total zero values
while ((i_last>0)&(abs(l(i_last))==1)&(i_trailing<3))
    level(idx)=l(i_last);
    i_total = i_total + 1;
    i_trailing = i_trailing + 1;
    
    if l(i_last)==-1
        sign = [sign '1'];
    else
        sign = [sign '0'];
    end
    
    run(idx) = 0;
    i_last = i_last - 1;
    while ((i_last>0)&(l(i_last)==0))
        run(idx) = run(idx) + 1;
        i_total_zero = i_total_zero + 1;
        i_last = i_last - 1;
    end
    idx = idx + 1;
    
end

while (i_last>0)
    level(idx)=l(i_last);
    i_total = i_total + 1;
    
    
    run(idx) = 0;
    i_last = i_last - 1;
    while ((i_last>0)&(l(i_last)==0))
        run(idx) = run(idx) + 1;
        i_total_zero = i_total_zero + 1;
        i_last = i_last - 1;
    end
    idx = idx + 1;
    
end

%% Write coeff_token

% find n parameter (context adaptive)
if (nL>0)&(nU>0)
    n = (nL + nU)/2;
elseif (nL>0)|(nU>0)
    n = nL + nU;
else
    n = 0;
end

% Coeff_token mapping
% Rows are the total coefficient(0-16) and columns are the trailing ones(0-3)
% TABLE_COEFF0,1,2,3 ARE STORED IN TABLE.MAT OR CAVLC_TABLES.M FILE
% Choose proper table_coeff based on n value
if 0<=n<2
    Table_coeff = Table_coeff0;
elseif 2<=n<4
    Table_coeff = Table_coeff1;
elseif 4<=n<8
    Table_coeff = Table_coeff2;
elseif 8<=n
    Table_coeff = Table_coeff3;
end

% Assign coeff_token and append it to output bits
% Here coeff_token is cell array so needs to be coverted to char

%==========================================================================
% JEC: XORing
%==========================================================================
% AES_count=AES_count+1;
% AESBit=AES_Bits(AES_count)*1;
% if (i_total==0)
%     Enc_i_total=bitxor(i_total,AESBit);
% else
%     Enc_i_total=bitxor(i_total-1,AESBit);
% end
%==========================================================================

coeff_token = Table_coeff(i_total + 1,i_trailing + 1);
bits = [bits char(coeff_token)];

% If the total coefficients == 0 exit from this function
if i_total==0
    return;
end

% Append sign of trailing ones to bits
if i_trailing>0
    bits = [bits sign];
end

%% Encode the levels of remaining non-zero coefficients

% find the suffix length
if (i_total>10)&(i_trailing<3)
    i_sufx_len = 1;
else
    i_sufx_len = 0;
end

% loop
for i=(i_trailing + 1):i_total
    
    if level(i)<0
        i_level_code = -2*level(i) - 1;
    else
        i_level_code = 2*level(i) - 2;
    end
    
    if (i == i_trailing + 1)&(i_trailing<3)
        i_level_code = i_level_code - 2;
    end
    
    if bitshift(i_level_code,-i_sufx_len)<14
        level_prfx = bitshift(i_level_code,-i_sufx_len);
        while(level_prfx>0)
            bits = [bits '0'];
            level_prfx = level_prfx - 1;
        end
        bits = [bits '1'];
        
        if i_sufx_len>0
            level_sufx = dec2bin(i_level_code,i_sufx_len);
            x = length(level_sufx);
            if x>i_sufx_len
                level_sufx = level_sufx(x-i_sufx_len+1:x);
            end
            bits = [bits level_sufx];
        end
    elseif (i_sufx_len==0)&(i_level_code<30)
        level_prfx = 14;
        while(level_prfx>0)
            bits = [bits '0'];
            level_prfx = level_prfx - 1;
        end
        bits = [bits '1'];
        
        level_sufx = dec2bin(i_level_code-14,4);
        x = length(level_sufx);
        if x>4
            level_sufx = level_sufx(x-4+1:x);
        end
        bits = [bits level_sufx];
        
    elseif (i_sufx_len>0)&(bitshift(i_level_code,-i_sufx_len)==14)
        level_prfx = 14;
        while(level_prfx>0)
            bits = [bits '0'];
            level_prfx = level_prfx - 1;
        end
        bits = [bits '1'];
        
        level_sufx = dec2bin(i_level_code,i_sufx_len);
        x = length(level_sufx);
        if x>i_sufx_len
            level_sufx = level_sufx(x-i_sufx_len+1:x);
        end
        bits = [bits level_sufx];
    else
        level_prfx = 15;
        while(level_prfx>0)
            bits = [bits '0'];
            level_prfx = level_prfx - 1;
        end
        bits = [bits '1'];
        
        i_level_code = i_level_code - bitshift(15,i_sufx_len);
        
        if i_sufx_len==0
            i_level_code = i_level_code - 15;
        end
        
        if (i_level_code>=bitshift(1,12))|(i_level_code<0)
            disp('Overflow occured');
        end
        
        level_sufx = dec2bin(i_level_code,12);
        x = length(level_sufx);
        if x>12
            level_sufx = level_sufx(x-12+1:x);
        end
        bits = [bits level_sufx];
    end
    
    if i_sufx_len==0
        i_sufx_len = i_sufx_len + 1;
    end
    if ((abs(level(i)))>bitshift(3,i_sufx_len - 1))&(i_sufx_len<6)
        i_sufx_len = i_sufx_len + 1;
    end
    
end

%% Encode Total zeros

% Here Rows(1-16) are Total coefficient and colums(0-15) are total zeros
% Rearranged from the standard for simplicity
% Table_zeros is located in table.mat or cavlc_tables.m file

if i_total<16
    total_zeros = Table_zeros(i_total,i_total_zero + 1);
    bits = [bits char(total_zeros)];
end

%% Encode each run of zeros
% Rows are the run before, and columns are zeros left
% Table_run is located in table.mat or cavlc_tables.m file

i_zero_left = i_total_zero;
if i_zero_left>=1
    for i=1:i_total
        if (i_zero_left>0)&(i==i_total)
            break;
        end
        if i_zero_left>=1
            i_zl = min(i_zero_left,7);
            run_before = Table_run(1 + run(i),i_zl);
            bits = [bits char(run_before)];
            i_zero_left = i_zero_left - run(i);
        end
    end
end

