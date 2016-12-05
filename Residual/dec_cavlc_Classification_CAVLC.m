function [data,i] = dec_cavlc_Classification_CAVLC(bits,nL,nU)

%% CAVLC Decoder
% By A. A. Muhit
% It takes bitstream and decodes 4x4 block of data


% Load the table containing all the tables
% load table.mat;
global Table_coeff0 Table_coeff1 Table_coeff2 Table_coeff3
global Table_run Table_zeros

global QP
global AES_count AES_Bits


% test data
% bits = '00000000001000101000000011001011111001011101010000000010000111111111111111001110';
% nL = 0;
% nU = 0;

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
% Choose proper Table_coeff based on n value
if 0<=n<2
    Table_coeff = Table_coeff0;
elseif 2<=n<4
    Table_coeff = Table_coeff1;
elseif 4<=n<8
    Table_coeff = Table_coeff2;
elseif 8<=n
    Table_coeff = Table_coeff3;
end

i = 1;
coeff_token = '';

% Find total coefficients and trailing ones
while (i<=length(bits))
    coeff_token = [coeff_token bits(i)];
    x = strcmp(Table_coeff,coeff_token);
    [r,c]=find(x==1);
    i = i + 1;
    if (r>0)&(c>0)
        break;
    end
end

% Find total coefficients and trailing ones
i_total = r - 1;
i_trailing = c - 1;

%==========================================================================
% JEC: XORing
%==========================================================================
% AES_count=AES_count+1;
% AESBit=AES_Bits(AES_count)*1;
% Enc_i_total=i_total;
% if (Enc_i_total==0)
%     i_total=bitxor(Enc_i_total,AESBit);
% else
%     i_total=bitxor(Enc_i_total+1,AESBit)+1;
% end

%==========================================================================




% if no coefficients return 4x4 empty blocks of data
if i_total==0
    data = zeros(4,4);
    return;
end

k = 1;
m = i_trailing;

while m>0
    if bits(i)=='0'
        level(k)=1;
    elseif bits(i)=='1'
        level(k)=-1;
    end
    k = k + 1;
    m = m - 1;
    i = i + 1;
end

%% Decode the non-zero coefficient/level values

if (i_total>10)&(i_trailing<3)
    i_sufx_len = 1;
else
    i_sufx_len = 0;
end

while k<=i_total
    % Decode level prefix
    [level_prfx,i]= dec_prfx(bits,i);
    
    % Decode level suffix
    level_sufx_size = 0;
    
    if (i_sufx_len>0)||(level_prfx>=14)
        if (level_prfx==14)&(i_sufx_len==0)
            level_sufx_size = 4;
        elseif level_prfx>=15
            level_sufx_size = level_prfx - 3;
        else
            level_sufx_size = i_sufx_len;
        end
    end
    
    if level_sufx_size==0
        level_sufx = 0;
    else
        sufx = bits(i : i + level_sufx_size -1);
        level_sufx = bin2dec(sufx);
        i = i + level_sufx_size;
    end
    
    i_level_code = bitshift(min(15,level_prfx),i_sufx_len,'int64') + level_sufx;
    
    if (level_prfx>=15)&(i_sufx_len==0)
        i_level_code = i_level_code + 15;
    end
    if level_prfx>=16
        i_level_code = i_level_code + (bitshift(1,level_prfx - 3,'int64') - 4096);
    end
    
    if (k == i_trailing + 1)&(i_trailing<3)
        i_level_code = i_level_code + 2;
    end
    
    if rem(i_level_code,2)==0 % i_level_code is even
        level(k) = bitshift(i_level_code + 2,-1,'int64');
    else % odd number
        level(k) = bitshift(-i_level_code - 1, -1,'int64');
    end
    
    if i_sufx_len==0
        i_sufx_len = 1;
    end
    
    if ((abs(level(k)))>bitshift(3,i_sufx_len - 1,'int64'))&(i_sufx_len<6)
        i_sufx_len = i_sufx_len + 1;
    end
    
    k = k + 1;
end

%% Decode total zeros

s='';
i_total_zero = 0;

if i_total==16
    i_zero_left = 0;
else
    while (i<=length(bits))
        s = [s bits(i)];
        x = strcmp(Table_zeros(i_total,:),s);
        r = find(x==1);
        i = i + 1;
        if r>0
            i_total_zero = r-1;
            break;
        end
    end
end


%% Decode run information

i_zero_left = i_total_zero;
j=1;
ss = '';
run = zeros(1,length(level));

while i_zero_left>0
    while (j<i_total)&(i_zero_left>0)
        ss = [ss bits(i)];
        i_zl = min(i_zero_left,7);
        x = strcmp(Table_run(:,i_zl),ss);
        r = find(x==1);
        i = i + 1;
        if r>0
            run(j)=r-1;
            i_zero_left = i_zero_left - run(j);
            j = j + 1;
            ss = '';
        end
    end
    if i_zero_left>0
        run(j)=i_zero_left;
        i_zero_left = 0;
    end
end

%% Combine level and run information

k = i_total + i_total_zero;
l = zeros(1,16);

while k>0
    for j=1:length(level)
        l(k)=level(j);
        k = k - 1;
        k = k - run(j);
    end
end

%% Reorder the data into 4x4 block


scan = [1,1;1,2;2,1;3,1;2,2;1,3;1,4;2,3;3,2;4,1;4,2;3,3;2,4;3,4;4,3;4,4];

for k=16:-1:1
    m=scan(k,1);
    n=scan(k,2);
    data(m,n)=l(k); % l contains the reordered data
    
    
end