function [Seq_r,bits_frame] = encode_i_frame_Secure(Seq,Quant)

%---------------------------------------------------------
%% H.264 Encoder - 
% it works on monochrome image sequence
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This script/function was created by 
% Abdullah Al Muhit
% contact - almuhit@gmail.com
% website - https://sites.google.com/site/almuhit/
% Please use it at your own risk. Also, Please cite the following paper:
% A A Muhit, M R Pickering, M R Frater and J F Arnold, “Video Coding using Elastic Motion Model and Larger Blocks,” IEEE Trans. Circ. And Syst. for Video Technology, vol. 20, no. 5, pp. 661-672, 2010. [Impact factor – 3.18] [PDF]
% A A Muhit, M R Pickering, M R Frater and J F Arnold, “Video Coding using Geometry Partitioning and an Elastic Motion Model,” accepted for publication in Journal of Visual Communication and Image Representation. [Impact factor – 1.33] [PDF]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

global QP h w AES_count
global Table_coeff0 Table_coeff1 Table_coeff2 Table_coeff3
global Table_run Table_zeros
load table.mat

[h,w,u] = size(Seq); 
QP = Quant;
AES_count = 0;
[bits_frame1,Seq_r1,sae1] = intra_4(Seq);
% AES_count = 0;
% [bits_frame2,Seq_r2,sae2] = intra_16(Seq);

% disp(['intra_4:',num2str(sae1)])
% disp(['intra_16:',num2str(sae2)])

% if sae1<sae2
    bits_frame = bits_frame1;
    Seq_r = Seq_r1;
    
% else
%     bits_frame = bits_frame2;
%     Seq_r = Seq_r2;  
% end

%--------------------------------------------------------------
function [bits_frame,Seq_r,total_sae] = intra_4(Seq)

global h w

bits_frame = '';
total_sae = 0;
mode_prev = 0;;
mode = 0;
mb_type = 1;            % 0 denotes 16x16, 1 denotes 4x4

for m = 1:16:h
    for n = 1:16:w
        bits_frame = [bits_frame dec2bin(mb_type)];    % mb header
        for i = m:4:m+15
            for j = n:4:n+15
                if (i==1)&(j==1)    % No prediciton
                    mode = 9;       % Special mode to describe no prediction in 4x4 
                    bits = enc_golomb(mode - mode_prev, 1);
                    mode_prev = mode;
                    bits_frame = [bits_frame bits];
                    [Seq_r(i:i+3,j:j+3,1),bits] = code_block_Secure(Seq(i:i+3,j:j+3,1));
                    bits_frame = [bits_frame bits];
                elseif (i==1)       % Horz prediction
                    mode = 1;       
                    bits = enc_golomb(mode - mode_prev, 1);
                    mode_prev = mode;
                    bits_frame = [bits_frame bits];

                    [icp,pred,sae] = pred_horz_4(Seq,Seq_r,i,j);
                    [icp_r,bits] = code_block_Secure(icp);
                    bits_frame = [bits_frame bits];
                    Seq_r(i:i+3,j:j+3,1)= icp_r + pred;
                    total_sae = total_sae + sae;
                elseif (j==1)       % Vert prediction
                    mode = 0;       
                    bits = enc_golomb(mode - mode_prev, 1);
                    mode_prev = mode;
                    bits_frame = [bits_frame bits];

                    [icp,pred,sae] = pred_vert_4(Seq,Seq_r,i,j);
                    [icp_r,bits] = code_block_Secure(icp);
                    bits_frame = [bits_frame bits];
                    Seq_r(i:i+3,j:j+3,1)= icp_r + pred;
                    total_sae = total_sae + sae;
                else                % Try all different prediction
                    [icp,pred,sae,mode] = mode_select_4(Seq,Seq_r,i,j); 
                    bits = enc_golomb(mode - mode_prev, 1);
                    mode_prev = mode;
                    bits_frame = [bits_frame bits];

                    [icp_r,bits] = code_block_Secure(icp);
                    bits_frame = [bits_frame bits];
                    Seq_r(i:i+3,j:j+3,1)= icp_r + pred;
                    total_sae = total_sae + sae;
                end

            end
        end
    end
end

%% 4x4 Horizontal prediciton

function [icp,pred,sae] = pred_horz_4(Seq,Seq_r,i,j)

pred = Seq_r(i:i+3,j-1)*ones(1,4);
icp = Seq(i:i+3,j:j+3,1) - pred;
sae = sum(sum(abs(icp)));

%-------------------------------------------------------
%% 4x4 Vertical Prediciton

function [icp,pred,sae] = pred_vert_4(Seq,Seq_r,i,j)

pred = ones(4,1)*Seq_r(i-1,j:j+3);
icp = Seq(i:i+3,j:j+3,1) - pred;
sae = sum(sum(abs(icp)));

%-------------------------------------------------------
%% 4x4 DC prediction

function [icp,pred,sae] = pred_dc_4(Seq,Seq_r,i,j)

pred = bitshift(sum(Seq_r(i-1,j:j+3))+ sum(Seq_r(i:i+3,j-1))+4,-3);
icp = Seq(i:i+3,j:j+3,1) - pred;
sae = sum(sum(abs(icp)));

%--------------------------------------------------------
function [icp,pred,sae] = pred_ddl_4(Seq,Seq_r,i,j)

Seq_r(i+4:i+7,j-1)=Seq_r(i+3,j-1);

for x = 0:3
    for y = 0:3
        if (x==3)&(y==3)
            pred(x+1,y+1) = bitshift(Seq_r(i+6,j-1) + 3*Seq_r(i+7,j-1) + 2,-2);
        else
            pred(x+1,y+1) = bitshift(Seq_r(i+x+y,j-1) + 2*Seq_r(i+x+y+1,j-1) + Seq_r(i+x+y+2,j-1) + 2,-2);
        end
    end
end

icp = Seq(i:i+3,j:j+3,1) - pred;
sae = sum(sum(abs(icp)));

%--------------------------------------------------------
function [icp,pred,sae] = pred_ddr_4(Seq,Seq_r,i,j)

for x = 0:3
    for y = 0:3
        if (x>y)
            pred(x+1,y+1) = bitshift(Seq_r(i+x-y-2,j-1) + 2*Seq_r(i+x-y-2,j-1) + Seq_r(i+x-y,j-1) + 2,-2);
        elseif (x<y)
            pred(x+1,y+1) = bitshift(Seq_r(i-1,j+y-x-2) + 2*Seq_r(i-1,j+y-x-1) + Seq_r(i-1,j+y-x) + 2,-2);
        else
            pred(x+1,y+1) = bitshift(Seq_r(i,j-1) + 2*Seq_r(i-1,j-1) + Seq_r(i-1,j) + 2,-2);
        end
    end
end

icp = Seq(i:i+3,j:j+3,1) - pred;
sae = sum(sum(abs(icp)));


%--------------------------------------------------------
function [icp,pred,sae] = pred_vr_4(Seq,Seq_r,i,j)

for x = 0:3
    for y = 0:3
        z = 2*x-y;
        w = bitshift(y,-1);
        if rem(z,2)==0
            pred(x+1,y+1)= bitshift(Seq_r(i+x-w-1,j-1) + Seq_r(i+x-w,j-1) + 1,-1);
        elseif rem(z,2)==1
            pred(x+1,y+1)= bitshift(Seq_r(i+x-w-2,j-1) + 2*Seq_r(i+x-w-1,j-1) + Seq_r(i+x-w,j-1) + 2,-2);
        elseif z==-1
            pred(x+1,y+1)= bitshift(Seq_r(i-1,j)+ 2*Seq_r(i-1,j-1) + Seq_r(i,j-1) + 2,-2);
        else
            pred(x+1,y+1) = bitshift(Seq_r(i-1,j+y-1)+ 2*Seq_r(i-1,j+y-2) + Seq_r(i-1,j+y-3) + 2,-2);
        end
    end
end

icp = Seq(i:i+3,j:j+3,1) - pred;
sae = sum(sum(abs(icp)));

%--------------------------------------------------------
function [icp,pred,sae] = pred_hd_4(Seq,Seq_r,i,j)

for x = 0:3
    for y = 0:3
        z = 2*y-x;
        w = bitshift(x,-1);
        if rem(z,2)==0
            pred(x+1,y+1)= bitshift(Seq_r(i-1,j+y-w-1) + Seq_r(i-1,j+y-w) + 1,-1);
        elseif rem(z,2)==1
            pred(x+1,y+1)= bitshift(Seq_r(i-1,j+y-w-2) + 2*Seq_r(i-1,j+y-w-1) + Seq_r(i-1,j+y-w) + 2,-2);
        elseif z==-1
            pred(x+1,y+1)= bitshift(Seq_r(i-1,j)+ 2*Seq_r(i-1,j-1) + Seq_r(i,j-1) + 2,-2);
        else
            pred(x+1,y+1) = bitshift(Seq_r(i+x-1,j-1)+ 2*Seq_r(i+x-2,j-1) + Seq_r(i+x-3,j-1) + 2,-2);
        end
    end
end

icp = Seq(i:i+3,j:j+3,1) - pred;
sae = sum(sum(abs(icp)));

%--------------------------------------------------------
function [icp,pred,sae] = pred_vl_4(Seq,Seq_r,i,j)

Seq_r(i+4:i+7,j-1)=Seq_r(i+3,j-1);

for x = 0:3
    for y = 0:3
        w = bitshift(y,-1);
        if rem(y,2)==0
            pred(x+1,y+1) = bitshift(Seq_r(i+x+w,j-1) + Seq_r(i+x+w+1,j-1) + 1,-1);
        else
            pred(x+1,y+1) = bitshift(Seq_r(i+x+w,j-1) + 2*Seq_r(i+x+w+1,j-1) + Seq_r(i+x+w+2,j-1) + 2,-2);
        end
    end
end

icp = Seq(i:i+3,j:j+3,1) - pred;
sae = sum(sum(abs(icp)));

%--------------------------------------------------------
function [icp,pred,sae] = pred_hu_4(Seq,Seq_r,i,j)

for x = 0:3
    for y = 0:3
        z = 2*y+x;
        w = bitshift(x,-1);
        if (z==0)|(z==2)|(z==4)
            pred(x+1,y+1)= bitshift(Seq_r(i-1,j+y+w) + Seq_r(i-1,j+y+w+1) + 1,-1);
        elseif (z==1)|(z==3)
            pred(x+1,y+1)= bitshift(Seq_r(i-1,j+y+w) + 2*Seq_r(i-1,j+y+w+1) + Seq_r(i-1,j+y+w+2) + 2,-2);
        elseif z==5
            pred(x+1,y+1)= bitshift(Seq_r(i-1,j+2)+ 3*Seq_r(i-1,j+3) + 2,-2);
        else
            pred(x+1,y+1) = Seq_r(i-1,j+3);
        end
    end
end

icp = Seq(i:i+3,j:j+3,1) - pred;
sae = sum(sum(abs(icp)));

%---------------------------------------------------------
%% Mode selection for 4x4 prediciton

function [icp,pred,sae,mode] = mode_select_4(Seq,Seq_r,i,j)

[icp1,pred1,sae1] = pred_vert_4(Seq,Seq_r,i,j);
[icp2,pred2,sae2] = pred_horz_4(Seq,Seq_r,i,j);
[icp3,pred3,sae3] = pred_dc_4(Seq,Seq_r,i,j);
[icp4,pred4,sae4] = pred_ddl_4(Seq,Seq_r,i,j);
[icp5,pred5,sae5] = pred_ddr_4(Seq,Seq_r,i,j);
[icp6,pred6,sae6] = pred_vr_4(Seq,Seq_r,i,j);
[icp7,pred7,sae7] = pred_hd_4(Seq,Seq_r,i,j);
[icp8,pred8,sae8] = pred_vl_4(Seq,Seq_r,i,j);
[icp9,pred9,sae9] = pred_hu_4(Seq,Seq_r,i,j);

[val,idx]=min([sae1 sae2 sae3 sae4 sae5 sae6 sae7 sae8 sae9]);

switch idx
    case 1
        sae = sae1;
        icp = icp1;
        pred = pred1; 
        mode = 0;
    case 2
        sae = sae2;
        icp = icp2;
        pred = pred2;
        mode = 1;
    case 3
        sae = sae3;
        icp = icp3;
        pred = pred3;
        mode = 2;
    case 4
        sae = sae4;
        icp = icp4;
        pred = pred4; 
        mode = 3;
    case 5
        sae = sae5;
        icp = icp5;
        pred = pred5; 
        mode = 4;
    case 6
        sae = sae6;
        icp = icp6;
        pred = pred6; 
        mode = 5;
    case 7
        sae = sae7;
        icp = icp7;
        pred = pred7; 
        mode = 6;
    case 8
        sae = sae8;
        icp = icp8;
        pred = pred8; 
        mode = 7;
    case 9
        sae = sae9;
        icp = icp9;
        pred = pred9; 
        mode = 8;
end


%--------------------------------------------------------------
function [bits_frame,Seq_r,total_sae] = intra_16(Seq)

global h w

bits_frame = '';
total_sae = 0;
mode_prev = 0;
mode = 0;
mb_type = 0;            % 0 denotes 16x16, 1 denotes 4x4

for i = 1:16:h
    for j = 1:16:w
        if (i==1)&(j==1)    % No prediciton
            mode = 4;       % Special mode to describe no prediction
            [bits,mode_prev]= mb_header(mb_type,mode,mode_prev);
            bits_frame = [bits_frame bits];
            
            [Seq_r(i:i+15,j:j+15,1),bits] = code_block_Secure(Seq(i:i+15,j:j+15,1));
            bits_frame = [bits_frame bits];
        elseif (i==1)       % Horz prediction
            mode = 1;       
            [bits,mode_prev]= mb_header(mb_type,mode,mode_prev);
            bits_frame = [bits_frame bits];
            
            [icp,pred,sae] = pred_horz_16(Seq,Seq_r,16,i,j);
            [icp_r,bits] = code_block_Secure(icp);
            bits_frame = [bits_frame bits];
            Seq_r(i:i+15,j:j+15,1)= icp_r + pred;
            total_sae = total_sae + sae;
        elseif (j==1)       % Vert prediction
            mode = 0;       
            [bits,mode_prev]= mb_header(mb_type,mode,mode_prev);
            bits_frame = [bits_frame bits];
            
            [icp,pred,sae] = pred_vert_16(Seq,Seq_r,16,i,j);
            [icp_r,bits] = code_block_Secure(icp);
            bits_frame = [bits_frame bits];
            Seq_r(i:i+15,j:j+15,1)= icp_r + pred;
            total_sae = total_sae + sae;
        else                % Try all different prediction
            [icp,pred,sae,mode] = mode_select_16(Seq,Seq_r,16,i,j);
            [bits,mode_prev]= mb_header(mb_type,mode,mode_prev);
            bits_frame = [bits_frame bits];
            
            [icp_r,bits] = code_block_Secure(icp);
            bits_frame = [bits_frame bits];
            Seq_r(i:i+15,j:j+15,1)= icp_r + pred;
            total_sae = total_sae + sae;
        end
        
    end
end



%--------------------------------------------------------
%% Transform, Quantization, Entropy coding
% transform = Integer transform
% Quantization = h.264 
% VLC = CAVLC (H.264)

function [err_r,bits_mb] = code_block(err)

global QP
[n,m] = size(err);

bits_mb = '';

for i = 1:4:n
    for j = 1:4:m
        c(i:i+3,j:j+3) = integer_transform(err(i:i+3,j:j+3));
        cq(i:i+3,j:j+3) = quantization(c(i:i+3,j:j+3),QP);
        %==================================================================

        %==================================================================
        [bits_b] = enc_cavlc(cq(i:i+3,j:j+3), 0, 0);       
        bits_mb = [bits_mb bits_b];
        %==================================================================
        
        %==================================================================
        Wi = inv_quantization(cq(i:i+3,j:j+3),QP);
        Y = inv_integer_transform(Wi);        
        err_r(i:i+3,j:j+3) = round(Y/64);        
    end
end

function [err_r,bits_mb] = code_block_Secure(err)
global QP
[n,m] = size(err);

bits_mb = '';

for i = 1:4:n
    for j = 1:4:m
        c(i:i+3,j:j+3) = integer_transform(err(i:i+3,j:j+3));
        cq(i:i+3,j:j+3) = quantization(c(i:i+3,j:j+3),QP);
        %==================================================================

        %==================================================================
        [bits_b] = enc_cavlc_Secure(cq(i:i+3,j:j+3), 0, 0);       
        bits_mb = [bits_mb bits_b];
        %==================================================================
        
        %==================================================================
        Wi = inv_quantization(cq(i:i+3,j:j+3),QP);
        Y = inv_integer_transform(Wi);        
        err_r(i:i+3,j:j+3) = round(Y/64);        
    end
end

%-------------------------------------------------------
%% 16x16 Horizontal prediciton

function [icp,pred,sae] = pred_horz_16(Seq,Seq_r,bs,i,j)

pred = Seq_r(i:i+15,j-1)*ones(1,bs);
icp = Seq(i:i+15,j:j+15,1) - pred;
sae = sum(sum(abs(icp)));

%-------------------------------------------------------
%% 16x16 Vertical Prediciton

function [icp,pred,sae] = pred_vert_16(Seq,Seq_r,bs,i,j)

pred = ones(bs,1)*Seq_r(i-1,j:j+15);
icp = Seq(i:i+15,j:j+15,1) - pred;
sae = sum(sum(abs(icp)));

%-------------------------------------------------------
%% 16x16 DC prediction

function [icp,pred,sae] = pred_dc_16(Seq,Seq_r,bs,i,j)

pred = bitshift(sum(Seq_r(i-1,j:j+15))+ sum(Seq_r(i:i+15,j-1))+16,-5);
icp = Seq(i:i+15,j:j+15,1) - pred;
sae = sum(sum(abs(icp)));

%------------------------------------------------------
%% 16x16 Plane prediction

function [icp,pred,sae] = pred_plane_16(Seq,Seq_r,bs,i,j)

x = 0:7;
H = sum((x+1)*(Seq_r(i+x+8,j-1)-Seq_r(i+6-x,j-1)));
y = 0:7;
V = sum((y+1)*(Seq_r(i-1,j+8+y)'-Seq_r(i-1,j+6-y)'));

a = 16*(Seq_r(i-1,j+15) + Seq_r(i+15,j-1));
b = bitshift(5*H + 32,-6,'int64');
c = bitshift(5*V + 32,-6,'int64');

% pred = clipy() << refer to the standard
for m = 1:16
    for n = 1:16
        d = bitshift(a + b*(m-8)+ c*(n-8) + 16, -5,'int64');
        if d <0
            pred(m,n) = 0;
        elseif d>255
            pred(m,n) = 255;
        else
            pred(m,n) = d;
        end
    end
end

icp = Seq(i:i+15,j:j+15,1) - pred;
sae = sum(sum(abs(icp)));

%---------------------------------------------------------
%% Mode selection for 16x16 prediciton

function [icp,pred,sae,mode] = mode_select_16(Seq,Seq_r,bs,i,j)

[icp1,pred1,sae1] = pred_vert_16(Seq,Seq_r,bs,i,j);
[icp2,pred2,sae2] = pred_horz_16(Seq,Seq_r,bs,i,j);
[icp3,pred3,sae3] = pred_dc_16(Seq,Seq_r,bs,i,j);
[icp4,pred4,sae4] = pred_plane_16(Seq,Seq_r,bs,i,j);

[val,idx]=min([sae1 sae2 sae3 sae4]);

switch idx
    case 1
        sae = sae1;
        icp = icp1;
        pred = pred1; 
        mode = 0;
    case 2
        sae = sae2;
        icp = icp2;
        pred = pred2;
        mode = 1;
    case 3
        sae = sae3;
        icp = icp3;
        pred = pred3;
        mode = 2;
    case 4
        sae = sae4;
        icp = icp4;
        pred = pred4; 
        mode = 3;
end

