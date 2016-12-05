function [Sr,idx] = decode_p_frame(idx,bitstream,Si)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This script/function was created by 
% Abdullah Al Muhit
% contact - almuhit@gmail.com
% website - https://sites.google.com/site/almuhit/
% Please use it at your own risk. Also, Please cite the following paper:
% A A Muhit, M R Pickering, M R Frater and J F Arnold, “Video Coding using Elastic Motion Model and Larger Blocks,” IEEE Trans. Circ. And Syst. for Video Technology, vol. 20, no. 5, pp. 661-672, 2010. [Impact factor – 3.18] [PDF]
% A A Muhit, M R Pickering, M R Frater and J F Arnold, “Video Coding using Geometry Partitioning and an Elastic Motion Model,” accepted for publication in Journal of Visual Communication and Image Representation. [Impact factor – 1.33] [PDF]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Initialization
global QP block_size
global mvm S h w 

global Table_coeff0 Table_coeff1 Table_coeff2 Table_coeff3
global Table_run Table_zeros
load table.mat

% Array for motion vectors
mvm = inf.*ones(h+1,w+2,2);
% reconstructed frame
Sr = zeros(h,w);
S = Si;

for i = 1:block_size:h
    for j = 1:block_size:w

        % bs - macroblock size
        bs = block_size;
        % decode MB header
        if (bitstream(idx)=='0')
%             disp('Decoding 16x16 MB')
            idx = idx + 1;
        end

        % determine the 16x16 partition mode
        % ind_16 - partition index (0 = full, 1 = vert, 2 = horiz, 3 = quad)
        [mode,idx]= dec_golomb(idx,bitstream,0);
        
        switch mode
            case 0
                [Sr(i:i+bs-1,j:j+bs-1),idx,mvx,mvy]= dec_block(idx,bitstream,bs,bs,i,j,'m');
                
                mvm(i:i+bs-1,j:j+bs-1,1) = mvx.*ones(bs,bs);
                mvm(i:i+bs-1,j:j+bs-1,2) = mvy.*ones(bs,bs);
            case 1
                [Sr(i:i+bs-1,j:j+bs/2-1),idx,mvx0,mvy0]= dec_block(idx,bitstream,bs/2,bs,i,j,'a');
                [Sr(i:i+bs-1,j+bs/2:j+bs-1),idx,mvx1,mvy1]= dec_block(idx,bitstream,bs/2,bs,i,j+bs/2,'c');
                
                mvm(i:i+bs-1,j:j+bs-1,1) = [mvx0.*ones(bs,bs/2) mvx1.*ones(bs,bs/2)];
                mvm(i:i+bs-1,j:j+bs-1,2) = [mvy0.*ones(bs,bs/2) mvy1.*ones(bs,bs/2)];
            case 2
                [Sr(i:i+bs/2-1,j:j+bs-1),idx,mvx0,mvy0]= dec_block(idx,bitstream,bs,bs/2,i,j,'b');
                [Sr(i+bs/2:i+bs-1,j:j+bs-1),idx,mvx1,mvy1]= dec_block(idx,bitstream,bs,bs/2,i+bs/2,j,'a');
                
                mvm(i:i+bs-1,j:j+bs-1,1) = [mvx0.*ones(bs/2,bs); mvx1.*ones(bs/2,bs)];
                mvm(i:i+bs-1,j:j+bs-1,2) = [mvy0.*ones(bs/2,bs); mvy1.*ones(bs/2,bs)];
            case 3
                [Sr(i:i+bs/2-1,j:j+bs/2-1),idx]= dec_block_8(idx,bitstream,bs/2,i,j);
                [Sr(i:i+bs/2-1,j+bs/2:j+bs-1),idx]= dec_block_8(idx,bitstream,bs/2,i,j+bs/2);
                [Sr(i+bs/2:i+bs-1,j:j+bs/2-1),idx]= dec_block_8(idx,bitstream,bs/2,i+bs/2,j);
                [Sr(i+bs/2:i+bs-1,j+bs/2:j+bs-1),idx]= dec_block_8(idx,bitstream,bs/2,i+bs/2,j+bs/2);
                
        end
        
%         % display the decoded frame
%         figure(3)
%         image(Sr)
%         colormap(gray(256))
%         truesize([2*h 2*w])
%         drawnow
    end
end

%--------------------------------------------------------------------------
function [Cr,idx] = dec_block_8(idx,bitstream,bs,i,j)

global mvm
% determine the 8x8 partition mode
% partition index (0 = full, 1 = vert, 2 = horiz, 3 = quad)
[mode,idx]= dec_golomb(idx,bitstream,0);

Cr = zeros(bs,bs);
x = 1;
y = 1;

switch mode
    case 0
        [Cr(x:x+bs-1,y:y+bs-1),idx,mvx,mvy]= dec_block(idx,bitstream,bs,bs,i,j,'m');
        
        mvm(i:i+bs-1,j:j+bs-1,1) = mvx.*ones(bs,bs);
        mvm(i:i+bs-1,j:j+bs-1,2) = mvy.*ones(bs,bs);
    case 1
        [Cr(x:x+bs-1,y:y+bs/2-1),idx,mvx0,mvy0]= dec_block(idx,bitstream,bs/2,bs,i,j,'a');
        [Cr(x:x+bs-1,y+bs/2:y+bs-1),idx,mvx1,mvy1]= dec_block(idx,bitstream,bs/2,bs,i,j+bs/2,'c');
        
        mvm(i:i+bs-1,j:j+bs-1,1) = [mvx0.*ones(bs,bs/2) mvx1.*ones(bs,bs/2)];
        mvm(i:i+bs-1,j:j+bs-1,2) = [mvy0.*ones(bs,bs/2) mvy1.*ones(bs,bs/2)];
    case 2
        [Cr(x:x+bs/2-1,y:y+bs-1),idx,mvx0,mvy0]= dec_block(idx,bitstream,bs,bs/2,i,j,'b');
        [Cr(x+bs/2:x+bs-1,y:y+bs-1),idx,mvx1,mvy1]= dec_block(idx,bitstream,bs,bs/2,i+bs/2,j,'a');
        
        mvm(i:i+bs-1,j:j+bs-1,1) = [mvx0.*ones(bs/2,bs); mvx1.*ones(bs/2,bs)];
        mvm(i:i+bs-1,j:j+bs-1,2) = [mvy0.*ones(bs/2,bs); mvy1.*ones(bs/2,bs)];
    case 3
        [Cr(x:x+bs/2-1,y:y+bs/2-1),idx,mvx0,mvy0]= dec_block(idx,bitstream,bs/2,bs/2,i,j,'m');
        [Cr(x:x+bs/2-1,y+bs/2:y+bs-1),idx,mvx1,mvy1]= dec_block(idx,bitstream,bs/2,bs/2,i,j+bs/2,'m');
        [Cr(x+bs/2:x+bs-1,y:y+bs/2-1),idx,mvx2,mvy2]= dec_block(idx,bitstream,bs/2,bs/2,i+bs/2,j,'m');
        [Cr(x+bs/2:x+bs-1,y+bs/2:y+bs-1),idx,mvx3,mvy3]= dec_block(idx,bitstream,bs/2,bs/2,i+bs/2,j+bs/2,'m');
        
        mvm(i:i+bs-1,j:j+bs-1,1) = [mvx0.*ones(bs/2,bs/2) mvx1.*ones(bs/2,bs/2);mvx2.*ones(bs/2,bs/2) mvx3.*ones(bs/2,bs/2)];
        mvm(i:i+bs-1,j:j+bs-1,2) = [mvy0.*ones(bs/2,bs/2) mvy1.*ones(bs/2,bs/2);mvy2.*ones(bs/2,bs/2) mvy3.*ones(bs/2,bs/2)];
end

%--------------------------------------------------------------------------
function [rec,k,mvx,mvy] = dec_block(k,bits,bsx,bsy,i,j,mvp_mode)

global S h w 

p = 8;
s = 2;
eb = 2;

% calculate the motion vector predictions mvxp and mvyp
switch mvp_mode
    case 'm'
        [mvxp,mvyp] = find_median_mv_pred(i,j,bsx);
    case 'a'
        [mvxp,mvyp] = find_mv_pred_a(i,j,bsx);
    case 'b'
        [mvxp,mvyp] = find_mv_pred_b(i,j,bsx);
    case 'c'
        [mvxp,mvyp] = find_mv_pred_c(i,j,bsx);
end

if bits(k)=='0'
%     disp('Translational MV')
    k = k + 1;
    % get the mvx and mvy difference from the bitstream
    [mvx_diff,k] = dec_golomb(k,bits,1);
    [mvy_diff,k] = dec_golomb(k,bits,1);

    % original mvx and mvy
    mvx = mvx_diff/4 + mvxp;
    mvy = mvy_diff/4 + mvyp;

    % finding mr is important for finding the prediction
    mr = find_mr(mvx,mvy);
    mr = mr + 1;
    y = 4*mvx + 4*mr +1;
    x = 4*mvy + 4*mr +1;

    % extract the search window from the previous reconstructed frame
    swr = extract_object(S,[j-mr i-mr bsx+2*mr bsy+2*mr],[1 1 w h]);
    [n,m] = size(swr);
    swr_sup = interp2([1:m],[1:n]',mirror_pad(swr),[1:.25:m],[1:.25:n]','bicubic');


    % extract motion compensated prediction block
    pred = swr_sup(x:4:x+4*bsy-1,y:4:y+4*bsx-1);

    % the residual data
    [mcp,k] = code_block(k,bits,bsx,bsy);
    % Reconstructed block
    rec = pred + mcp;
else
    disp('Extended MV')
    k = k + 1;
    % get the mvx or M(1) from the bitstream
    [M(1),k] = dec_golomb(k,bits,1);
    [M(2),k] = dec_golomb(k,bits,1);
    [M(3),k] = dec_golomb(k,bits,1);
    [M(4),k] = dec_golomb(k,bits,1);
    % get the mvy or M(5) from the bitstream
    [M(5),k] = dec_golomb(k,bits,1);
    [M(6),k] = dec_golomb(k,bits,1);
    [M(7),k] = dec_golomb(k,bits,1);
    [M(8),k] = dec_golomb(k,bits,1);
    
    M(1) = M(1)/p + mvxp;
    M(2) = M(2)/p;
    M(3) = M(3)/p;
    M(4) = M(4)/p;
    M(5) = M(5)/p + mvyp;
    M(6) = M(6)/p;
    M(7) = M(7)/p;
    M(8) = M(8)/p;
    
    mvx = M(1);
    mvy = M(5);
    
    % finding mr is important for finding the prediction
    mr = find_mr(mvx,mvy);
    mr = mr + 1;
    % extract the search window from the previous reconstructed frame
    prev_w = mirror_pad(extract_object(S,[j-mr i-mr bsx+2*mr bsy+2*mr],[1 1 w h]));
    I = prev_w(mr+1-eb:mr+bsy+eb,mr+1-eb:mr+bsx+eb);
    Iw = warp_2D_wav_cub_sa(I,M,s);
    
    % get prediction
    pred = Iw(eb+1:eb+bsy,eb+1:eb+bsx);
    
    % the residual data
    [mcp,k] = code_block(k,bits,bsx,bsy);
    % Reconstructed block
    rec = pred + mcp;
end



%--------------------------------------------------------------------------
function [mcp,k] = code_block(k,bits,bsx,bsy)

global QP debugInd

for i=1:4:bsy
    for j=1:4:bsx
        [Z1(i:i+3,j:j+3,1),m] = dec_cavlc(bits(k:length(bits)),0,0);
        Wi = inv_quantization(Z1(i:i+3,j:j+3,1),QP);
        Y = inv_integer_transform(Wi);
        X = round(Y/64);
        mcp(i:i+3,j:j+3,1) = X;
        k = k + m - 1;
    end
end

% #################################################################
% STEGANOGRAPHY!!!

if debugInd == 1
    % Show cq with & without steganography!
    fprintf('\n\n--------- Steg Block Decoded (Z1) ----------\n\n')
    Z1
    debugInd = 1000;
end

% #################################################################


%--------------------------------------------------------------------------
function [mvxp,mvyp] = find_median_mv_pred(i,j,bs);

global mvm

for k = 1:2

    mvpa = mvm(i+1,j,k);
    mvpb = mvm(i,j+1,k);
    mvpc = mvm(i,j+bs+1,k);
    mvpd = mvm(i,j,k);

    
    switch bin2dec([num2str(isinf(mvpa)) num2str(isinf(mvpb)) num2str(isinf(mvpc)) num2str(isinf(mvpd))])
        case 0
            mvp(k) = median([mvpa mvpb mvpc]);
        case 1
            mvp(k) = median([mvpa mvpb mvpc]);
        case 2
            mvp(k) = median([mvpa mvpb mvpd]);
        case 3
            mvp(k) = median([mvpa mvpb]);
        case 4
            mvp(k) = median([mvpa mvpc mvpd]);
        case 5
            mvp(k) = median([mvpa mvpc]);
        case 6
            mvp(k) = median([mvpa mvpd]);
        case 7
            mvp(k) = mvpa;
        case 8
            mvp(k) = median([mvpb mvpc mvpd]);
        case 9
            mvp(k) = median([mvpb mvpc]);
        case 10
            mvp(k) = median([mvpb mvpd]);
        case 11
            mvp(k) = mvpb;
        case 12
            mvp(k) = median([mvpc mvpd]);
        case 13
            mvp(k) = mvpc;
        case 14
            mvp(k) = mvpd;
        case 15
            mvp(k) = 0;
    end
end

mvxp = mvp(1);
mvyp = mvp(2);

% function [mvxp,mvyp] = find_median_mv_pred(i,j,bs);
% 
% global mvm
% 
% for k = 1:2
% 
%     mvpa = mvm(i+1,j,k);
%     mvpb = mvm(i,j+1,k);
%     mvpc = mvm(i,j+bs+1,k);
% 
%     
%     switch bin2dec([num2str(isinf(mvpa)) num2str(isinf(mvpb)) num2str(isinf(mvpc))])
%         case 0
%             mvp(k) = median([mvpa mvpb mvpc]);
%         case 1
%             mvp(k) = median([mvpa mvpb]);
%         case 2
%             mvp(k) = median([mvpa mvpc]);
%         case 3
%             mvp(k) = mvpa;
%         case 4
%             mvp(k) = median([mvpc mvpd]);
%         case 5
%             mvp(k) = mvpb;
%         case 6
%             mvp(k) = mvpc;
%         case 7
%             mvp(k) = 0;
%     end
% end
% 
% mvxp = mvp(1);
% mvyp = mvp(2);


%--------------------------------------------------------------------------
function [mvxp,mvyp] = find_mv_pred_a(i,j,bs);

global mvm

for k = 1:2

    mvpa = mvm(i+1,j,k);

    switch isinf(mvpa)
        case 0
            mvp(k) = mvpa;
        case 1
            mvp(k) = 0;
    end
end

mvxp = mvp(1);
mvyp = mvp(2);


%--------------------------------------------------------------------------
function [mvxp,mvyp] = find_mv_pred_b(i,j,bs);

global mvm

for k = 1:2

    mvpb = mvm(i,j+1,k);

    switch isinf(mvpb)
        case 0
            mvp(k) = mvpb;
        case 1
            mvp(k) = 0;
    end
end

mvxp = mvp(1);
mvyp = mvp(2);


%--------------------------------------------------------------------------
function [mvxp,mvyp] = find_mv_pred_c(i,j,bs);

global mvm

for k = 1:2

    mvpc = mvm(i,j+bs+1,k);

    switch isinf(mvpc)
        case 0
            mvp(k) = mvpc;
        case 1
            mvp(k) = 0;
    end
end

mvxp = mvp(1);
mvyp = mvp(2);


%------------------------------------------------------------------
function [mr] = find_mr(mvx,mvy)

mr = ceil(max(abs(mvx),abs(mvy)));