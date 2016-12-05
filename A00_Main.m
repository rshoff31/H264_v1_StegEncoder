%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Written by:  Rich Shoff
% Source:      https://www.mathworks.com/matlabcentral/fileexchange/45151-h-264-baseline-codec--reading-and-writing-videos-
% Date:        12/3/2016
% Course:      EN.525.759: Image Compression, Packet Video, & Video
%                  Processing.
% Assigment:   Final Project: 'Steganography and Its Effects on H.264
%                  Compressed Video.'
% USAGE:       - Modify the variables below under 'Variables to Modify!'
%              - Run with Matlab.

%% Clean all and Close all previuse compution
clc
clear all;
close all;
system_dependent('DirChangeHandleWarn', 'Never');
addpath(genpath('.'));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Variables to Modify!!!
% Encoding Video Paramiters
Frame_start = 1;       % I frame
Frame_end = 6;         % Following P frames
% Operation Type:
operation = 'decode';  % Choices are encode or decode.
% Video (uncomment only 1):
%videoName = '.\Videos\CoastGuard_Raw_jpg.avi';
%videoName = '.\Videos\formancif.mov';
videoName = '.\Videos\Hall_Raw_jpg.avi';
% Secret Message:
global secretMessFile; secretMessFile = 'davinci_lc_over100mod';  % SECRET MESSAGE!!!
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Initialize Steganography:
% Program follows JSTEG model.
% Program  sets 0 as -2or2 && 1 as -1or3!
global steg; steg = 1;  % Set this to toggle steganography!!! 1=enabled, else disabled!
global debugInd; debugInd = 1;  % Set to 1 if you want to debug!
% #################################################################
if steg == 1
    global secretMess secretMessBin1r secretMessLen stegCapCount;
    global zeroCount; zeroCount = 0;  % Tracks # of -2 and 2 coefficients (indicating 0). steg and non-steg.
    global secretMessInd; secretMessInd = 1;
    secretMess = fileread(['.\SecretMessages\' secretMessFile '.txt']);  % Read secret message from text!
    secretMessLen = length(secretMess) * 7;
    secretMessBin = dec2bin(double(secretMess));
    secretMessBin1r = reshape(secretMessBin', 1, numel(secretMessBin));
end
% #################################################################


if strcmp('encode',operation);
    
    NumberOfFrames=Frame_start+Frame_end;
    
    %% Read video and save it in 00OriginalVideo.mat
    A01Read(videoName,NumberOfFrames); % You have to do it once, then comment it
    load 00OriginalVideo.mat;
    %% Option: playing & resizing 
    A01Play(mov,Obj);
    A01Resize(mov,Frame_start,Frame_end,128,128)
    %% Encoding Inputs/paramiters
    load 01ResizedFrames.mat;   % 1. Input video sequance (VideoSeq_Input)
    QP = 27;                    % 2. Quality Praramter, QP values

    %% Encoding Setup
    global h w              % Hight and Width
    block_size = 16;        % Macroblock size for P frames
    ext = 0;                % switch for extended ME
    [h,w,N] = size(VideoSeq_Input); % Extract input video sequance diminsion

    %% Encoding Outputs intialization
    Frames_PSNR = zeros(N,1);       % Initialize PSNR for each frame
    Frames_BR =  zeros(N,1);        % Initialize BitRate for each frame
    VideoSeq_Rec = zeros(h,w,N);    % Initialize Reconstructed video sequance
    bitstream = '';                 % Initialize output bitstream
    
    %% Encoding Processes
    %% 1. Saving the header: Encoding Inputs/paramiters (1x byte each)
    [bits] = header(h,w,QP,Frame_start,Frame_end);
    bitstream = bits;
    %% 2. Encode I-Frame
    disp(['Encoding I-Frame: ',num2str(Frame_start)]);      % Dispaly
    bitstream = [bitstream '1111'];     % Appending I-Frame header ('1111')
    Seq(:,:,1) = double(VideoSeq_Input(:,:,1)); % Extracting I-Frame
    [Seq_r(:,:,1),bits] = encode_i_frame(Seq(:,:,1),QP);    % Encoding
    %imshow(Seq_r)  % DEBUG
    Frames_Rec(:,:,1)=Seq(:,:,1);   % Storing Rec. Frame
    bitstream = [bitstream bits];   % Appending I-Frame bitstream ('11100...')
    X(:,:,1) = Seq_r(:,:,1);        % X[1]: Refrence Frame (for P-Frames)
    %% 3. Encoding P-Frams

    % #################################################################
    % STEGANOGRAPHY!!!
    if steg == 1
        stegCapCount = 0;  % This tracks # of bits that can be stored in video!
    end
    % #################################################################

    for K = 2:Frame_end
        k = K-1;
        disp(['Encoding P Frame: ',num2str(K)]);            % Dispaly
        bitstream = [bitstream '0000'];     % Appending P-Frame header ('0000')
        Seq(:,:,2) = double(VideoSeq_Input(:,:,K)); % Extracting P-Frame
        X(:,:,2) = Seq(:,:,2); % X[2]: Rnter coded Frame (P-Frames)
        [Seq_r(:,:,2),bits] = encode_p_frame(X,QP,ext,block_size); % Encoding
        Frames_Rec(:,:,K)=Seq_r(:,:,2); % Storing Rec. Frame
        bitstream = [bitstream bits]; % Appending P-Frame bitstream ('1100...')
        X(:,:,1) = Seq_r(:,:,2); % X[1]: Refrence Frame (for the next P-Frames)
    end

    % #################################################################
    % STEGANOGRAPHY!!!
    if steg == 1
        fprintf('############ Steganography Usage Stats ############\n\n' )

        fprintf('Video used:                     %s\n', videoName)
        [x, y, z] =  size(Seq);
        fprintf('Used resolution:                %i x %i\n', x, y)
        fprintf('Number of P frames:             %i\n\n', Frame_end-1)

        fprintf('Secret Message Used:            %s\n\n', secretMessFile)

        fprintf('bit capacity for steg:          %i\n', stegCapCount )
        fprintf('# bits in secret message:       %i\n', secretMessLen )
        fprintf('# bits used for steg:           %i\n\n', secretMessInd-1 )

        fprintf('character capacitity for steg:  %i\n', floor(stegCapCount/7) )
        fprintf('# characters in secret message: %i\n', floor(secretMessLen/7) )
        fprintf('# characters used for steg:     %i\n\n', floor(secretMessInd/7) )

        fprintf('Percent capacity used:          %i\n', floor((secretMessInd/stegCapCount)*100) )
        fprintf('Percent 0-bit Coefficients:     %f\n\n', (zeroCount/stegCapCount)*100 )

        fprintf('Video Size (KB):                %f\n\n', ((size(bitstream,2))/8) / (1024) );

        fprintf('###################################################\n' )
    end
    % #################################################################

    %% Storage
    save('02BitStream.mat','bitstream');                % 1. Bitstream
    save('02ReconstructedVideo.mat','Frames_Rec');      % 2. Reconstracted Frames
    save('02LastVideoName.mat','videoName');            % 3. Video encoded...used for decoding.
    
    % END ENCODING...

else
    
    %% H.264 Decoding
    clc
    clear all;
    close all;
    load 01ResizedFrames.mat;
    
    %% Decoding Bitsteam
    load 02BitStream.mat
    load 02LastVideoName.mat
    
    %% load 02ReconstructedVideo.mat
    global h w QP block_size
    idx = 1;
    block_size = 16;

    steg = 1;
    % #################################################################
    % STEGANOGRAPHY!!!
    if steg == 1
        global debugInd; debugInd = 1;  % Set to 1 if you want to debug!
    end
    % #################################################################

    %---------------------------------------------------------
    %% Decode header
    [h,w,QP,Frame_start,Frame_end,m] = dec_header(bitstream);
    idx = idx + m - 1;
    N = 1 + (Frame_end - Frame_start);
    %% Decode I-Frame
    if (bitstream(idx:idx+3)=='1111')
       disp(['Decoding I Frame: ',num2str(Frame_start)])
       idx = idx + 4;
       [Ceq_r(:,:,1),idx]=decode_i_frame(idx,bitstream);
       Frames_Dec(:,:,1)=Ceq_r(:,:,1);
    end
    %% Decode the following P-Frames
    for k = 2:N       
        if (bitstream(idx:idx+3)=='0000')
            disp(['Decoding P Frame: ', num2str(k)])
            idx = idx + 4;
            [Ceq_r(:,:,k),idx]= decode_p_frame(idx,bitstream,Ceq_r(:,:,k-1));
            Frames_Dec(:,:,k)=Ceq_r(:,:,k);
        end  
    end
    % End the decoding
    %% Option: Save the decoded video and paly it
    save('03DecodedVideo.mat','Frames_Dec');
    A01SaveAndPlay('HaveDecoded.avi',Frames_Dec);

    %Frames_Dec_CGorig = Frames_Dec;  % TEMP.
    %save('03DecodedVideo_CGorig.mat','Frames_Dec_CGorig');
    %Frames_Dec_FMorig = Frames_Dec;  % TEMP.
    %save('03DecodedVideo_FMorig.mat','Frames_Dec_FMorig');
    %Frames_Dec_HWorig = Frames_Dec;  % TEMP.
    %save('03DecodedVideo_HWorig.mat','Frames_Dec_HWorig');

    % #################################################################
    % STEGANOGRAPHY!!!
    if steg == 1
        % Load 3 videos for analysis:
        load 01ResizedFrames.mat;

        % Load unmodified Quantized video that corresponds with the chosen one:
        if strcmp(videoName, '.\Videos\CoastGuard_Raw_jpg.avi')
            load 03DecodedVideo_CGorig.mat;
            Frames_Dec_Orig = Frames_Dec_CGorig;
        elseif strcmp(videoName, '.\Videos\formancif.mov')
            load 03DecodedVideo_FMorig.mat;
            Frames_Dec_Orig = Frames_Dec_FMorig;
        else strcmp(videoName, '.\Videos\Hall_Raw_jpg.avi')
            load 03DecodedVideo_HWorig.mat;
            Frames_Dec_Orig = Frames_Dec_HWorig;
        end
           
        mse = 0; psnrat = 0; plotInd = 0;

        % ##### CALCS #####:
        for k = 2:N
            plotInd = plotInd + 1;

            frameOrig = uint8(VideoSeq_Input(:,:,k));
            frameQNoSteg = uint8(Frames_Dec_Orig(:,:,k));
            frameQSteg = uint8(Frames_Dec(:,:,k));

            mse = mse + immse( frameOrig, frameQSteg );
            psnrat = psnrat + psnr( frameOrig, frameQSteg );

            %title('Quantized (steg. modified) VS. Quantized (original)')
            subplot(3,5,plotInd), imshow(frameOrig)  % Orig Image (Unquantized/modified).
            subplot(3,5,plotInd+5), imshow(frameQNoSteg)  % Orig Quant. Image.
            subplot(3,5,plotInd+10), imshow(frameQSteg)  % Stegged Image.

        end

        fprintf('############ Steganography Usage Stats ############\n\n' )

        fprintf('Cumulative MSE (2-6):           %i\n', mse )
        fprintf('Cumulative PSNR (2-6):          %i\n', psnrat )

        fprintf('###################################################\n' )
    end
    % #################################################################

end  % End encode/decode.