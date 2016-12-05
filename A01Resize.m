% Generates Frames
% clc
% clear all;
% close all;

%% Resizing 
function A01Resize(mov,Frame_start,Frame_end,Width,Height)
% Width=160;
% Height=120;
% Frame_start=1;
% Frame_end=100;
FrameNo=(Frame_end - Frame_start)+1;
VideoSeq_Input=zeros(Height,Width,FrameNo);

% VideoSeq_mov(1:FrameNo) = ...
%     struct('cdata', zeros(Width,Height,  1, 'uint8'),...
%            'colormap', gray);
writerObj = VideoWriter('ToBeEncod.avi','Uncompressed AVI');
open(writerObj);
for i=Frame_start:Frame_end
% Calculate the monochrome luminance by combining the RGB values according 
% to the NTSC standard, which applies coefficients related to the eye's 
% sensitivity to RGB colors:
 I= .2989*mov(1,i).cdata(1:Height,1:Width,1)...
	+.5870*mov(1,i).cdata(1:Height,1:Width,2)...
	+.1140*mov(1,i).cdata(1:Height,1:Width,3);
VideoSeq_Input(:,:,i)=I;
writeVideo(writerObj,I);
% VideoSeq_mov(1,i).cdata=I;
end
close(writerObj);

save('01ResizedFrames.mat','VideoSeq_Input');

VideoSeq_Obj = VideoReader('ToBeEncod.avi');
nFrames = VideoSeq_Obj.NumberOfFrames;
vidHeight = VideoSeq_Obj.Height;
vidWidth = VideoSeq_Obj.Width;
% Preallocate movie structure.
VideoSeq_mov(1:nFrames) = ...
    struct('cdata', zeros(vidHeight, vidWidth, 3, 'uint8'),...
           'colormap', []);
% Read one frame at a time.
for k = 1 : nFrames
    VideoSeq_mov(k).cdata = read(VideoSeq_Obj, k);
end
hf = figure;
set(hf, 'position', [150 150 Width Height])

% % Read one frame at a time.
% for k = 1 : nFrames
%     VideoSeq_Input(:,:)=VideoSeq_mov(1,i).cdata(:,:,1);
% end


% Play back the movie once at the video's frame rate.
movie(hf, VideoSeq_mov, 1, VideoSeq_Obj.FrameRate);
