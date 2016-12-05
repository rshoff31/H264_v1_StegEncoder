function A01SaveAndPlay(InputFile,Frames)
% Crate writer object
writerObj = VideoWriter(InputFile,'Uncompressed AVI');
% Open it
open(writerObj);
[Height, Width, Frame_end]=size(Frames);
% Store the Fremas in it
for i=1:Frame_end
    I=uint8(Frames(:,:,i));
    writeVideo(writerObj,I);
end
% Close/Store it
close(writerObj);

% Play it
VideoSeq_Obj = VideoReader(InputFile);
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
set(hf, 'position', [150 150 vidHeight vidWidth])

% Play back the movie once at the video's frame rate.
movie(hf, VideoSeq_mov, 1, VideoSeq_Obj.FrameRate);