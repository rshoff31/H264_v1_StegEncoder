%% Reading the Original Video sequence
function A01Read(Input,NumberOfFrames)
Obj = VideoReader(Input);
nFrames = NumberOfFrames;
vidHeight = Obj.Height;
vidWidth = Obj.Width;
% Preallocate movie structure.
mov(1:nFrames) = ...
    struct('cdata', zeros(vidHeight, vidWidth, 3, 'uint8'),...
           'colormap', []);
% Read one frame at a time.
for k = 1 : nFrames
    mov(k).cdata = read(Obj, k);
end
save('00OriginalVideo.mat','Obj','mov');
