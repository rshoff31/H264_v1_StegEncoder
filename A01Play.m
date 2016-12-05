%% Playing
function A01Play(mov,Obj)

hf = figure;
set(hf, 'position', [150 150 Obj.Width Obj.Height])

% Play back the movie once at the video's frame rate.
movie(hf, mov, 1, Obj.FrameRate);