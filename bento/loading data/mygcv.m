function [currentVideoUrl,VidReader] = mygcv(gui)
%% currentVideoUrl = mygcv(gui)
% MYDATAFOLDER = gui.mydatafolder;
currentVideoUrl=gui.data.io.movie.fid{1};% provisory
VidReader = gui.data.io.movie.reader{1}.reader;