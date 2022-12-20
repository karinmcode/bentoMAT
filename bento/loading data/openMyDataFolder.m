function openMyDataFolder(h0)
% openMyDataFolder(h0)
gui = guidata(h0);

currentVideoUrl = mygcv(gui);
fo = fileparts(currentVideoUrl);
mywinopen(fo);