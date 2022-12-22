function gui = rmChannel(gui,toDelete)
%
% (C) Ann Kennedy, 2019
% California Institute of Technology
% Licensing: https://github.com/annkennedy/bento/blob/master/LICENSE.txt



if(isempty(toDelete))
    return;
end
if ischar(toDelete)
    toDelete = {toDelete};
end

newChannels              = setdiff(gui.annot.channels,toDelete);
if(isempty(newChannels))
    newChannels = {''};
end
gui.annot.channels       = newChannels;
gui.ctrl.annot.ch.String = newChannels;

% remove channel from active annotations
for ich = 1:numel(toDelete)
    ch = toDelete{ich};
    gui.data.annot      = rmfield(gui.data.annot,ch);
end

% update the drop-down menu
gui.ctrl.annot.ch.String = newChannels;
if(any(strcmpi(toDelete,gui.annot.activeCh))) % need to set a new active channel
    gui.ctrl.annot.ch.Value  = 1;
    gui.annot.activeCh  = newChannels{1};
    gui = transferAnnot(gui,gui.data);
    updateSliderAnnot(gui);
else
    gui.ctrl.annot.ch.Value = find(strcmpi(gui.annot.channels,gui.annot.activeCh));
end