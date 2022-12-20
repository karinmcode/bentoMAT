function Import_events_as_annotation_channel(source)
% Import_events_as_annotation_channel(gcf)

try
    gui = guidata(source);
catch
    gui = source;
end

%% Get variables
info = findmydata(gui);

folder_proc = info.fo.vid;
url_events = info.url.events;
load(url_events,'events');%mywinopen(url_position)
newChannelName = 'events';
new_labels = fieldnames(events);

%% Add channel 
ChannelExists = ismember(newChannelName,gui.ctrl.annot.ch.String);
if ChannelExists
    toDelete = {newChannelName};
    gui      = rmChannel(gui,toDelete);
    gui.ctrl.annot.ch.Value = 1;
end
gui = addChannel(gui,newChannelName );
gui.enabled.fineAnnot(strcmp(gui.annot.channels,newChannelName))=1;
existingLabels = fieldnames(gui.annot.bhv);
for newlab = new_labels(:)'
    if ~ismember(newlab{:},existingLabels)
    gui = addLabel(gui,newlab{1},'');
    end
end

gui.enabled.annot        = [1 1]; % enable annots if they haven't been already
gui.enabled.legend       = [1 1];
gui.enabled.fineAnnot(1) = 1; % don't display fineAnnot by default?
gui = redrawPanels(gui);

updateLegend(gui,1);
guidata(gui.h0,gui);
updatePlot(gui.h0,[]);% 

%% create logical vector for gui.annot.bhv and add new annotation times to gui.data
nlab = numel(new_labels);
timestampsAnnot =gui.data.annoTime;
nTimestamps = numel(timestampsAnnot);
% load 
load(fullfile(folder_proc,sprintf("vid_CAM%s.mat",info.camID)),'video');%provisory

for ila = 1:nlab
    lab = new_labels{ila};
    if isfield(events,'time_vid_proc')
        times = events.(lab).time_vid_proc;
    else
        times = events.(lab).time- video.time_proc(1);% align to the begining of the processed video time
    end
    if any(times<=0)
        keyboard
    end
    t1 = times;
    if isfield(events.(lab),'time_end')
        t2 = events.(lab).time_end;
    elseif isfield(events.(lab),'dur')
        event_duration_s = events.(lab).dur;
        t2 = event_duration_s+t1;
    else
        t2 = t1+1;% assuming event was 1 sec long
    end
    dur = t2-t1;
    Bouts = [floor(t1*gui.data.annoFR) floor(t1*gui.data.annoFR)+round(dur*gui.data.annoFR)];%FLOOR t1 !!!
    rast = convertToRast(Bouts,numel(gui.data.annoTime));

    % update bhv
    gui.annot.bhv.(lab)=rast;
    gui.data.annot.(newChannelName).(lab)=Bouts;

end
gui.enabled.fineAnnot(gui.enabled.fineAnnot==0)=1;
gui.ctrl.slider.updateSliderAnnot(gui);

guidata(gui.h0,gui);

%% update slider annotation image

updatePlot(gui.h0,[]);%update annotations case


end




