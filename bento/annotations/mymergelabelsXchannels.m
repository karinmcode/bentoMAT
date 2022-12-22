function gui=mymergelabelsXchannels(gui,h)
%  gui=mymergelabelsXchannels(gui)

labelToMerge = h.bhv.String{h.bhv.Value};

allChannels = gui.annot.channels;
currChannel = gui.annot.activeCh;
activeChInd = find(strcmp(allChannels,currChannel));
data = gui.data;

% find out which channels are relevant
pChannels = gui.ctrl.annot.ch.String;
nch = numel(pChannels);
isRelevant = false(nch,1);
for ich =1:nch
    % check out data to see if channel contains annotation of label to merge
    ch = pChannels{ich};
    chdata = data.annot.(ch);
    ISFIELD = isfield(chdata,labelToMerge);
    if ISFIELD
        NOLABELS(ich) = isempty(chdata.(labelToMerge));
    else
        NOLABELS(ich) = true;
    end
    isRelevant(ich) = ISFIELD;

end
if all(NOLABELS)
    disp('Nothing to merge.')
    return;
end

nRelevant = sum(isRelevant);
if nRelevant==1
    disp('Nothing to merge.')
    return;
end

pChannels =  pChannels(isRelevant);
nch = numel(pChannels);
myUnfoldAnnotStruct(gui);
if nRelevant>2
    CHOICE_CHANNELS = listdlg("SelectionMode","multiple",'ListString',pChannels,'PromptString','Merging across channels not coded for more than 2 channels at the time. Choose 2 channels and repeat merging across channels');
    pChannels = pChannels(CHOICE_CHANNELS);
    nch = numel(pChannels);
end
if nch>2 || nch==1
    fprintf('\nYou chose %g instead of 2 channels to merge. Merging across channels not coded for more than 2 channels at the time.',nch);
    return;
end
%% ask how we should merge
listOptions = {[ pChannels{1} ' OR ' pChannels{2}]; [pChannels{1} ' OVERWRITES ' pChannels{2}] ; [pChannels{2} ' OVERWRITES ' pChannels{1} ]; ...
    ['ADD ' pChannels{2} ' TO ' pChannels{1} ];['ADD ' pChannels{1} ' TO ' pChannels{2} ]};
CHOICE_INDEX = listdlg("SelectionMode","single",'ListString',listOptions,'PromptString','How should we merge the annotations?','ListSize',[500 200]);
disp(listOptions{CHOICE_INDEX});
%% make logical data nrows = label iterations (nchannels), ncol= nframes
nfr = numel(data.annoTime);
HITS = ones(2,nfr);
annot = data.annot;
for ich = 1:length(pChannels)
    ch_annot = annot.(pChannels{ich});

    if(strcmpi(pChannels{ich},gui.annot.activeCh))
        hits = gui.annot.bhv.(labelToMerge);
    elseif(~isempty(ch_annot.(labelToMerge)))
        idx_startend = ch_annot.(labelToMerge);
        hits  = convertToRast(idx_startend,nfr);
    else
        hits = false(1,nfr);
        
    end
    HITS(ich,:)   = hits;
end

switch CHOICE_INDEX
    case 1 %OR
        newHITs =HITS(1,:) | HITS(2,:);
        for ich = 1:nch
            HITS(ich,:) = newHITs;
        end
    case 2 %ch1>ch2
        newHITs = HITS(1,:);
        for ich = 1:nch
            HITS(ich,:) = newHITs;
        end
    case 3 %ch1<ch2
        newHITs = HITS(2,:);
        for ich = 1:nch
            HITS(ich,:) = newHITs;
        end        
    case 4 % add ch2 to ch1
        newHITs =HITS(1,:) | HITS(2,:);
        HITS(1,:) = newHITs;

    case 5% add ch1 to ch2
        newHITs =HITS(1,:) | HITS(2,:);
        HITS(2,:) = newHITs;

end


%% UPDATE DATA in bhv, data, allData

% write data in current channel bhv
gui.annot.bhv.(labelToMerge) =newHITs;

ANNOT = gui.data.annot;% bout format data annot with start and end times
for ich = 1:nch
    newHITs = HITS(ich,:);
    bouts = convertToBouts(newHITs);
    ANNOT.(pChannels{ich}).(labelToMerge)=bouts;
end

% make sure to not have double labelled frames in same channel
OPT.OverwriteOtherLabels  =1;
if OPT.OverwriteOtherLabels
    for ich = 1:nch
        newHITs = HITS(ich,:);
        ch = pChannels{ich};

        pLabels = setdiff(fieldnames(ANNOT.(ch)),labelToMerge);
        nlab = numel(pLabels);

        for ilab =1:nlab
            thislab = pLabels{ilab};
            bouts = ANNOT.(ch).(thislab);
            rast = convertToRast(bouts,nfr);
            rast(newHITs==1)=false;
            bouts = convertToBouts(rast);
            ANNOT.(ch).(thislab) = bouts;
        end
    end
end

gui.data.annot = ANNOT;
%% update gui.allData
data    = gui.allData;
% get active m sess and trial
m    = gui.data.info.mouse;
sess = gui.data.info.session;
trial   = gui.data.info.trial;
  
% store new ANNOT;
data(m).(sess)(trial).annot = ANNOT;
gui.allData = data;%KM


%gui = mymergelabels(gui,labelToMerge,newName,toKill);

% update gui.annot.bhv
gui=my_update_gui_annot_bhv(gui);


%% update guidata
guidata(gui.h0,gui);% important because updateLegend calls guidata
% update the legend
updateLegend(gui,1)

%% UPDATE GUI
updatePlot(gui.h0,[]);%KM
slider = gui.ctrl.slider;%KM
slider.updateSliderAnnot(gui);%KM


% go back to original channel
gui.annot.activeCh = currChannel;
gui.ctrl.annot.ch.Value = find(strcmp(gui.ctrl.annot.ch.String,currChannel));
setChannel(gui.ctrl.annot.ch);


myUnfoldAnnotStruct(gui);

