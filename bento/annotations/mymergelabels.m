function gui = mymergelabels(gui,oldList,newName,toKill)
%  gui = mymergelabels(gui,hAn)
% if labels are different: Merges within active channels in gui.data and gui.all_Data

anno = gui.annot;% smallcaps with logical vector format
activeCh = anno.activeCh;
ANNOT = gui.data.annot;% bout format data annot with start and end times

oldIsNew = strcmp(newName,oldList);
if oldIsNew
    keyboard;
end

%% update gui.data.annot

%% create new label in channel if not existing
if(~isfield(ANNOT.(activeCh),newName))
    ANNOT.(activeCh).(newName) = [];
end

%% concatenate bouts 
newRast = []; bhvList = [newName;oldList];
for b = 1:length(bhvList)
    newRast = [newRast; ANNOT.(activeCh).(bhvList{b})];
    if(toKill && b>1)
        ANNOT.(activeCh) = rmfield(ANNOT.(activeCh),bhvList{b});
    end
end
newRast = cleanMergedRaster(newRast);

ANNOT.(activeCh).(newName) = newRast;

gui.data.annot = ANNOT;%KM

%% update gui.allData
data    = gui.allData;
% get active m sess and trial
m    = gui.data.info.mouse;
sess = gui.data.info.session;
trial   = gui.data.info.trial;
  

% store new ANNOT;
data(m).(sess)(trial).annot = ANNOT;
gui.allData = data;%KM

%% update gui.annot.bhv with active channel values

% remove old list labels
if(toKill && b>1)
    oldList2 = intersect(oldList,fieldnames(gui.annot.bhv));%KM
    gui.annot.bhv = rmfield( gui.annot.bhv,oldList2);

end
% update gui.annot.bhv
gui=my_update_gui_annot_bhv(gui);


%% update guidata
guidata(gui.h0,gui);% important because updateLegend calls guidata
% update the legend
updateLegend(gui,1)