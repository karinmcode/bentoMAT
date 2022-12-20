%% Prelabel_all_similar_frames
% Karin Morandell 2022 11 24
% src: edit G:\My Drive\Projects\20190620 Isometric force task\GUIs\GUI_AnalysisBehPhys_211122\code\KM video functions\get_beh_states_220608.m
function Prelabel_all_similar_frames(source)
try
    gui = guidata(source);
catch
    gui = source;
end
%% Get options 
% prelabel_options function outputting opt
% Process frames steps: HOG?
% Dim red steps: PCA
% Clustering steps: uMAP, kMeans, ...
edit get_prelabel_options.m;

[opt,isloading_necessary,urlstep] = get_prelabel_options(gui);


%% Process data and save
% - use prewritten code

[new_labels,new_annotations,prelab_data,opt] = do_prelabelling_steps(gui,opt,isloading_necessary,urlstep);

figs = do_prelabelling_figs(gui,opt,prelab_data);

%% Add new Behavior label with  "Prelabel_" prefix such that user can merge them with pre existing labels
% - find function that translates annotations formats to .annot
% - find how to display all new behavior labels in GUI
% - find how to edit current .annot file and reload it
newChannelName = 'prelabels_Karin';
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

% add new annotation

clu_ids = prelab_data.cluids(:)';
pid = unique(clu_ids);
nid  = numel(pid);


for iid=1:nid%iid =1
    label_name = new_labels{iid};

    % add logical vector to behavior structure
    i4 = clu_ids==iid;
    gui.annot.bhv.(label_name)=i4;


    % store annotations start and end indexes
    annIdx = convertToBouts(gui.annot.bhv.(label_name));
    gui.data.annot.(newChannelName).(label_name)=annIdx;
    m = gui.data.info.mouse;
    s = gui.data.info.session;
    tr = gui.data.info.trial;
    gui.allData(m).(s)(tr).annot.(newChannelName).(label_name) = annIdx;

%     % add Annotations
%     nAn = size(annIdx,1);
%     key = gui.annot.hotkeysDef.(lab);
% 
%     gui.annot.activeBeh = lab;
%     
%     gui.annot.modified = 1;
%     gui.annot.highlightStart = [];
%     gui.annot.highlighting = 0;
% 
%     for iAn = 1:nAn%iAn=3
%         startIdx = annIdx(iAn,1);
%         startTime = startIdx/gui.data.annoFR;
%         endIdx = annIdx(iAn,2);
%         endTime = endIdx/gui.data.annoFR;
%         if startTime-endTime>0.5
%             keyboard;
%         end
% 
%         % start annot
%         gui.annot.highlightStart = startIdx;
%         gui.ctrl.slider.Value = startTime;% in seconds!!!
%         gui = toggleAnnot(gui,'start',key);
% 
%         % end annot
%         lastKey = gui.annot.highlighting;
%         gui.ctrl.slider.Value = endTime;
%         gui = toggleAnnot(gui,'stop',key,lastKey);
% 
%         gui.annot.modified = 1;
%         gui.annot.highlightStart = [];
%         gui.annot.highlighting = 0;
% 
%     end
end
guidata(gui.h0,gui);

gui.enabled.legend       = [1 1];
gui.enabled.fineAnnot(1) = 1; % display fineAnnot by default
gui = redrawPanels(gui);

updatePlot(gui.h0,[]);%update annotations case
updateSliderAnnot(gui);

guidata(gui.h0,gui);
