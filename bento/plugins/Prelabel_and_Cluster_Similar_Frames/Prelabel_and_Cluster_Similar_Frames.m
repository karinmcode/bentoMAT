%% Prelabel_and_Cluster_Similar_Frames
% Karin Morandell 2022 11 24
% src: edit G:\My Drive\Projects\20190620 Isometric force task\GUIs\GUI_AnalysisBehPhys_211122\code\KM video functions\get_beh_states_220608.m
function Prelabel_and_Cluster_Similar_Frames(source,varargin)

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
current_beh = gui.annot.activeBeh
ind_current_beh = gui.annot.bhv.(current_beh);
xy_current_beh = prelab_data.clu_xy(ind_current_beh,:);


%% recluster with fcm or kmeans
nclu = 6;
[Centroids,U] = fcm(xy_current_beh,nclu,[2 100 1e-5 false]);

% U , Fuzzy partition matrix, returned as a matrix with Ncenters rows and Nd columns.
% Element U(i,j) indicates the degree of membership of the jth data point in the ith cluster.
% For a given data point, the sum of the membership values for all clusters is one
maxU = max(U);
nobs = size(xy_current_beh,1);
clusterIdentifiers = nan(nobs,1);
for iclu = 1:nclu
    clusterIdentifiers(U(iclu,:) == maxU) = iclu;
end
new_annotations = cellfun(@(x) sprintf('prelabel_%s_%g',current_beh,x) ,num2cell(clusterIdentifiers),'uniformoutput',false);
new_labels = unique(new_annotations);

 %plotting
makegoodfig('fuzzy partition matrix');
CM = jet(nclu);
for iclu = 1:nclu
i4clu = clusterIdentifiers==iclu;
p(1)=plot(xy_current_beh(i4clu,1),xy_current_beh(i4clu,2),'.b');
hold on

p(2)=plot(Centroids(iclu,1),Centroids(iclu,2),'xb','MarkerSize',15,'LineWidth',3);
set(p,'color',CM(iclu,:));
end
hold off


%% Add new Behavior label with  "Prelabel_" prefix such that user can merge them with pre existing labels
% - find function that translates annotations formats to .annot
% - find how to display all new behavior labels in GUI
% - find how to edit current .annot file and reload it
newChannelName = ['prelabel_' current_beh];
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
nfr = numel(ind_current_beh);
clu_ids = nan(1,nfr);
clu_ids(ind_current_beh) =clusterIdentifiers;
pid = unique(clu_ids(~isnan(clu_ids)));
nid  = numel(pid);
annIdx = struct();
sz = 0;
for iid=1:nid%iid =1
    b = new_labels{iid};

    % add logical vector to behavior structure
    gui.annot.bhv.(b)=clu_ids==iid;

    % store annotations start and end indexes
    annIdx = convertToBouts(gui.annot.bhv.(b));
    gui.data.annot.(newChannelName).(b)=annIdx;

    % add Annotations
    nAn = size(annIdx,1);
    key = gui.annot.hotkeysDef.(b);

    gui.annot.activeBeh = b;
    
    gui.annot.modified = 1;
    gui.annot.highlightStart = [];
    gui.annot.highlighting = 0;

    for iAn = 1:nAn%iAn=3
        startIdx = annIdx(iAn,1);
        startTime = startIdx/gui.data.annoFR;
        endIdx = annIdx(iAn,2);
        endTime = endIdx/gui.data.annoFR;
        if startTime-endTime>0.5
            keyboard;
        end

        % start annot
        gui.annot.highlightStart = startIdx;
        gui.ctrl.slider.Value = startTime;% in seconds!!!
        gui = toggleAnnot(gui,'start',key);

        % end annot
        lastKey = gui.annot.highlighting;
        gui.ctrl.slider.Value = endTime;
        gui = toggleAnnot(gui,'stop',key,lastKey);

        gui.annot.modified = 1;
        gui.annot.highlightStart = [];
        gui.annot.highlighting = 0;

    end
end

updatePlot(gui.h0,[]);%update annotations case
updateSliderAnnot(gui);

guidata(gui.h0,gui);

launchAnnotEditor(gui.h0);

disp done
