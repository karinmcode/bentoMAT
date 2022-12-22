function         [gui,hManager] = mySplitLabel(gui,hManager)
%%  gui = mySplitLabel(gui,hManager)
% Karin Morandell 2022 12 21

%% get current selected label
current_beh = hManager.bhv.String{ hManager.bhv.Value};
ind_current_beh = gui.annot.bhv.(current_beh);

%% Get options
% prelabel_options function outputting opt
% Process frames steps: HOG?
% Dim red steps: PCA
% Clustering steps: uMAP, kMeans, ...
edit get_prelabel_options.m;
opt.steps = {'vidmotion' 'HOG','PCA','tSNE','kmeans'};
opt.steps = {'vidmotion' 'HOG','PCA','uMAP','kmeans'};

[opt,isloading_necessary,urlstep] = get_prelabel_options(gui,opt);

%% overwrite some options
opt.SelectedFrames=find(ind_current_beh);
opt.kmeans.nclusters = 6;
opt.fcm.nclusters = 6;
nclu = 6;


%% Process data and save
% - use prewritten code

[new_labels,new_annotations,prelab_data,opt] = do_prelabelling_steps(gui,opt,isloading_necessary,urlstep);
xy_current_beh = prelab_data.clu_xy;
clusterIdentifiers = prelab_data.cluids;

new_annotations = cellfun(@(x) sprintf('prelabel_%s_%g',current_beh,x) ,num2cell(clusterIdentifiers),'uniformoutput',false);
new_labels = unique(new_annotations);

%plotting
makegoodfig('fuzzy partition matrix');
CM = jet(nclu);
for iclu = 1:nclu
    i4clu = clusterIdentifiers==iclu;
    p(1)=plot(xy_current_beh(i4clu,1),xy_current_beh(i4clu,2),'.b');
    hold on

    set(p,'color',CM(iclu,:));
end
    title(join(opt.steps));

set(gca,'colormap',CM);
cb = colorbar(gca,"eastoutside");
cbstep = 0.5/nclu;
set(cb,"Limits",[0 1],'Ticks',linspace(cbstep,1-cbstep,nclu),'TickLabels',1:nclu)
ylabel(cb,'Cluster IDs')
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

end
guidata(gui.h0,gui);

gui.enabled.legend       = [1 1];
gui.enabled.fineAnnot(1) = 1; % display fineAnnot by default
gui = redrawPanels(gui);

updatePlot(gui.h0,[]);%update annotations case
updateSliderAnnot(gui);

guidata(gui.h0,gui);


hManager=launchAnnotEditor(gui.h0);
guidata(gui.h0,gui);

disp done

