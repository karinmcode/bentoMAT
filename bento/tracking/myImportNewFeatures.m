function gui=myImportNewFeatures(gui,newFeatures)
%% gui2=myImportNewFeatures(gui1,newFeatures)
% function:
% - add features to data
% - adds new features to gui menu
% INPUTS:
% newFeatures: structure with fields being new feature names

%% INITIALIZE data.TRACKING
if ~isfield(gui.data,'tracking')
    %edit unpackExperiment.m
    gui.data.tracking=struct();
    gui.data.tracking.args = cell(1,1);% cell with structure inside
    gui.data.tracking.features = cell(1,1);
    gui.data.tracking.active = cell(1,1);
    gui.data.tracking.inactive = cell(1,1);
    gui.data.tracking.crop = [];
end

fNames = fieldnames(newFeatures);
nfe = numel(fNames);
%% LOOP NEW FEATURES
for ife = 1:nfe
    fName = fNames{ife};
    fValues = newFeatures.(fName);

    %%    - add to gui data
    channelNb = 1;% not sure how channel numbers are defined

    gui.data.tracking.args{channelNb}.(fName)=fValues;

    allFeatures = fieldnames(gui.data.tracking.args{channelNb});

    i4feat = find(strcmp(allFeatures,fName));

    gui.data.tracking.features{i4feat}=fValues;%#ok

    gui.data.tracking.fun = 'generic_timeseries';% provisory


    %%    - add to gui menu
    if ~ismember(fName,gui.features.menu.String)
        gui.features.menu.String = vertcat(gui.features.menu.String,fName);
    end



end

%% Update all data
updateAllData(gui);

%% Update guidata in case
guidata(gui.h0,gui);

%% Update GUI elements
gui.enabled.features = [1 1];
gui = redrawPanels(gui);
redrawFeaturePlots(gui);
guidata(gui.h0,gui);

