function [mask,gui] = myAddCondStatementForFeature(gui,thisFeat,mask,source)
%%  [mask,gui] = myAddCondStatementForFeature(gui,thisFeat,mask,source)



ifButton=thisFeat.condStat;
%% Check value of ifButton
allFeatH = [gui.features.feat(:).featNum];
if ifButton.Value==0%  off
    gui.features.feat(allFeatH==thisFeat.featNum).condStat.BackgroundColor = [1 1 1]*0.7;
    guidata(gui.h0,gui);
    return;
else% turn on
    gui.features.feat(allFeatH==thisFeat.featNum).condStat.BackgroundColor = [0 1 0];
end


%% Ask user to input if statement
if isempty(ifButton.UserData)
    defaultString = 'resting==1 & adjusting==1';
else
    defaultString = ifButton.UserData;
end
if strcmp(source.Style,'togglebutton')

    myUnfoldAnnotStruct(gui);
    newStatement=inputdlg('Input if statetement to include bout as behavior (example running==1)','Add/Edit conditional statement',1,{defaultString});
    if isempty(newStatement)
        return;
    end

    newStatement = newStatement{1};
else
    newStatement = defaultString;

end
gui.features.feat(allFeatH==thisFeat.featNum).condStat.UserData = newStatement;
guidata(gui.h0,gui);
%% Extract new mask from newStatement
statements = split(newStatement,'&');
nstat = numel(statements);
nframes = numel(gui.data.annoTime);
newrast = false(1,nframes);
currChan= 'beh';% provisory
for istat = 1:nstat
    thisstat = statements{istat};
    str = split(thisstat,'==');

    label = replace(str{1},' ','');
    value = str2double(str{2});
    if isfield(gui.data.annot.(currChan),label)
    thisrast = convertToRast(gui.data.annot.(currChan).(label),nframes);
    newrast = newrast | thisrast==value;
    end

end

%% Transform existing mask

mask = mask & newrast;

% gui.annot.bhv.unsaved_feature = rast;
%% 

fprintf('   \n>>> Added if condition to unsaved_feature: %s',newStatement);



