function [mask,gui] = myAddCondStatementForFeature(gui,thisFeat,mask,source)
%%  [mask,gui] = myAddCondStatementForFeature(gui,thisFeat,mask,source)


ifButton=thisFeat.condStat;

%% Check value of ifButton
allFeatH = [gui.features.feat(:).featNum];


if ifButton.Value==0%  off
    gui.features.feat(allFeatH==thisFeat.featNum).condStat.BackgroundColor = [1 1 1]*0.7;
    return;
else% turn on
    gui.features.feat(allFeatH==thisFeat.featNum).condStat.BackgroundColor = [0 1 0];
end


%% Ask user to input if statement
if isempty(ifButton.UserData)
    defaultString = 'grooming==0 & running==0';
else
    defaultString = ifButton.UserData;
end
if strcmp(source.Style,'togglebutton')

    myUnfoldAnnotStruct(gui);
    newStatement=inputdlg('Input if statetement to include bout as behavior (example running==1)','Add/Edit conditional statement',1,{defaultString});
    baseTooltip = 'HANDLE=feat.condStat: add conditional statement as inclusion criteria.';


    if isempty(newStatement)
        return;
    end
    newStatement = newStatement{1};
    ifButton.Tooltip= sprintf('%s\nCONDITION = %s',baseTooltip,newStatement);% update Tooltip button with statement

else
    newStatement = defaultString;

end
gui.features.feat(allFeatH==thisFeat.featNum).condStat.UserData = newStatement;

%% Extract new mask from newStatement: 
statements = split(newStatement, '&&');
statements = split(statements,'&' );
statements = split(statements,'||' );
statements = split(statements,'|' );
statements = split(statements,'AND' );
statements = split(statements,'OR' );


nstat = numel(statements);
nframes = numel(gui.data.annoTime);
currChan= 'beh';% provisory, in the future stored previous active channel somewhere
ALLRASTERS = false(nstat,nframes);
for istat = 1:nstat
    thisstat = statements{istat};
    str = split(thisstat,'==');

    label = replace(str{1},' ','');
    value = str2double(str{2});
    if isfield(gui.data.annot.(currChan),label)
        bouts = gui.data.annot.(currChan).(label);
        thisrast = convertToRast(bouts,nframes) == value;

        ALLRASTERS(istat,:)=thisrast;

    end

end

% manage multiple conditions
if contains(newStatement,{'&' 'AND'})
    newrast = all(ALLRASTERS,1);
elseif contains(newStatement,{'|' 'OR'})
    newrast = any(ALLRASTERS,1);
else
    newrast = ALLRASTERS;
end

%% Transform existing mask
mask = mask & newrast;
gui.annot.bhv.unsaved_feature = mask;

fprintf('   \n>>> Added condition to unsaved_feature: %s\n',newStatement);



