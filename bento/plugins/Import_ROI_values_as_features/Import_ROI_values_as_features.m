%% Import_ROI_values_as_features
% Karin Morandell 2022 11 30
% src: edit 'G:\My Drive\Projects\20190620 Isometric force task\GUIs\GUI_AnalysisBehPhys_211122\code\KM video functions\get_beh_states_220608.m'
function Import_ROI_values_as_features(source)


%% Get variables

%gui
try
    gui = guidata(source);
catch
    gui = source;
end

% video
[~,VidReader] = mygcv(gui);

% folders
info = findmydata(gui);
url_save = info.url.features;


edit Import_ROI_values_as_features.m;




%% Set ROI value


list = {'snout' 'frame' 'mouth' 'whiskers' 'wheel' 'other' 'keyboard'};
ANSidx=listdlg('ListString',list);
if isempty(ANSidx)
    return
end

if numel(ANSidx)>1
disp('code for multiple ROIs at once not written yet')
    keyboard
    return
end
featureName = list{ANSidx};

if ismember(featureName,{'other'})
    featureName=inputdlg('Input name of the new feature.');
    featureName = featureName{1};
end
if ismember(featureName,{ 'keyboard'})
    keyboard
end



%% NORMAL CASES
varname.ROI_dim = sprintf('ROI_%s_dim',featureName);
varname.ROI_values = sprintf('ROI_%s_values',featureName);

% if exist(url_save,'file')==2
%     load(url_save,varname.ROI_values);
% else

sz=0;
% draw ROI
VidReader.CurrentTime=0;
%ROI_dim = [670.7274  288.5000   28.5481   41.9825 ];%provisory
if ~exist('ROI_dim','var')
    ftemp = makegoodfig('Draw ROI');
    ax = axes(ftemp);

    frame = VidReader.readFrame();
    im=imagesc(ax,frame);
    im.CDataMapping = "direct";
    axis(ax,'equal','tight');
    title(ax,featureName)
    ROI_dim= drawrectangle(ax);
    wait(ROI_dim);
    if ~isvalid(ftemp)
        return
    end
    ROI_dim =ROI_dim.Position;
    close(ftemp);
end

ROI_values = nan(1,VidReader.NumFrames);
VidReader.CurrentTime = 0;
computationOpt = {'abs(diff(red - green))' 'sum all pixels across channels' 'sum red channel' 'sum green channel' 'sum red and green channel' 'min value across ROI pixels' 'max value across ROI pixels' 'sum white pixels'  'other'};
choiceIdx = listdlg('ListString',computationOpt,'SelectionMode','single');
Choice = computationOpt{choiceIdx};
%% GO ACROSS FRAMES
for ifr = 1:VidReader.NumFrames

    if isnan(ROI_values(ifr))% compute ME only for frames that have not been compute before

        cleanline(sz);
        sz=fprintf('%g/%g  (%.0f%%) : DO NOT TOUCH GUI WHILE THIS IS RUNNING (VIDEOR READER GETS TRIGGERED BY GUI)',ifr,VidReader.NumFrames,100*ifr/VidReader.NumFrames);
        frame = VidReader.readFrame;
        switch Choice
            case 'sum all pixels across channels'
                I = imcrop(frame,ROI_dim);
                ROI_values(ifr) =sum(I(:));
            case 'sum red channel'
                I = imcrop(frame(:,:,1),ROI_dim);
                ROI_values(ifr) =sum(I(:));
            case 'sum green channel'
                I = imcrop(frame(:,:,2),ROI_dim);
                ROI_values(ifr) =sum(I(:));        
            case 'sum red and green channel'
                I = imcrop(frame(:,:,2),ROI_dim);
                I = imcrop(frame(:,:,1),ROI_dim)+I;
                ROI_values(ifr) =sum(I(:)); 

            case 'abs(diff(red - green))'
                I = imcrop(frame,ROI_dim);
                I = abs(diff(I(:,:,1:2),1,3));
                ROI_values(ifr) =mean(I(:));
            case 'min value across ROI pixels'
                keyboard
            case 'max value across ROI pixels'
                keyboard
            case 'sum white pixels'
                I = imcrop(frame,ROI_dim);
                I = I/255;
                I =round(double(I),1);
                I = all(I==1,3);
                ROI_values(ifr) =sum(I(:));
 
            otherwise
                keyboard
        end
    end
end

ROI_values =zscore(ROI_values);
%ROI_values = -ROI_values;
if startsWith(Choice,'sum')
ROI_values = ROI_values-mymin(ROI_values);
end
%url = fullfile(V.Path,V.Name);
eval(sprintf('%s=ROI_values;',varname.ROI_values));
eval(sprintf('%s=ROI_dim;',varname.ROI_dim));

if exist(url_save,'file')==0
    mymkdir(url_save);
    save(url_save,varname.ROI_values);%load(url_save)winopen(folder)
    
    %% update experiment sheet
    params.Tracking = url_save;
    my_set_current_expt_params(gui,params);
else
    save(url_save,varname.ROI_values,'-append');%load(url_save)winopen(folder)
end

%% update data and gui menu
newFeatures.(varname.ROI_values)=ROI_values;

myImportNewFeatures(gui,newFeatures)

% end

% % color low motion energy
% MEmat=cell2mat(ROI_values);
% MEavg = mean(MEmat);
% MEstd = std(MEmat);
% grp_edges = [fliplr(MEavg:-MEstd:0)  MEavg:MEstd:10*MEstd];
% MEmat(MEmat>max(grp_edges))=max(grp_edges);
% MEmat(MEmat<min(grp_edges))=min(grp_edges);
% ngr = numel(grp_edges)-1;
% CM= cool(ngr)+0.2;
% CM(CM>1)=1;
% grp_ME = discretize(MEmat,grp_edges);
% tab.BackgroundColor = CM(grp_ME,:);
%
% % hist of ME
% fme = makegoodfig('ME');
% ax = axes(fme);
% histogram(ax,MEmat,grp_edges);
%
% %% Get ROI Value
%
% %% Set ROI threshold
%
% %% Add new Behavior label with  "Prelabel_ROI" prefix such that user can merge them with pre existing labels
% % - find function that translates annotations formats to .annot
% % - find how to display all new behavior labels in GUI
% % - find how to edit current .annot file and reload it
% newChannelName = 'prelabel_ROI';
% ChannelExists = ismember(newChannelName,gui.ctrl.annot.ch.String);
% if ChannelExists
%     toDelete = {newChannelName};
%     gui      = rmChannel(gui,toDelete);
%     gui.ctrl.annot.ch.Value = 1;
% end
% gui = addChannel(gui,newChannelName );
% gui.enabled.fineAnnot(strcmp(gui.annot.channels,newChannelName))=1;
% existingLabels = fieldnames(gui.annot.bhv);
%
% gui = addLabel(gui,'ROI1','');
%
% gui.enabled.annot        = [1 1]; % enable annots if they haven't been already
% gui.enabled.legend       = [1 1];
% gui.enabled.fineAnnot(1) = 1; % display fineAnnot by default?
% gui = redrawPanels(gui);
%
% updateLegend(gui,1);
% guidata(gui.h0,gui);
% updatePlot(gui.h0,[]);%
%
% % add new annotation
%
% pid = 1;
% nid  = numel(pid);
% annIdx = struct();
% for iid=1%iid =1
%     b = ROI1;
%
%     % add logical vector to behavior structure
%     gui.annot.bhv.(b)=new_ids==iid;
%
%     % store annotations start and end indexes
%     annIdx = convertToBouts(gui.annot.bhv.(b));
%     gui.data.annot.(newChannelName).(b)=annIdx;
%
%     % add Annotations
%     nAn = size(annIdx,1);
%     key = gui.annot.hotkeysDef.(b);
%
%     gui.annot.activeBeh = b;
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
% end
% updatePlot(gui.h0,[]);%update annotations case
% updateSliderAnnot(gui);
%
% guidata(gui.h0,gui);
