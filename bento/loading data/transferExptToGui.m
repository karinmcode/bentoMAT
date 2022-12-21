function transferExptToGui(sheet_url,gui)
% data is either the path to an excel file, or is already the contents of
% that excel file (if those contents were already loaded using the data
% builder gui)
%
% (C) Ann Kennedy, 2019
% California Institute of Technology
% Licensing: https://github.com/annkennedy/bento/blob/master/LICENSE.txt



if(iscell(sheet_url))
    raw = sheet_url;
    if(isempty(raw{1,1}))
        raw{1,1} = uigetdir(pwd,'Please provide the path to your data''s parent directory!');
    end
else
    try
        if exist(sheet_url,'file')==2
            [~,~,raw]   = xlsread(sheet_url,'Sheet1');
        else
            sheet_url = replace(sheet_url,'MyDrive','My Drive');%KM
            [~,~,raw]   = xlsread(sheet_url,'Sheet1');
        end

    catch
        [~,~,raw]   = xlsread([sheet_url 'x'],'Sheet1');% KM
    end
    if(isempty(raw{1,1})|isnan(raw{1,1}))
        pth = fileparts(sheet_url);
        raw{1,1} = [pth filesep]; %blank path means read from the directory the sheet is in
    end
end
[allData,enabled,pth,hotkeys,params] = unpackExperiment(raw);
gui.pth = pth;
gui.mydatafolder =params.My_Data_Folder;%KM
gui.expt_sheet = sheet_url;%KM
gui.annot = mergeHotkeys(gui.annot,hotkeys);

% toggle window visiblity
gui.enabled         = enabled;
gui.enabled.welcome = [0 0];
gui.enabled.ctrl    = [1 1];

gui.traces.toPlot = 'rast'; %default plot type for  cell traces/image

if(isfield(gui,'data')), gui = rmfield(gui,'data'); end

gui.allData         = allData;                %stores all mice!
gui.allPopulated    = cell2mat(raw(3:end,1:3)); %keeps a list of which mouse/sess/trials are populated
mouseList           = unique(gui.allPopulated(:,1));
use                 = gui.allPopulated(:,1) == mouseList(1);
guidata(gui.h0,gui); % just in case it crashes after this

% set the active mouse/session/trial
sessList    = cellstr(num2str(unique(gui.allPopulated(use,2))));
trialList   = strtrim(cellstr(num2str(unique(gui.allPopulated(use,3)))));
gui         = setActiveMouse(gui,mouseList(1),['session' sessList{1}],str2num(trialList{1}),1);

% update the control bar that lets us browse through mice, sesions, and
% trials:
set(gui.ctrl.expt.mouse,   'String',cellstr(num2str(mouseList)));
set(gui.ctrl.expt.session, 'Enable','on','String',sessList);
set(gui.ctrl.expt.trial,   'Enable','on','String',trialList);

% update everything
gui = redrawPanels(gui);%  figure out which panels are enabled, and set their visibility
guidata(gui.h0,gui);
updatePlot(gui.h0,[]); % update values in panels

