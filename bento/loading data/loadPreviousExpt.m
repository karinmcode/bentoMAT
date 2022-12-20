function loadPreviousExpt(source,~)
% loadPreviousExpt(source,~)


if(strcmpi(class(source),'matlab.ui.Figure'))
    useSource = source;
else
    useSource = source.Parent.Parent;
end

% remember where the last bento file you opened was?
if myisfile("PreviousExperiment.txt")
    fid=fopen("PreviousExperiment.txt");
    txt=fscanf(fid,'%s');
    fclose(fid);
    URL = replace(txt,'''','''');
    [PathName,FileName]=fileparts(URL);
else
    [FileName,PathName] = uigetfile('*.xls;*.xlsx');
    URL = fullfile(PathName,FileName);
    fid=fopen("PreviousExperiment.txt",'w');
    fprintf(fid,'%s',URL);
    fclose(fid);
end
if(~FileName)
    return;
end

gui = guidata(useSource);
transferExptToGui(URL,gui);
%% Set some preferences here

% put gui on the most right monitor
p = get(0,'MonitorPositions');
p = p(p(:,1)==max(p(:,1)),:);%most right
gui.h0.Position(1) = p(1)+p(3)*0.01;
gui.h0.Position(2) = p(2)+p(4)-gui.h0.OuterPosition(4);

% enable some panels
gui.enabled.fineAnnot = [1 1];
gui.enabled.features = [1 1];
eventdata.Source= gui.ctrl.track.win;
updatePlot(gui.h0,eventdata);

% set prefered