function h=launchAnnotEditor(source,~)
%
% (C) Ann Kennedy, 2019
% California Institute of Technology
% Licensing: https://github.com/annkennedy/bento/blob/master/LICENSE.txt


gui = guidata(source);

% close previously opened manager %KM
previousAnnotEditor = findobj('tag','Annotation manager','type','figure');
if ~isempty(previousAnnotEditor)
    close(previousAnnotEditor)
end
h.fig = figure('dockcontrols','off','menubar','none',...
    'Tag','Annotation manager','name','Annotation Manager','NumberTitle','off');
gui.browser = h.fig;
p = get(h.fig,'Position');
pGUI = gui.h0.Position;
set(h.fig,'Position',[pGUI(1) pGUI(2)-330  500 300]);
myaddeditfilemenu(h.fig,'launchAnnotEditor');%KM
h.bhvMasterList = getAllBehaviors(gui);
gui = synchronizeMice(gui,h.bhvMasterList);   % get everyone to the same starting point

uicontrol('Style','text','horizontalalign','center',...
    'String','List of annotations',...
    'Units','normalized','Position',[.05 .9 .45 .05],...
    'fontsize',11);
h.bhv = uicontrol('Style','listbox',...
    'String',h.bhvMasterList,'Max',length(h.bhvMasterList),...
    'Units','normalized','Position',[.05 .1 .45 .775],...
    'fontsize',11);
h.gui = gui;
guidata(h.fig,h);

% button positioning variables
ButHeight = 0.08;
vOff = 0.1;%vertical offset
y = 0.95;

%% Play
y = y-vOff;
hbut=uicontrol('Style','pushbutton','String','Play selected',...
    'Units','normalized','Position',[.6 y .3 ButHeight],...
    'fontsize',11,'callback',{@Manager_callback,'play'});
myaddeditfilemenu(hbut,'Manager_callback');

%% Merge 
y = y-vOff;
hbut=uicontrol('Style','pushbutton','String','Merge selected',...
    'Units','normalized','Position',[.6 y .3 ButHeight],...
    'fontsize',11,'callback',{@Manager_callback,'merge'},'Tooltip','Merges within channel if two labels are selected. Merges across channels if one label selected.');
myaddeditfilemenu(hbut,'Manager_callback');

%% Split selected 
y = y-vOff;
hbut=uicontrol('Style','pushbutton','String','Split selected',...
    'Units','normalized','Position',[.6 y .3 ButHeight],...
    'fontsize',11,'callback',{@Manager_callback,'split'},'Tooltip','Splits labels in 6 sublabel based on frame similarity.');
myaddeditfilemenu(hbut,'Manager_callback');

%% Rename
y = y-vOff;
hbut=uicontrol('Style','pushbutton','String','Rename',...
    'Units','normalized','Position',[.6 y .3 ButHeight],...
    'fontsize',11,'callback',{@Manager_callback,'rename'});
myaddeditfilemenu(hbut,'Manager_callback');

%% Delete
y = y-vOff;
hbut=uicontrol('Style','pushbutton','String','Delete',...
    'Units','normalized','Position',[.6 y .3 ButHeight],...
    'fontsize',11,'callback',{@Manager_callback,'delete'});
myaddeditfilemenu(hbut,'Manager_callback');

%% 'Add new'
y = y-vOff;
hbut=uicontrol('Style','pushbutton','String','Add new',...
    'Units','normalized','Position',[.6 y .3 ButHeight],...
    'fontsize',11,'callback',{@Manager_callback,'add'});
myaddeditfilemenu(hbut,'Manager_callback');

%% Save changes
y = y-vOff*2;
hbut=uicontrol('Style','pushbutton','String','Save changes',...
    'Units','normalized','Position',[.6 y .3 .15],...
    'fontsize',11,'backgroundcolor',[.65 1 .65],...
    'callback',@Manager_applyLabelChanges);
myaddeditfilemenu(hbut,'Manager_applyLabelChanges');
