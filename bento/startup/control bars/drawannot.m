function annot = drawannot(gui,row,scale,nRows)
%
% (C) Ann Kennedy, 2019
% California Institute of Technology
% Licensing: https://github.com/annkennedy/bento/blob/master/LICENSE.txt


%% panel
annot.panel = uipanel('parent',gui.ctrl.panel,...
        'position',[0.01 (row-.5)/(nRows+1) 0.98 scale/(nRows+1)],'bordertype','none');

%% Channel text
uicontrol('parent',annot.panel,'Style','text',...
            'String','Channel','units','normalized','Position', [.025 0 .075 .7]);

%% Channel popupmenu
annot.ch = uicontrol('parent',annot.panel,'Style','popupmenu',...
            'String',{''},'units','normalized','Position', [.1 0 .15 .8],...
            'Callback',@setChannel);
myaddeditfilemenu(annot.ch,'drawannot');
myaddeditfilemenu(annot.ch,'setChannel');

%% 
h=uicontrol('parent',annot.panel,'Style','pushbutton',...
            'String','Add','units','normalized','Position', [.255 0 .065 .8],...
            'Tag','add channel','Callback',@annotButtonHandler);
myaddeditfilemenu(h,'drawannot');
myaddeditfilemenu(h,'annotButtonHandler');

%% 
h=uicontrol('parent',annot.panel,'Style','pushbutton',...
            'String','Delete','units','normalized','Position', [.32 0 .065 .8],...
            'Tag','remove channel','Callback',@annotButtonHandler);
myaddeditfilemenu(h,'drawannot');
myaddeditfilemenu(h,'annotButtonHandler');

%%
h=uicontrol('parent',annot.panel,'Style','pushbutton',...
            'String','Duplicate','units','normalized','Position', [.385 0 .065 .8],...
            'Tag','copy channel','Callback',@annotButtonHandler);
myaddeditfilemenu(h,'drawannot');
myaddeditfilemenu(h,'annotButtonHandler');
myaddeditfilemenu(h,'copyChannel');

  %%      
h=uicontrol('parent',annot.panel,'Style','text',...
            'String','Behaviors','units','normalized','Position', [.465 0 .1 .7]);
myaddeditfilemenu(h,'drawannot');

%%
h=uicontrol('parent',annot.panel,'Style','pushbutton',...
            'String','Add','units','normalized','Position', [.55 0 .065 .8],...
            'Tag','add behavior','Callback',@annotButtonHandler);
myaddeditfilemenu(h,'drawannot');
myaddeditfilemenu(h,'annotButtonHandler');

%% Delete behavior button
h=uicontrol('parent',annot.panel,'Style','pushbutton',...
            'String','Delete','units','normalized','Position', [.615 0 .065 .8],...
            'Tag','remove behavior','Callback',@annotButtonHandler);
myaddeditfilemenu(h,'drawannot');
myaddeditfilemenu(h,'annotButtonHandler');

%% Fast Edit button
annot.fastEdit = uicontrol('parent',annot.panel,'Style','togglebutton',...
            'String','Fast Edit','units','normalized','Position', [.68 0 .065 .8],...
            'Tag','edit');
myaddeditfilemenu(annot.fastEdit,'drawannot');
myaddeditfilemenu(annot.fastEdit,'annotButtonHandler');

%% Save button
load('icons.mat');%#ok
annot.save        = uicontrol(annot.panel,'Style','pushbutton','CData',imgs.f3,'units','normalized',...
               'position',[.79 .05 .05 .8],'tooltip','Save annotations','callback',@quickSave);
myaddeditfilemenu(annot.save,'drawannot');
myaddeditfilemenu(annot.save,'quickSave');

%% MARS button
annot.MARS        = uicontrol(annot.panel,'Style','pushbutton','string','MARS','units','normalized',...
               'position',[.85 0 .15 1],'backgroundcolor',[1 .75 .85],'callback',@MARSoptions);
myaddeditfilemenu(annot.MARS,'drawannot');
myaddeditfilemenu(annot.MARS,'MARSoptions');