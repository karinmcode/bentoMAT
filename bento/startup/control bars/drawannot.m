function annot = drawannot(gui,row,scale,nRows)
%
% (C) Ann Kennedy, 2019
% California Institute of Technology
% Licensing: https://github.com/annkennedy/bento/blob/master/LICENSE.txt


%% panel
annot.panel = uipanel('parent',gui.ctrl.panel,...
        'position',[0.01 (row-.5)/(nRows+1) 0.98 scale/(nRows+1)],'bordertype','none');

%% x locations of GUI elements
x= .005;
xstep =0.002;

%% Channel text
h=uicontrol('parent',annot.panel,'Style','text',...
            'String','Channel','units','normalized','Position', [x 0 .065 .7]);

%% Channel popupmenu
x = x+h.Position(3)+xstep;
annot.ch = uicontrol('parent',annot.panel,'Style','popupmenu',...
            'String',{''},'units','normalized','Position', [x 0 .15 .8],...
            'Callback',@setChannel);
myaddeditfilemenu(annot.ch,{'drawannot' 'setChannel' }  );
h = annot.ch;
%% 
x = x+h.Position(3)+xstep;
h=uicontrol('parent',annot.panel,'Style','pushbutton',...
            'String','Add','units','normalized','Position', [x 0 .065 .8],...
            'Tag','add channel','Callback',@annotButtonHandler);
myaddeditfilemenu(h,{'drawannot' 'annotButtonHandler' });

%% 
x = x+h.Position(3)+xstep;
h=uicontrol('parent',annot.panel,'Style','pushbutton',...
            'String','Delete','units','normalized','Position', [x 0 .065 .8],...
            'Tag','remove channel','Callback',@annotButtonHandler);
myaddeditfilemenu(h,{'drawannot' 'annotButtonHandler' });

%% Duplicate
x = x+h.Position(3)+xstep;
h=uicontrol('parent',annot.panel,'Style','pushbutton',...
            'String','Duplicate','units','normalized','Position', [x 0 .065 .8],...
            'Tag','copy channel','Callback',@annotButtonHandler);
myaddeditfilemenu(h,'drawannot');
myaddeditfilemenu(h,'annotButtonHandler');
myaddeditfilemenu(h,'copyChannel');

%% Rename
x = x+h.Position(3)+xstep;
h=uicontrol('parent',annot.panel,'Style','pushbutton',...
            'String','Rename','units','normalized','Position', [x 0 .065 .8],...
            'Tag','rename channel','Callback',@annotButtonHandler);
myaddeditfilemenu(h,{'drawannot' 'annotButtonHandler' 'renameChannel'});


%%  Labels 
x = x+h.Position(3)+xstep;
h=uicontrol('parent',annot.panel,'Style','text',...
            'String','Labels','units','normalized','Position', [x 0 .08 .7]);
myaddeditfilemenu(h,'drawannot');

%% Add
x = x+h.Position(3)+xstep;
h=uicontrol('parent',annot.panel,'Style','pushbutton',...
            'String','Add','units','normalized','Position', [x 0 .065 .8],...
            'Tag','add behavior','Callback',@annotButtonHandler);
myaddeditfilemenu(h,{'drawannot' 'annotButtonHandler' });

%% Delete behavior button
x = x+h.Position(3)+xstep;
h=uicontrol('parent',annot.panel,'Style','pushbutton',...
            'String','Delete','units','normalized','Position', [x 0 .065 .8],...
            'Tag','remove behavior','Callback',@annotButtonHandler);
myaddeditfilemenu(h,{'drawannot' 'annotButtonHandler' });

%% Fast Edit button
x = x+h.Position(3)+xstep;
annot.fastEdit = uicontrol('parent',annot.panel,'Style','togglebutton',...
            'String','Fast Edit','units','normalized','Position', [x 0 .065 .8],...
            'Tag','edit');
h = annot.fastEdit;
myaddeditfilemenu(h,{'drawannot' 'annotButtonHandler' });

%% Save button
x = x+h.Position(3)+xstep*2;
load('icons.mat');%#ok
annot.save        = uicontrol(annot.panel,'Style','pushbutton','CData',imgs.f3,'units','normalized',...
               'position',[x .05 .05 .8],'tooltip','Save annotations','callback',@quickSave);
myaddeditfilemenu(annot.save,{'drawannot' 'quickSave'} );

%% MARS button
annot.MARS        = uicontrol(annot.panel,'Style','pushbutton','string','MARS','units','normalized',...
               'position',[.85 0 .1 1],'backgroundcolor',[1 .75 .85],'callback',@MARSoptions);
myaddeditfilemenu(annot.MARS,{'drawannot' 'MARSoptions'});
