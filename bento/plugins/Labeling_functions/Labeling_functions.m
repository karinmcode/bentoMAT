function Labeling_functions(source)
%
% (C) Ann Kennedy, 2019
% California Institute of Technology
% Licensing: https://github.com/annkennedy/bento/blob/master/LICENSE.txt


% create a GUI to build new labeling functions:
% gui = guidata(source);
h.fig = figure('dockcontrols','off','menubar','none',...
    'Tag','Labeling Functions','name','Create labeling function','NumberTitle','off');
gui.browser = h.fig; %I'm not actually using gui.browser right now...
p = get(h.fig,'Position');
set(h.fig,'Position',[p(1:2)-[240 0] 680 360]);

h.ctrlsBox     = uipanel('parent',h.fig,'position',[0.8 0 0.2 0.75],'bordertype','none');
h.featsBox     = uipanel('parent',h.fig,'position',[0.0125 0.0125 0.8-.0125 0.75]);
h.lblfnBox     = uipanel('parent',h.fig,'position',[0.0125 0.75+.0125 1-.025 0.25-.025],'bordertype','none');

% the features:
h.count = 0;
featList = gui.data.tracking.args{1}.features;
% no features by default

% the labeling function:
uicontrol('Parent',h.lblfnBox,'Style','edit','horizontalalign','left',...
            'String','',...
            'Units','normalized',...
            'Position',[0.1 0.2 0.8 0.6],...
            'fontsize',14,...
            'callback',{@evalFn,gui});

uicontrol('Parent',h.ctrlsBox,'Style','pushbutton','horizontalalign','left',...
            'String','Add Feature',...
            'Units','normalized',...
            'Position',[0.1 .8 0.8 .15],...
            'fontsize',10,...
            'callback',{@addFeature,featList,gui});

guidata(gcf,h);
end

function addFeature(source,~,featList,gui)
    h = guidata(source);
    h.count = h.count+1;
    varName = char(64+h.count);
    feat.rm   = uicontrol('Parent',h.featsBox,'Style','pushbutton','horizontalalign','left',...
                        'String','X','Units','normalized',...
                        'Position',[0.0125 0.95-0.125*h.count 0.05 .1],...
                        'fontsize',8,'callback',{@rmFeature,h.count});
    feat.id   = uicontrol('Parent',h.featsBox,'Style','text','horizontalalign','left',...
                        'String',varName,'Units','normalized',...
                        'Position',[0.09 0.95-0.125*h.count 0.05 .1],...
                        'fontsize',16);
    feat.menu = uicontrol('Parent',h.featsBox,'Style','popupmenu','horizontalalign','left',...
                        'String',featList,'Units','normalized',...
                        'Position',[0.15 0.95-0.125*h.count 0.4 .1],...
                        'fontsize',10,'callback',{@setFeature,gui});
    feat.ineq = uicontrol('Parent',h.featsBox,'Style','pushbutton','horizontalalign','left',...
                        'String','>','Units','normalized',...
                        'Position',[0.575 0.95-0.125*h.count 0.075 .1],...
                        'fontsize',10,'callback',@flipSign);
	feat.thr  = uicontrol('Parent',h.featsBox,'Style','edit','horizontalalign','left',...
                        'String','100','Units','normalized',...
                        'Position',[0.675 0.95-0.125*h.count 0.08 .1],...
                        'fontsize',10,'callback',@setSlider);
    feat.thrS = uicontrol('Parent',h.featsBox,'Style','slider','horizontalalign','left',...
                        'Units','normalized',...
                        'Position',[0.775 0.95-0.125*h.count 0.2 .1],...
                        'fontsize',10,'callback',@setSlider);
    
    h.feats(h.count) = feat;
    guidata(source, h);
end

function flipSign(source,~)
    if(source.String=='>')
        source.String = '<';
    else
        source.String = '>';
    end
end

function setFeature(source,~,gui)

end

function evalFn(source,~,gui)

end