function addFeature(source,~)
%
% (C) Ann Kennedy, 2019
% California Institute of Technology
% Licensing: https://github.com/annkennedy/bento/blob/master/LICENSE.txt


gui = guidata(source);

if strcmpi(source.String,'Add')
    i4CurrentFeat = gui.features.menu.Value;%gui.features.menu.String
    i4CurrentCh = gui.features.channels.Value;
else
    i4CurrentFeat = 1:size(gui.features.menu.String,1);
    i4CurrentCh   = gui.features.channels.Value*ones(1,size(gui.features.menu.String,1));
end

for i = 1:length(i4CurrentFeat)
    feat = struct();
    feat.featNum      = i4CurrentFeat(i);
    feat.ch           = i4CurrentCh(i);
    %tag = [gui.features.channels.String{gui.features.channels.Value} ': ' strtrim(gui.features.menu.String(feat.featNum,:))];% ORIGINAL CODE
    tag = [gui.features.channels.String{gui.features.channels.Value} ': ' strtrim(gui.features.menu.String{feat.featNum})];
    if(~isempty(gui.features.feat))
        if(any(strcmpi(tag,{gui.features.feat.tag})))
            return;
        end
    end

    % create a new feature object
    feat.axes         = axes('parent',gui.features.panel,'ytick',[]); hold on;
    feat.rmBtn        = uicontrol('parent',gui.features.panel,'Style','pushbutton',...
                            'String','X','units','normalized','Position', [.5 .5 .1 .05],...
                            'callback',{@rmFeature,tag});

    %% thresholdU slider
    feat.thresholdU    = uicontrol('parent',gui.features.panel,'Style','slider','tag','U',...
                            'units','normalized','position',[0 0 .05 .1],...
                            'min',0,'max',1,'value',1,...
                            'visible',gui.features.threshOff.Visible,...
                            'callback',{@setFeatThreshold,tag},'Tooltip','upperThreshold (handle = feat.thresholdU)');
    %% thresholdL slider
    feat.thresholdL    = uicontrol('parent',gui.features.panel,'Style','slider','tag','L',...
                            'units','normalized','position',[0 0 .05 .1],...
                            'min',0,'max',1,'value',0.25,...
                            'visible',gui.features.threshOff.Visible,...
                            'callback',{@setFeatThreshold,tag},'Tooltip','lower Threshold (handle = feat.thresholdL)');
    feat.thresholdL.BackgroundColor = [0.9 0.9 1];% KM

    feat.threshLineU   = patch([0 0 0 0],[0 0 0 0],[.25 .75 1],'facealpha',.35,'edgecolor','none',...
                            'hittest','off','visible',gui.features.threshOff.Visible);

    feat.threshLineL   = patch([0 0 0 0],[0 0 0 0],[.25 .75 1],'facealpha',.35,'edgecolor','none',...
                            'hittest','off','visible',gui.features.threshOff.Visible);
    %% thresholdU edit
    feat.threshValU   = uicontrol('parent',gui.features.panel,'Style','edit','tag','U',...
                            'units','normalized','position',[.5 .5 .1 .05],...
                            'visible',gui.features.threshOff.Visible,...
                            'callback',{@setFeatThreshold,tag},'Tooltip','upperThreshold (handle = feat.threshValU)');


    %% thresholdL edit
    feat.threshValL   = uicontrol('parent',gui.features.panel,'Style','edit','tag','L',...
                            'units','normalized','position',[.5 .5 .1 .05],...
                            'visible',gui.features.threshOff.Visible,...
                            'callback',{@setFeatThreshold,tag},'Tooltip','lowThreshold  (handle = feat.threshValL)');
    feat.thrBound     = 1;

    feat.tag          = tag;
    feat.label        = text(feat.axes,-gui.features.win*0.975,...
                             mean(reshape(gui.data.tracking.features{feat.ch}(:,feat.featNum),1,[])),...
                             strrep(tag,'_',' '),'fontsize',12);

    feat.trace        = plot(0,0,'color',[.1 .1 .1],'hittest','off');
    feat.zeroLine     = plot([0 0],get(gca,'ylim'),'k--','hittest','off');
    feat.bg           = image(ones(1,1,3),'hittest','off');
    uistack(feat.trace,'top');%gui.features.feat.trace;
    uistack(feat.bg,'bottom');
    feat.axes.ButtonDownFcn = {@figBoxCheck,'features'};
    myaddeditfilemenu(feat.axes,{'addFeature' 'getThresholdedFeatureMask'});
    xlim(gui.features.win*[-1 1]);

    %% KM  conditions statement togglebutton
    feat.condStat        = uicontrol('parent',gui.features.panel,'Style','togglebutton',...
                            'String','+if','units','normalized','Position', [.5 .4 .1 .05],...
                            'callback',{@setFeatThreshold,tag},'tooltipstring','feat.condStat: add conditional statement as inclusion criteria.','UserData','','Visible','off');
    myaddeditfilemenu(feat.condStat,{'addFeature' 'drawFeature' 'setFeatThreshold' 'myAddCondStatementForFeature' 'redrawFeaturePlots'});
    %add this feature to the list
    if(isempty(gui.features.feat))
        gui.features.feat        = feat;
    else
        gui.features.feat(end+1) = feat;
    end
end

gui = redrawFeaturePlots(gui);

guidata(gui.h0,gui);
updatePlot(gui.h0,[]);
