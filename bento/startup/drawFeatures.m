function gui = drawFeatures(gui)
%
% (C) Ann Kennedy, 2019
% California Institute of Technology
% Licensing: https://github.com/annkennedy/bento/blob/master/LICENSE.txt

    if(isfield(gui,'tracker'))
        if(isfield(gui.features,'hullPlot'))
            delete(gui.features.hullPlot);
        end
        delete(gui.features.panel);
    end

    features.panel        = uipanel('position',[0 0 1 1],'bordertype','none');
    features.win          = 20;
    features.clickPt      = 0;
    features.axes.XLim    = [-20 20];
    features.feat         = [];
    features.featsByChannel = {};
%% features.channel
    features.channels     = uicontrol('parent',features.panel,'Style','popupmenu',...
                            'String',{''},'units','normalized','Position', [.05 0.015 .125 0.05],...
                            'Callback',@setFeatMenu);
        myaddeditfilemenu(features.channels,{'drawFeatures' 'setFeatMenu'});
%% features.menu 
    features.menu         = uicontrol('parent',features.panel,'Style','popupmenu',...
                            'String',{''},'units','normalized','Position', [.175 0.015 .34 0.05]);
    %% add feature

    uicontrol('parent',features.panel,'Style','pushbutton','String','Add',...
                            'units','normalized','Position', [.525 0.01 .075 .08],...
                            'Callback',@addFeature);

    %% add feature
   h= uicontrol('parent',features.panel,'Style','pushbutton','String','<html>Add<br/>All</html>',...
                            'units','normalized','Position', [.6 0.01 .075 .08],...
                            'Callback',@addFeature);
    myaddeditfilemenu(h,{'drawFeatures' 'addFeature'});


    % for thresholding features and making new annotations
    features.threshOn = uicontrol('parent',features.panel,'Style','pushbutton','String','Edit Thresholds',...
                            'units','normalized','Position', [.675 0.01 .25 .08],...
                            'Callback',{@toggleFeatThresholds,1});
    myaddeditfilemenu(features.threshOn,{'drawFeatures' 'toggleFeatThresholds'});

    %% hidden buttons that are shown while we're working with feature thresholds
    %% Cancel
    features.threshOff = uicontrol('parent',features.panel,'Style','pushbutton','String','Cancel',...
                            'units','normalized','Position', [.8 0.01 .15 .08],...
                            'Callback',{@toggleFeatThresholds,0},'visible','off');
    myaddeditfilemenu(features.threshOff,{'drawFeatures' 'toggleFeatThresholds'});
    %% Save button
    features.threshSave = uicontrol('parent',features.panel,'Style','pushbutton','String','Save',...
                            'units','normalized','Position', [.645 0.01 .15 .08],...
                            'Callback',@saveFeatThresholds,'visible','off');
    myaddeditfilemenu(features.threshOff,{'drawFeatures' 'saveFeatThresholds'});

    gui.features = features;
end

function setFeatMenu(source,~)
    gui = guidata(source);
    chNum = source.Value;
    if(chNum>length(gui.features.featsByChannel) || isempty(gui.features.featsByChannel{chNum}))
        return
    end
    gui.features.menu.Value = 1;
    gui.features.menu.String = gui.features.featsByChannel{chNum};
end