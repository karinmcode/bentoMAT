function setFeatThreshold(source,~,tag)
%
% (C) Ann Kennedy, 2019
% California Institute of Technology
% Licensing: https://github.com/annkennedy/bento/blob/master/LICENSE.txt



gui = guidata(source);
try
    featInd = strcmpi({gui.features.feat(:).tag},tag);
    thisFeat = gui.features.feat(featInd);
    vals = gui.data.tracking.features{thisFeat.featNum};
    thisFeat.minmax = myminmax(vals);
    lim     = thisFeat.minmax ;
catch

    keyboard
    vals = gui.features.feat(featInd).trace.YData;
    lim     = myminmax(vals);
    lim     = get(gui.features.feat.axes,'ylim');
end

switch source.Style
    case 'slider'
        newThr  = source.Value*(lim(2)-lim(1)) + lim(1);
        gui.features.feat(featInd).(['threshVal' source.Tag]).String = num2str(newThr);
    case 'edit'
        newThr = str2num(source.String);
        if(newThr<lim(1))
            newThr = lim(1);
            source.String = num2str(lim(1));
        elseif(newThr>lim(2))
            newThr = lim(2);
            source.String = num2str(lim(2));
        end
        
        gui.features.feat(featInd).(['threshold' source.Tag]).Value = ...%gui.features.feat(featInd).(['thresholdU'])
            (newThr-lim(1))/(lim(2)-lim(1));
end

updateThresholdBoxes(gui,featInd);% blue ones

%% display effects of thresholding as annotations throughout the movie
if(strcmpi(gui.annot.activeCh,'thresholded_features'))
    
    mask = getThresholdedFeatureMask(gui);

    [mask,gui] = myAddCondStatementForFeature(gui,thisFeat,mask,source);

    % and display them
    gui.annot.bhv.unsaved_feature = mask;
    guidata(gui.h0,gui);
    updateSliderAnnot(gui);
    guidata(gui.h0,gui);
    updatePlot(gui.h0,[]);
end




