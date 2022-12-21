function [mask,params] = getThresholdedFeatureMask(gui,params)
%
% (C) Ann Kennedy, 2019
% California Institute of Technology
% Licensing: https://github.com/annkennedy/bento/blob/master/LICENSE.txt



% create mask parameters from current plots
if(~exist('params','var'))
    params = [];
    for ifeat=1:length(gui.features.feat) %iterate over all plotted features (assume we're ANDing them)

        ch      = gui.features.feat(ifeat).ch;
        featNum = gui.features.feat(ifeat).featNum;

        params(ifeat).ch        = ch;
        params(ifeat).featNum   = featNum;
        params(ifeat).limL      = str2num(gui.features.feat(ifeat).threshValL.String);
        params(ifeat).limU      = str2num(gui.features.feat(ifeat).threshValU.String);
    end
end

if(isempty(params))
    mask = false(1,numel(gui.data.tracking.features{1}));
    return
end

% apply parameters to tracking data
NFrames = length(gui.data.annoTime);
vals        = gui.data.tracking.features{:,params(ifeat).featNum};%(:,params(i).featNum);

mask = true(1,numel(vals));%KM fix
for ifeat=1:length(params)
    %vals        = gui.data.tracking.features{params(i).ch}(:,params(i).featNum);
    vals        = gui.data.tracking.features{:,params(ifeat).featNum};%(:,params(i).featNum);
    vals = vals(:)';
    limL = params(ifeat).limL;
    limU = params(ifeat).limU;
    if(limL<limU)
        mask = mask & (  (vals>=limL) & (vals<=limU)  );
    else
        mask = mask & (  (vals<=limL) | (vals>=limU)  );
    end
end



% corrrect for different frame rates?
% not sure why this is here? maybe for case of annotations being different fs from features
Bouts = convertToBouts(mask);
mask = convertToRast(floor(gui.data.trackTime(Bouts)*gui.data.annoFR),NFrames);



