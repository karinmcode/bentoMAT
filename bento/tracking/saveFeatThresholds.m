function saveFeatThresholds(source,~)
%
% (C) Ann Kennedy, 2019
% California Institute of Technology
% Licensing: https://github.com/annkennedy/bento/blob/master/LICENSE.txt



gui = guidata(source);

[mask,params] = getThresholdedFeatureMask(gui);
if numel(gui.features.feat)==1
    activeFeatHandle = gui.features.feat;
else
    disp 'multiple features saving not coded yet! Close all the features that are irrelevant the the thresholded_feature'
    return;
end
[mask,gui] = myAddCondStatementForFeature(gui,activeFeatHandle,mask,source);% KM


if(~isfield(gui.features,'savedFilters'))
    gui.features.savedFilters{1} = params;
else
    try
    gui.features.savedFilters{end+1} = params;
    catch
        keyboard
    end
end

if(strcmpi(gui.annot.activeCh,'thresholded_features'))
    prompt = ['Save to channel:'];
    prevChan = setdiff(gui.annot.channels,'thresholded_features');
    channels = [prevChan(:); {'[create new channel]'}];
    [chStr, tf] = listdlg('promptString',prompt,'listString',channels,'SelectionMode','single');
    if ~tf
        return;
    end
    
    if chStr==length(channels)
        newCh   = addNewFieldPopup('Name of new channel:',{''});
        temp            = gui.annot.bhv.unsaved_feature; % remove the temporary feature for channel creation
        gui.annot.bhv   = rmfield(gui.annot.bhv,'unsaved_feature');
        gui             = addChannel(gui, newCh);
        
        gui.annot.activeCh = 'thresholded_features'; % flip the active channel back to where it should be
        gui.annot.bhv.unsaved_feature = temp;
        chStr = newCh;
    else
        chStr = channels{chStr};
    end
        
    str     = ['feature_filter_' num2str(length(gui.features.savedFilters))];
    prompt  = ['Name of new annotation:'];
    newStr  = inputdlg(prompt,'',1,{str});
    if(isempty(newStr))
        return;
    end
    newStr  = [newStr{:}];
    newStr  = strrep(newStr,' ','_');
    
    gui.annot.activeCh = chStr;
    gui = addLabel(gui,newStr);
	gui.annot.activeCh = 'thresholded_features';
    gui.annot.bhv = rmfield(gui.annot.bhv,newStr);
    newStr(ismember(newStr,'?!@#$%^&*()+=-<>,./\[]}{')) = [];
    
    gui.data.annot.(chStr).(newStr) = convertToBouts(mask);
    
    m = gui.data.info.mouse;
    s = gui.data.info.session;
    tr = gui.data.info.trial;
    gui.allData(m).(s)(tr).annot.(chStr).(newStr) = convertToBouts(mask);
    
    gui.enabled.legend       = [1 1];
    gui.enabled.fineAnnot(1) = 1; % display fineAnnot by default
    gui = redrawPanels(gui);
end

guidata(gui.h0,gui);
updateLegend(gui,1);
disp('added new label ')
