function gui = addLabel(gui,newStr,varargin)
%
% (C) Ann Kennedy, 2019
% California Institute of Technology
% Licensing: https://github.com/annkennedy/bento/blob/master/LICENSE.txt

if numel(varargin)==1%KM
    newhotkey= varargin{1};
    if isempty(newhotkey)% generate random new hotkey
        existingHotkeys = cell2mat(fieldnames(gui.annot.hotkeys)');
        allHotkeys ='a':'z' ;
        possibleHotkeys = setdiff(allHotkeys,existingHotkeys);
        if isempty(possibleHotkeys)
            allHotkeys ='A':'Z' ;
            possibleHotkeys = setdiff(allHotkeys,existingHotkeys);
        end
        if isempty(possibleHotkeys)
            allHotkeys ='1':'8' ;
            possibleHotkeys = setdiff(allHotkeys,existingHotkeys);
        end

        if isempty(possibleHotkeys)
            %keyboard
            allHotkeys =mynum2str((1:12)','F%g','cellstr') ;
            possibleHotkeys = '-';
        end
        newhotkey = randsample(possibleHotkeys,1);

    end
end

if(isempty(newStr))
    return;
end
if(isfield(gui.annot.bhv,newStr))
    uiwait(msgbox(['A behavior called ' newStr ' already exists!']));
    return
end


% save the new label:
gui.annot.bhv.(newStr)  = false(size(gui.data.annoTime));
gui.annot.show.(newStr) = 1;
gui.annot.modified      = 1;


% set the new annotation color:
if(any(strcmpi(fieldnames(gui.annot.cmapDef),newStr)))
    newColor = gui.annot.cmapDef.(newStr);
else
    used     = [1 1 1; 0 0 0; .94 .94 .94];
    for f = fieldnames(gui.annot.cmap)'
        used(end+1,:) = gui.annot.cmap.(f{:});
    end
    newColor = distinguishable_colors(1,used)/2+.5;
end
gui.annot.cmap.(newStr)     = newColor;
gui.annot.cmapDef.(newStr)  = newColor;


% set the new hotkey:
if(any(strcmpi(fieldnames(gui.annot.hotkeysDef),newStr)) && ~strcmpi(gui.annot.hotkeysDef.(newStr),'_'))
    hotkey = gui.annot.hotkeysDef.(newStr);
    if(strcmpi(newStr,'unsaved_feature'))
        hotkey='z';
    end
    gui.annot.hotkeysDef.(newStr) = hotkey;
    gui.annot.hotkeys.(hotkey) = newStr;

elseif(~strcmpi(newStr,'unsaved_feature'))
    if exist('newhotkey','var')
        hotkey = {newhotkey};
    else
        hotkey = inputdlg(['Assign hotkey for ' strrep(newStr,'_',' ') '?']);% hotkey assignment
    end
    if(~isempty(hotkey))
        hotkey = regexprep(hotkey{:},'[^a-zA-z]','');
        if(length(hotkey)==1)
            gui.annot.hotkeys.(hotkey)    = newStr; % add the hotkey to the list


            for f = fieldnames(gui.annot.hotkeysDef)' %unassign that hotkey from other behaviors
                if(gui.annot.hotkeysDef.(f{:}) == hotkey) && hotkey~='-'
                    %fprintf('\n>unassign that hotkey from other behaviors: %s',f{:})
                    gui.annot.hotkeysDef.(f{:}) = '_';
                end
            end
            gui.annot.hotkeysDef.(newStr) = hotkey;

        end
    end
end


% add the new label to current data:
gui.data.annot.(gui.annot.activeCh).(newStr) = [];


% add the new label to gui.allData:
sessionList = fieldnames(gui.allData);
mask        = false(1,length(gui.allData));
for i=1:length(sessionList)
    mask = mask|(~cellfun(@isempty,{gui.allData.(sessionList{i})}));
end
mice = find(mask); %find all the mice that contain data
for m = mice
    for s = sessionList'
        for tr = 1:length(gui.allData(m).(s{:}))
            gui.allData(m).(s{:})(tr).annot.(gui.annot.activeCh).(newStr) = [];
        end
    end
end
