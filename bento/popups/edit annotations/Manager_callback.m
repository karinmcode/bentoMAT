function Manager_callback(source,~,action)
%
% (C) Ann Kennedy, 2019
% California Institute of Technology
% Licensing: https://github.com/annkennedy/bento/blob/master/LICENSE.txt


h   = guidata(source);
errstr=[];
MergeXchannels =0;

if(isempty(h.bhv.Value) && any(strcmpi({'rename','delete'},action)))
    errstr = ['Please select an annotation to ' action];
elseif(length(h.bhv.Value)>1 && strcmpi(action,'rename'))
    errstr = 'More than one annotation selected';
elseif(length(h.bhv.Value)<=1 && strcmpi(action,'merge'))
    queststr = 'You selected only one label. Do you want to merge across channels?';

    ANS = questdlg(queststr );
    if ANS(1)=='Y'
        MergeXchannels =1;
    else
        return;
    end
end
if(~isempty(errstr))
    errordlg(errstr,'Error');
    return;
end
gui = h.gui;

switch action
    case 'play'% KM
        inds    = h.bhv.Value;
        selectBhv  = h.bhv.String(inds);
        if isempty(selectBhv)
            fprintf('No behavior was selected.')
            return;
        end
        gui.ctrl.slider.playSelectBhv(gui,selectBhv);

        guidata(source,h);
        guidata(gui.h0,gui);
        return;
    case 'add'
        answer          = inputdlg('Annotation to be added:','New annotation');
        if(isempty(answer)) return; end
        [answer,err]    = checkValid(answer,h.bhvMasterList);
        if(err) return; end
        h.bhvMasterList = [h.bhvMasterList; answer];
        gui = applyToAllMice(gui,'add',answer{:});

    case 'delete'
        inds    = h.bhv.Value;
        toKill  = h.bhv.String(inds);
        answer      = questdlg(['Confirm deletion of:  ' strjoin(toKill,', ')],'Confirm','Yes','No','Yes');
        if(strcmpi(answer,'no')) return; end
        h.bhvMasterList(inds) = [];
        gui = applyToAllMice(gui,'delete',toKill);

    case 'merge'
        if MergeXchannels
            gui=mymergelabelsXchannels(gui,h);
        else
            currChannel = gui.annot.activeCh;
            inds    = h.bhv.Value;
            toMerge = h.bhv.String(inds);
            potentialNewName = toMerge(~contains(toMerge,'prelabel'));
            if isempty(potentialNewName)
                potentialNewName = {''};
            end
            answer  = inputdlg({'Enter name of new (merged) annotation','Keep old annotations? (y/n)'},...
                ['Merging behaviors ' strjoin(toMerge,' ')],1,{potentialNewName{1},'no'});
            if(isempty(answer)) return; end
            if(~any(strcmpi({'y','yes','n','no'},answer{2})))
                answer{2} = questdlg('Keep old annotations?', '???', 'Yes','No','No');
            end
            [answer{1},err] = checkValid(answer{1},setdiff(h.bhvMasterList,answer{1}));
            if(err) return; end

            inds    = inds(~strcmpi(toMerge,answer{1})); %don't delete the behavior we merge into
            toMerge = h.bhv.String(inds);

            if(~any(strcmpi(h.bhvMasterList,answer{1})))
                h.bhvMasterList = [h.bhvMasterList; answer(1)];
            end
            toKill = 0;
            if(any(strcmpi({'n','no'},answer{2})))
                h.bhvMasterList(inds)=[];
                toKill = 1;
            end
            %gui = applyToAllMice(gui,'merge',toMerge,answer{1},toKill);% ORIGINAL CODE
            newName = answer{1};
            gui = mymergelabels(gui,toMerge,newName,toKill);
            %gui=guidata(gui.h0,gui);
            updatePlot(gui.h0,[]);%KM
            slider = gui.ctrl.slider;%KM
            slider.updateSliderAnnot(gui);%KM

            % update list of annotations on Annotation Manager
            NewListOfAnnotations = fieldnames(gui.annot.bhv);
            h.bhv.String = NewListOfAnnotations;
            h.bhv.Value = find(strcmp(NewListOfAnnotations,newName));
            myUnfoldAnnotStruct(gui);

            % go back to original channel
            gui.annot.activeCh = currChannel;
            gui.ctrl.annot.ch.Value = find(strcmp(gui.ctrl.annot.ch.String,currChannel));
            setChannel(gui.ctrl.annot.ch);
        end

        guidata(gui.h0,gui);
        return

        
    case 'split'
        [gui,h] = mySplitLabel(gui,h);
        h.bhv.Value     = [];
        h.bhv.String    = h.bhvMasterList;
        h.gui           = gui;

        guidata(h.fig,h);
        guidata(gui.h0,gui);
        return
    case 'rename'
        currChannel = gui.annot.activeCh;
        ind             = h.bhv.Value;
        old_name          = h.bhv.String{ind};
        new_name          = inputdlg('Rename annotation:','Rename',1,{old_name});
        new_name = new_name{1};
        if(isempty(new_name)) return; end
        [new_name,err]    = checkValid(new_name,h.bhvMasterList);
        if(err) return; end
        h.bhvMasterList{ind} = new_name;
        nch = numel(gui.annot.channels);

        for ich =1:nch
            thisChan = gui.annot.channels{ich};
            if isfield(gui.data.annot.(thisChan),old_name)
                Bouts = gui.data.annot.(thisChan).(old_name);
                gui.data.annot.(thisChan).(new_name)=Bouts;
                gui.data.annot.(thisChan)= rmfield(gui.data.annot.(thisChan),old_name);
            end
        end

        if isfield(gui.annot.bhv,old_name)
            rast=gui.annot.bhv.(old_name);
            gui.annot.bhv.(new_name) = rast;
            gui.annot.bhv= rmfield(gui.annot.bhv,old_name);
        end
        % get active m sess and trial
        m    = gui.data.info.mouse;
        sess = gui.data.info.session;
        trial   = gui.data.info.trial;

        gui.allData(m).(sess)(trial).annot = gui.data.annot;

        gui=my_update_gui_annot_bhv(gui);
        guidata(gui.h0,gui);% important because updateLegend calls guidata
        % update the legend
        updateLegend(gui,1)

        %gui = applyToAllMice(gui,'rename',toEdit{:},answer{:});% ORIGINAL CODE

        % UPDATE GUI
        updatePlot(gui.h0,[]);%KM
        slider = gui.ctrl.slider;%KM
        slider.updateSliderAnnot(gui);%KM


        % go back to original channel
        gui.annot.activeCh = currChannel;
        gui.ctrl.annot.ch.Value = find(strcmp(gui.ctrl.annot.ch.String,currChannel));
        setChannel(gui.ctrl.annot.ch);
    case 'summary'

        myUnfoldAnnotStruct(gui);
        return;
end

h.bhv.Value     = [];
h.bhv.String    = h.bhvMasterList;
h.gui           = gui;

guidata(source,h);
guidata(gui.h0,gui);

%% function 
    function [answer,err] = checkValid(answer,bhvList)
        answer(ismember(answer,'?!@#$%^&*()+=-<>,./\[]}{')) = [];
        answer = strrep(strtrim(answer),' ','_');

        err=0;
        if(isempty(answer))
            msgbox('Name must be at least one character long (special characters are removed.)');
            err=1;
        end
        if(any(strcmpi(bhvList,answer)))
            msgbox('A label with that name already exists.');
            err=1;
        end
    end

end

