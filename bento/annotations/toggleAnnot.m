function gui = toggleAnnot(gui,toggle,key,lastKey)
%
% (C) Ann Kennedy, 2019
% California Institute of Technology
% Licensing: https://github.com/annkennedy/bento/blob/master/LICENSE.txt


    gui.annot.prev = gui.annot.bhv;
    switch toggle
        case 'start'
            
            toggleAnnotBoxes(gui,'on',key);
            gui.annot.highlightStart = max(floor(...
                (gui.ctrl.slider.Value - gui.ctrl.slider.Min + 1/gui.data.annoFR)*gui.data.annoFR)-1,1);

        case 'cancel'
            toggleAnnotBoxes(gui,'off');
            gui.annot.highlightStart = [];
            gui.annot.highlighting = 0;
            
        case 'switch'
            toggleAnnotBoxes(gui,'on',key);
            inds    = getAnnotInds(gui);
            str     = gui.annot.hotkeys.(lastKey);
            gui.annot.bhv.(str)(inds) = true;
            
            gui.annot.highlightStart = max(floor(...
                (gui.ctrl.slider.Value - gui.ctrl.slider.Min + 1/gui.data.annoFR)*gui.data.annoFR)-1,1);
            
        case 'end'
            toggleAnnotBoxes(gui,'off');
            inds    = getAnnotInds(gui);
            if(any(strcmpi(lastKey,{'delete','backspace'}))) %erase everything
                for str = fieldnames(gui.annot.bhv)'
                    gui.annot.bhv.(str{:})(inds) = false;
                end
            elseif(isfield(gui.annot.hotkeys,lastKey))
                str = gui.annot.hotkeys.(lastKey);
                gui.annot.bhv.(str)(inds) = ~any(strcmpi(key,{'delete','backspace'}));
            end
            
            gui.annot.modified = 1;
            gui.annot.highlightStart = [];
            gui.annot.highlighting = 0;
            
        case 'fast'% replace label
            [inds,bhv] = getFastEditInds(gui); %figure out which bout/bhvr is just starting
            if ~isempty(bhv)
                gui.annot.bhv.(bhv)(inds) = 0;
                if(~any(strcmpi(key,{'delete','backspace'})))
                    str     = gui.annot.hotkeys.(key);
                    gui.annot.bhv.(str)(inds) = 1;
                end
                gui.annot.activeBeh = bhv;
                
                % go to next bout
                dummyEvent.Key = 'downarrow';
                evalKeyUp(gui.h0,dummyEvent);
            end
            gui.annot.modified = 1;

        case 'startOverwrite'%KM

            toggleAnnotBoxes(gui,'on',key);
            gui.annot.highlightStart = max(floor(...
                (gui.ctrl.slider.Value - gui.ctrl.slider.Min + 1/gui.data.annoFR)*gui.data.annoFR)-1,1);


        case 'endOverwrite'%KM

            toggleAnnotBoxes(gui,'off');
            inds    = getAnnotInds(gui);

            if(any(strcmpi(lastKey,{'delete','backspace'}))) %erase everything
                for str = fieldnames(gui.annot.bhv)'
                    gui.annot.bhv.(str{:})(inds) = false;
                end
            elseif(isfield(gui.annot.hotkeys,lastKey))
                str = gui.annot.hotkeys.(lastKey);
                gui.annot.bhv.(str)(inds) = ~any(strcmpi(key,{'delete','backspace'}));

                % put zeros for other labels% KM
                all_labels = fieldnames(gui.annot.bhv);
                labels2zero = setdiff(all_labels,str);
                for ila = 1:numel(labels2zero)
                    thislabel = labels2zero{ila};
                    gui.annot.bhv.(thislabel)(inds) = false;
                end

            end
            
            gui.annot.modified = 1;
            gui.annot.highlightStart = [];
            gui.annot.highlighting = 0;

    end
    
end

function toggleAnnotBoxes(gui, toggle, key)
    if(strcmpi(toggle,'on'))
        for f = fieldnames(gui.annot.Box)'
            gui.annot.Box.(f{:}).Visible = 'on';
            uistack(gui.annot.Box.(f{:}),'up');
            set(gui.annot.Box.(f{:}),'Ydata',kron(get(gui.(f{:}).axes,'ylim'), [1 1]));
            set(gui.annot.Box.(f{:}),'Xdata',[0 0 0 0]);

            %set box color
            if(any(strcmpi(key,{'delete','backspace'})))
                set(gui.annot.Box.(f{:}),'FaceColor',[.8 .8 .8]);
            else
                str = gui.annot.hotkeys.(key);
                set(gui.annot.Box.(f{:}),'FaceColor',gui.annot.cmap.(str)/3+2/3);
            end
        end
        
    else
        for f = fieldnames(gui.annot.Box)'
            gui.annot.Box.(f{:}).Visible = 'off';
            uistack(gui.annot.Box.(f{:}),'bottom');
            set(gui.annot.Box.(f{:}),'Ydata',[0 0 0 0]);
            set(gui.annot.Box.(f{:}),'FaceColor',[.8 .8 .8])
        end
    end
end

function inds = getAnnotInds(gui)
    p1 = gui.annot.highlightStart;
    p2 = ceil((gui.ctrl.slider.Value - gui.ctrl.slider.Min)*gui.data.annoFR);% KM ceil?
    if(p1<p2)
        inds = p1:p2;
    else
        inds = p2:p1;
    end
end

function [inds,bhv] = getFastEditInds(gui)
    t = floor((gui.ctrl.slider.Value - gui.ctrl.slider.Min)*gui.data.annoFR);% KM floor
    p1=t;p2=t;
    dtThresh = -inf;
    bhv = [];
    for b = fieldnames(gui.annot.bhv)'
        if(any(gui.annot.bhv.(b{:})(t:t+1)))
            dt = [1 gui.annot.bhv.(b{:})(2:(t+1))-gui.annot.bhv.(b{:})(1:t)];
            p1 = find(dt==1,1,'last');
            if(p1>dtThresh)
                p2 = t - 1 + find(gui.annot.bhv.(b{:})(t+1:end)==0,1,'first');
                bhv = b{:};
                dtThresh = p1;
            end
        end
    end
    inds = p1:p2;
end







