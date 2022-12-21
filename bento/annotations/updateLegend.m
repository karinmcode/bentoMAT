function updateLegend(gui,showHidden)
% updateLegend(gui,showHidden)
% (C) Ann Kennedy, 2019
% California Institute of Technology
% Licensing: https://github.com/annkennedy/bento/blob/master/LICENSE.txt


    delete(gui.legend.axes.Children);
    if(~exist('showHidden'))
        showHidden = 0;
    end

    % get the list of behavior labels to show
    bhvList = fieldnames(gui.annot.bhv);
    keep = [];
    for b = 1:length(bhvList)
        keep(b) = gui.annot.show.(bhvList{b})==1;
    end
    hiddenList = bhvList(~keep);
    bhvList = bhvList(find(keep));

    hotkeys = struct();
    for f = fieldnames(gui.annot.hotkeys)'
        hotkeys.(gui.annot.hotkeys.(f{:})) = f{:};
    end

    % now display swatch, text, and hotkey for each
    bump = 0;
    if(~isempty(bhvList))
        text(.05,.025+bump,'Annotations','parent',gui.legend.axes);
        bump = bump+.05;
        for i = 1:length(bhvList)
            patch([.01 .125 .125 .01],[.01 .01 .05 .05]+bump,gui.annot.cmap.(bhvList{i}),...
                'parent',gui.legend.axes,'ButtonDownFcn',{@updateAnnot,gui.h0,bhvList{i}});

            text(.3, .025+bump,strrep(bhvList{i},'_',' '),'parent',gui.legend.axes,...
                'ButtonDownFcn',{@updateAnnot,gui.h0,bhvList{i}});
            if(isfield(hotkeys,bhvList{i}))
                text(.215, .025+bump,['[' hotkeys.(bhvList{i}) ']'],...
                        'horizontalalignment','center','parent',gui.legend.axes,...
                        'ButtonDownFcn',{@updateAnnot,gui.h0,bhvList{i}});
            end

            bump=bump+.05;
        end
        bump = bump+.05;
    end
    
    if(showHidden && ~isempty(hiddenList))
        text(.05,.025+bump,'Hidden annotations','parent',gui.legend.axes);
        bump = bump+.05;
            for i = 1:length(hiddenList)
                patch([.05 .165 .165 .05],[.01 .01 .05 .05]+bump,gui.annot.cmap.(hiddenList{i})/2+.5,...
                    'edgecolor','none','parent',gui.legend.axes,'ButtonDownFcn',{@updateAnnot,gui.h0,hiddenList{i}});
                text(.19, .025+bump,['\it ' strrep(hiddenList{i},'_',' ')],'color',[.25 .25 .25],...
                    'parent',gui.legend.axes,'ButtonDownFcn',{@updateAnnot,gui.h0,hiddenList{i}});
                bump=bump+.05;
            end
    end
    bump = bump + 0.05;
    ylim(gui.legend.axes,[-max(1-bump,0)/2 max(1-(1-bump)/2,bump)]);
end

function updateAnnot(~,~,h0,str)
    gui = guidata(h0);
    gui.annot.show.(str) = ~gui.annot.show.(str);
    guidata(gui.h0,gui);

    updateSliderAnnot(gui);
    updateLegend(gui,1);
    updatePlot(gui.h0,[]);
end