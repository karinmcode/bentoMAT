function gui = applyToAllMice(gui,action,varargin)
%
% (C) Ann Kennedy, 2019
% California Institute of Technology
% Licensing: https://github.com/annkennedy/bento/blob/master/LICENSE.txt



allData    = gui.allData;
inds    = gui.allPopulated;

switch action
    case 'add'
        newStr = varargin{1};
    case 'delete'
        killStr = varargin{1};
    case 'merge'
        oldList = varargin{1};
        newName = varargin{2};
        toKill  = varargin{3};
    case 'rename'
        oldName = varargin{1};
        newName = varargin{2};
end

for i = 1:size(inds,1)
    m       = inds(i,1);
    sess    = ['session' num2str(inds(i,2))];
    trial   = inds(i,3);
    anno    = allData(m).(sess)(trial).annot;
    
    channels = fieldnames(anno);
    for ch = 1:length(channels)
        switch action
            case 'add'
                anno.(channels{ch}).(newStr) = [];
                
            case 'delete'
                for b = 1:length(killStr)
                    anno.(channels{ch}) = rmfield(anno.(channels{ch}),killStr{b});
                end
                
            case 'merge'


                if(~isfield(anno.(channels{ch}),newName))
                    anno.(channels{ch}).(newName) = [];
                end
                newRast = []; bhvList = [newName;oldList];
                for b = 1:length(bhvList)
                    newRast = [newRast; anno.(channels{ch}).(bhvList{b})];
                    if(toKill && b>1)
                        anno.(channels{ch}) = rmfield(anno.(channels{ch}),bhvList{b});
                    end
                end
                newRast = cleanMergedRaster(newRast);
                anno.(channels{ch}).(newName) = newRast;
                gui.data.annot = anno;%KM
                if(toKill && b>1)
                    oldList2 = intersect(oldList,fieldnames(gui.annot.bhv));%KM
                    gui.annot.bhv = rmfield( gui.annot.bhv,oldList2);
                    
                end
                % update gui.annot.bhv
                gui=my_update_gui_annot_bhv(gui);

                updateLegend(gui,1)
            case 'rename'
                anno.(channels{ch}).(newName) = anno.(channels{ch}).(oldName);
                anno.(channels{ch}) = rmfield(anno.(channels{ch}),oldName);
        end
    end
    allData(m).(sess)(trial).annot = anno;
end

gui.allData = allData;










