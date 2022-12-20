function gui=updateAllData(gui)
%%  gui=updateAllData(gui)

data = gui.data;
m = gui.data.info.mouse;
session = gui.data.info.session;
trial = gui.data.info.trial;
try
    gui.allData(m).(session)(trial) = data;
catch err

    if strcmp(err.message,'Subscripted assignment between dissimilar structures.')
        oldFields = fieldnames(gui.allData(m).(session)(trial));
        newFields = fieldnames(data);
        diffFields = setdiff(oldFields,newFields);
        diffFields = setdiff(newFields,oldFields);


        for fi = oldFields(:)'
            gui.allData(m).(session)(trial).(fi{:})= data.(fi{:});
        end
    end
end

guidata(gui.h0,gui);

