function toggleEnabled(source)
%
% (C) Ann Kennedy, 2019
% California Institute of Technology
% Licensing: https://github.com/annkennedy/bento/blob/master/LICENSE.txt



gui = guidata(source);

h = dialog('Name','Toggle object visibility');
p = h.Position;
    
count = 0;
for f = fieldnames(gui.enabled)'
    if(strcmpi(f{:},'ctrl')|strcmpi(f{:},'welcome'))
        continue;
    end
    if(gui.enabled.(f{:})(1))
        count = count+1;
        uicontrol(h,'Style','checkbox','Value',gui.enabled.(f{:})(2),...
                'String',['  ' f{:}],'fontsize',12,...
                'Position',[75 count*30+10 100 30],'callback',@toggleEnabledCallback);
    end
end

% open dialog box where my mouse cursor is
p(3) = 250;
p(4) = 30*count+60;
p(1:2) = get(0,'PointerLocation');
%MonitorPositions = get(0,'MonitorPositions');
h.Position = p;

guidata(h,gui.h0);
end

function toggleEnabledCallback(source,~)

h0 = guidata(source);
gui = guidata(h0);
gui.enabled.(strtrim(source.String))(2) = source.Value;

if(gui.enabled.annot(2))
    gui.ctrl.annot.panel.Visible = 'on';
else
    gui.ctrl.annot.panel.Visible = 'off';
end
updateSliderAnnot(gui);
gui = redrawPanels(gui);
guidata(gui.h0,gui);
updatePlot(gui.h0,[]);
end