function clearAction(source,~)
% KM heavily modified
% (C) Ann Kennedy, 2019
% California Institute of Technology
% Licensing: https://github.com/annkennedy/bento/blob/master/LICENSE.txt



gui = guidata(source);
if(any(gui.Action~=0) )


    gui = guidata(gui.h0);
    slider = gui.ctrl.slider;
    if(any(gui.Action))
        test = toc(slider.timer);
        if(test>0.15)
            if(~isempty(strfind(gui.Action,'drag')))
                gui = slider.setSliderVal(gui);%KM
            else
                gui = slider.incSliderVal(gui);%KM
            end
          
        end
    end
    drawnow;


end

gui.Action = 0;
guidata(source,gui);