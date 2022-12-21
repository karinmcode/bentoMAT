function h = incSliderVal(source)
% this function just implements normal slider behavior (arrow and bar
% clicks).
%
% (C) Ann Kennedy, 2019
% California Institute of Technology
% Licensing: https://github.com/annkennedy/bento/blob/master/LICENSE.txt


gui = guidata(source);
h = gui.ctrl.slider;

% check if out-of-range
if(((h.Value + gui.Action(1))>h.Max)|((h.Value + gui.Action(1))<h.Min))
    return
end


% check for overshoot
p       = get(0,'PointerLocation');         %click location
marker  = getpixelposition(h.Marker,true);
marker  = gui.h0.Position(1) + marker(1) + marker(3)/2; %marker coordinates
try
if(gui.Action(2) & (gui.h0==gcf) && (sign(p(1) - marker) ~= sign(gui.Action(1))))
    return;
end
end

h.Value = h.Value + gui.Action(1);