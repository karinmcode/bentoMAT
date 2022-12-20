function setAnnotOverwriteMode(source)
%
% (C) Karin Morandell


gui=guidata(source);
answer = questdlg('Set annotation overwrite modes', ...
	'Annotation overwrite modes', ...
	'Overwrite=On','Overwrite=Off','Cancel','Overwrite=On');

switch answer
    case 'Overwrite=On'
        gui.annot.OverwriteMode = true;
    case 'Overwrite=Off'
        gui.annot.OverwriteMode = false;
    otherwise
        return;
end
guidata(gui.h0,gui);