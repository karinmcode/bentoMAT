function  gui=my_update_gui_annot_bhv(gui)
%  gui=my_update_gui_annot_bhv(gui)
data = gui.data;
gui = transferAnnot(gui,data);

slider = gui.ctrl.slider;
slider.updateSliderAnnot(gui);
%img = makeAllChannelBhvImage(gui,data,cmap,inds,tmax,showAnnot)
% img = makeBhvImage
% 
%myUnfoldAnnotStruct(gui)
%myUnfoldAnnotStruct(guidata(gui.h0))

