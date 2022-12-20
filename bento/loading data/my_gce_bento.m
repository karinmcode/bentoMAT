function [a , d , rec,info] = my_gce_bento(gui)
%% [a , d , rec,info] = my_gce_bento(gui);


% get current experiment

ctrl = gui.ctrl;
expt = gui.ctrl.expt;
a = str2double(expt.mouse.String{expt.mouse.Value});
d = str2double(expt.session.String{expt.session.Value});
rec = str2double(expt.trial.String{expt.trial.Value});

info.a =expt.mouse.String{expt.mouse.Value};
info.d =expt.session.String{expt.session.Value};
info.rec =num2str(rec,'%03.f');

[params,raw,ind4currentexpt]=my_get_current_expt_params(gui);
currentExpt = [a d rec];
ind4currentexpt = find(ismember(gui.allPopulated,currentExpt,'rows'));
irow = ind4currentexpt+2;

[~,~,raw]   = xlsread(gui.expt_sheet,'Sheet1');

info.index = ind4currentexpt;
info.fo.proc = fullfile(gui.mydatafolder,['m' info.a],info.d,info.rec);
info.url.vid = fullfile(gui.mydatafolder,info.a,info.d,info.rec,replace(params.Behavior_movie,'.findme','.mp4'));

info.camID = mycamID(params.Behavior_movie);