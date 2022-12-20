function [params,raw,ind4currentexpt]=my_get_current_expt_params(gui)
%  my_get_expt_params(gui,new_params)
%edit my_set_current_expt_params.m
%edit my_get_expt_params.m
[~,~,raw]   = xlsread(gui.expt_sheet,'Sheet1');

params_row = raw(2,:);
i4keep = ~cell2mat(cellfun(@(x) (numel(x)==0 | all(isnan(x))),params_row,'uniformoutput',false));
params_row = params_row(i4keep);
params_names = params_row;

%find current experiments
thisMouse = str2double(gui.ctrl.expt.mouse.String{gui.ctrl.expt.mouse.Value});
thisSession = str2double(gui.ctrl.expt.session.String{gui.ctrl.expt.session.Value});
thisTrial = str2double(gui.ctrl.expt.trial.String{gui.ctrl.expt.trial.Value});

currentExpt = [thisMouse thisSession thisTrial];
ind4currentexpt = find(ismember(gui.allPopulated,currentExpt,'rows'));
irow = ind4currentexpt+2;
npa = numel(params_names);

params_values = raw(irow,1:npa);
params = struct();


for ipa = 1:npa

    this_param = params_names{ipa};
    pa  = replace(this_param,{' ' '\'},'_');
    pa  = replace(pa,num2cell(':;=()'),'');
    params.(pa) = params_values{ipa};

end




