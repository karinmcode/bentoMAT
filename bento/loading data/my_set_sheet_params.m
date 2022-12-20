function my_set_expt_params(gui,new_params)
%  my_set_expt_params(gui,new_params)
[params,params_name_column,raw] = my_get_expt_params(gui)
[~,~,raw]   = xlsread(gui.expt_sheet,'Sheet1');

params_row = raw(1,:);
params_names = params_row(2:2:end);
params_values = params_row(3:2:end);

% format names:



% i4empty = isempty(params_names) | isempty(params_values);
% params_names=params_names(~i4empty);
% params_values=params_values(~i4empty);

pnew_params=fieldnames(new_params);
npa = numel(pnew_params)
for ipa = 1:npa
    this_new_param = pnew_params{ipa};

    param_col = params_name_column.(this_new_param);
    raw{1,param_col+1}=new_params.(this_new_param);

end
xlswrite(gui.expt_sheet,raw,'Sheet1');

%disp('got params')



