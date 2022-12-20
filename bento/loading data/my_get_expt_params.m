function [params,params_name_column,raw] = my_get_expt_params(gui)
% params = my_get_expt_params(gui)


[~,~,raw]   = xlsread(gui.expt_sheet,'Sheet1');

params_row = raw(1,:);
params_names = params_row(2:2:end);
params_values = params_row(3:2:end);

% format names:



% i4empty = isempty(params_names) | isempty(params_values);
% params_names=params_names(~i4empty);
% params_values=params_values(~i4empty);
param_col = 0;
for ipa = 1:numel(params_names)
    param_col = param_col+2;
    pa_raw = params_names{ipa};
    pa  = replace(pa_raw,{' ' '\'},'_');
    pa  = replace(pa,num2cell(':;=()'),'');

    params.(pa)=params_values{ipa};
    params_name_column.(pa) = param_col
end

%disp('got params')



