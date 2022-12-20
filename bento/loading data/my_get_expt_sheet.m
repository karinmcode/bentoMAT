function [raw, header_row,header_fields,all_sessions,params] = my_get_expt_sheet(gui)
% [raw, param_row,param_fields] = my_get_expt_sheet(gui)
[~,~,raw]   = xlsread(gui.expt_sheet,'Sheet1');

%% get data
header_row = raw(2,:);
i4keep =~cell2mat(cellfun(@(x) numel(x)==0 | all(isnan(x)),header_row,'UniformOutput',false));
header_row = header_row(i4keep);
header_fields = header_row;
ncol = numel(header_fields);
for ihe = 1:ncol
    h = header_fields{ihe};
    h  = replace(h,{' ' '\'},'_');
    h  = replace(h,num2cell(':;=()'),'');
    header_fields{ihe} = h;
end
all_sessions = cell2struct(raw(3:end,1:ncol)',header_fields);


%% get params
params_row = raw(1,:);
params_names = params_row(2:2:end);
params_values = params_row(3:2:end);

for ipa = 1:numel(params_names)

    pa_raw = params_names{ipa};
    if isnan(pa_raw)
        continue
    end
    pa  = replace(pa_raw,{' ' '\' '/'},'_');
    pa  = replace(pa,num2cell(':;=()'),'');

    params.(pa)=params_values{ipa};
end

%disp('got params')

