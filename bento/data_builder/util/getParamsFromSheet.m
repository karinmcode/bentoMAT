function params = getParamsFromSheet(raw)

params_row = raw(1,:);
% i4nan = cell2mat(cellfun(@(x) isnan(x) , params_row,'UniformOutput',false));
params_row = params_row(1:19);
params_names = params_row(2:2:end);% first cell for bento data folder
params_values = params_row(3:2:end);

% format names:



% i4empty = isempty(params_names) | isempty(params_values);
% params_names=params_names(~i4empty);
% params_values=params_values(~i4empty);

for ipa = 1:numel(params_names)

    pa_raw = params_names{ipa};
    pa  = replace(pa_raw,{' ' '\'},'_');
    pa  = replace(pa,num2cell(':;=()'),'');

    params.(pa)=params_values{ipa};
end

%disp('got params')

if isfield(params,'My_Data_Folder')
    params.mydatafolder = params.My_Data_Folder;
end