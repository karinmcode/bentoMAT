function my_set_current_expt_params(gui,new_col_value)
%%  my_set_expt_params(gui,new_col_value)
% new_col_value: structure with fields corresponding to headers (if empty, displays possible fields)
% 
[raw, header_row,header_fields,data,params] = my_get_expt_sheet(gui);
if exist('new_col_value','var')==0
    new_col_value=[];
end
if isempty(new_col_value)
    disp('                                    ')
    disp('====================================')
    disp('POSSIBLE FIELDS FOR INPUT STRUCTURE:')
    disp([header_row;header_fields])
    disp('====================================')
    return;
end
thisMouse = str2double(gui.ctrl.expt.mouse.String{gui.ctrl.expt.mouse.Value});
thisSession = str2double(gui.ctrl.expt.session.String{gui.ctrl.expt.session.Value});
thisTrial = str2double(gui.ctrl.expt.trial.String{gui.ctrl.expt.trial.Value});

currentExpt = [thisMouse thisSession thisTrial];
ind4currentexpt = find(ismember(gui.allPopulated,currentExpt,'rows'));
irow = ind4currentexpt+2;
% i4empty = isempty(params_names) | isempty(params_values);
% params_names=params_names(~i4empty);
% params_values=params_values(~i4empty);

new_param_names=fieldnames(new_col_value);
npa = numel(new_param_names);

for ipa = 1:npa
    this_new_param = new_param_names{ipa};
    param_value = new_col_value.(this_new_param);

    i4col = find(strcmp(header_fields,this_new_param));
    if isempty(i4col)
        keyboard
    end
    raw{irow,i4col}=param_value;

end
xlswrite(gui.expt_sheet,raw,'Sheet1');

%disp('got params')



