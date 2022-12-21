function Import_features_from_data_folder(source)
% Import_features_from_data_folder(gcf)
edit Import_features_from_data_folder
%% Get variables

% gui
try
    gui = guidata(source);
catch
    gui = source;
end

% folders
info = findmydata(gui);
url = info.url;

folder_proc = info.fo.vid;% where is the data folder with encoder, speed, position data
url_save = url.features;% url of mat file where bento features are saved mywinopen(url_save)
mymkdir(url_save);%mywinopen(url_save)


%% Select features
if exist(url_save,'file')==2%mywinopen(url_save)
    list_var= whos('-file', url_save);
    list_var = {list_var.name};
    list2import = horzcat(list_var,setdiff({'position' 'speed' 'other' 'keyboard'},list_var));
else
    list2import = {'position' 'speed' 'other' 'keyboard'};
end

ANSidx=listdlg('ListString',list2import,'InitialValue',[1 2]);
if isempty(ANSidx)
    disp 'did not select new feature'
    return
end

featureNames = list2import(ANSidx);


for ife = 1:numel(featureNames)%ife=2
    featureName = featureNames{ife};
    if ismember(featureName,{'other' })
        featureName = inputdlg('input name of variable to import');
        disp not coded yet
        keyboard
    end

    if ismember(featureName,{'keyboard'})

        disp not coded yet
        keyboard
    end

    switch featureName
    %% position case
        case 'position'

        this_url =  fullfile(folder_proc,sprintf('encoder_CAM%g.mat',info.num.camID));
        temp=load(this_url,'position');%mywinopen(url_position)
        position=temp.position.values_vidproc;
        if exist(url_save,'file')==0
            save(url_save,'position');%load(url_save)winopen(folder)
        else
            save(url_save,'position','-append');%load(url_save)winopen(folder)
        end
        newFeatures.position = position;
    

    %% speed case
        case 'speed'

        this_url =  fullfile(folder_proc,sprintf('encoder_CAM%g.mat',info.num.camID));
        temp=load(this_url,featureName);%mywinopen(url_position)
        speed=temp.(featureName).values_vidproc;
        if exist(url_save,'file')==0
            save(url_save,'speed');%load(url_save)winopen(folder)
        else
            save(url_save,'speed','-append');%load(url_save)winopen(folder)
        end
        newFeatures.speed = speed;


        otherwise
            load(url_save,featureName)
            newFeatures.(featureName) = eval([featureName ';']);

    end

end

%% check numel of features
 [newFeatures,Nmax]=myCheckFeaturesDataValidity(newFeatures);
% update cell in excel/google sheet with filename
%ANS = questdlg('do you want to update experiment sheet with feature file location ?');
ANS = 'Y';% provisory
if ANS(1)=='Y'
    params.Tracking = replace(url_save,gui.pth,'');
    my_set_current_expt_params(gui,params);

end

%gui = myUploadFeaturesToData(gui,data_feat); future function 

% add features to current data
gui=myImportNewFeatures(gui,newFeatures);

guidata(gui.h0,gui)
updatePlot(gui.h0,[])


end




