function filename = saveAnnotSheetTxt(movieNames,trial,suggestedName,promptOverride,saveAsTime)
%
% (C) Ann Kennedy, 2019
% California Institute of Technology
% Licensing: https://github.com/annkennedy/bento/blob/master/LICENSE.txt


DoesNotExistYet = (~isempty(trial.io.annot.fid)) && ~all(trial.io.annot.fid{1}==0);
if DoesNotExistYet
    fid     = trial.io.annot.fid{:};
else
    fid = [];
end
tmin    = trial.io.annot.tmin;
tmax    = trial.io.annot.tmax;
FR      = trial.io.annot.FR;

% should prompt channels to save (if it's not quicksave)? also, should save
% channels to the appropriate source file

%always prompt the save path?
if isempty(fid)
    DO_CREATE_FILE = true;
else
    if endsWith(fid,{'blank.annot' '.annot'})
        DO_CREATE_FILE = true;
    else
        DO_CREATE_FILE = false;
    end
end

if DO_CREATE_FILE %need to create a new file
    if(isempty(suggestedName))
        suggestedName = [pwd filesep 'annotations'];
    end
    if ~promptOverride
        [fname,pth] = uiputfile(suggestedName);
        if fname==0
            return
        end
    else
        [pth,fname,ext] = fileparts(suggestedName);
    end
    filename = fullfile(pth ,[ fname ext]);
elseif(exist('promptOverride','var') && ~promptOverride)
    [fname,pth] = uiputfile(fid);
    filename = fullfile(pth , fname);
else
    filename = fid;
end
if filename(1)==0
    return;
end
% set a default value of frameFlag
if(~exist('saveAsTime','var'))
    saveAsTime = false;
end


%check to see if the framerate was changed for display, and change back if so
if(FR~=trial.annoFR)
    for ch = fieldnames(trial.annot)'
        for beh = fieldnames(trial.annot.(ch{:}))'
            trial.annot.(ch{:}).(beh{:}) = round(trial.annot.(ch{:}).(beh{:}) * FR/trial.annoFR);
        end
    end
end

annot = trial.annot;
stim  = trial.stim;

if ~endsWith(filename,'.annot')
    filename = [filename '.annot'];
end
saveAnnot(filename,annot,tmin,tmax,FR,movieNames,stim,saveAsTime);



