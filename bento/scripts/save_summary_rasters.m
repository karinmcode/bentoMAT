function save_summary_rasters(pth,filter,suffix,FR,includeMARSOutput,overridePrompt)
%
% (C) Ann Kennedy, 2019
% California Institute of Technology
% Licensing: https://github.com/annkennedy/bento/blob/master/LICENSE.txt



[cmapDef,hotkeys] = loadPreferredCmap();
doPDF = exist('export_fig.m','file');
if(~exist('suffix','var')|isempty(suffix))
    suffix='';
elseif(suffix(1)~='_')
    suffix = ['_' suffix];
end

% set a default framerate if none specified
if(~exist('FR','var'))
    FR = 30;
elseif(isempty(FR))
    FR = 30;
end

% find files matching any of the provided filters
if(~iscell(filter))
    filter = {filter};
end
if(ischar(pth))
    if(pth(end)~=filesep)
        pth = [pth filesep];
    end
    disp(['Locating annotation files [' pth '**' filesep strjoin(filter,';') ']...'])
    files = rdir([pth '**' filesep filter{1}]);
    for i=2:length(filter)
        files = [files; rdir([pth '**' filesep filter{i}])];
    end
else
    files = pth; %secret mode: pass struct containing output of dir instead of just directory name
end
[~,i] = sort({files.name});
files = files(i);

% remove raw MARS output files (if not overridden)
if(~exist('includeMARSOutput','var') || ~includeMARSOutput)
    drop = ~cellfun(@isempty,strfind({files.folder},'output_v1_')) ...
            & ~cellfun(@isempty,strfind({files.name},'actions_pred')) ...
            & strcmpi(cellfun(@(x) x(end-2:end),{files.name},'uniformoutput',false),'txt');
    files(drop) = [];
end

% prompt to continue (if not overriden)
if(~exist('overridePrompt','var') || ~overridePrompt)
    ans = questdlg(['Found ' num2str(length(files)) ' files. Continue with summary?'],'Done!',...
                    'Yes','No','Yes');
    if(strcmpi(ans,'No'))
        return;
    end
else
    disp(['Found ' num2str(length(files)) ' files. Saving plots...']);
end

% remove previously generated pdfs
if(ls([pth 'all_mice_summary' suffix '.pdf']))
    delete([pth 'all_mice_summary' suffix '.pdf']);
end

disp('Please do not resize/close figure while summary images are being generated!')
% make/save plots!
img = [];
h   = figure(2);clf;
for f = 1:length(files)
%     disp([num2str(f) '/' num2str(length(files))]);
    clf;
    fname           = files(f).name;
    [~,name,~] = fileparts(fname);
    disp(name)
    [annot,maxTime] = loadAnyAnnot(fname);

    if(~isempty(annot))
        make_behavior_raster_summary(annot,cmapDef,hotkeys,maxTime,FR,h,strrep(fname,pth,''));
        im = getframe(h);
        img = cat(1,img,im.cdata);
        if(doPDF)
            export_fig([pth 'all_mice_summary' suffix '.pdf'],'-pdf','-append','-painters');
        end
    end
end
imwrite(img,[pth 'all_mice_summary' suffix '.png']);
close(h);

