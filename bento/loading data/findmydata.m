function info=findmydata(varargin)
%% info=findmydata(FilenameOrURL,[MYDATAFOLDER])
%% info=findmydata(gui)

% this function could be user specific.
% This function explains how the data is organized and outputs useful paths and filenames in the structure info as follows:
%   info.a : animal ID 'm543'
%   info.d : date 'YYMMDD'
%   info.rec : video recording number or trial number '001'
%   info.fn : filenames with extensions
%   info.fn0 : filenames without extensions
%   info.fo : folders
%   info.url : urls for various type of variables
%% SUBFIELDS more or less correspond to panels
%   .annot : annotation related
%   .features :
%   .phys :
%   .vid :
%   .tracking :
%   .events  :
%   .bento_proc : processing data from plugins

%% INPUT:
% The input can be a filename or an URL.
%%  DATA STRUCTURE : 
% info.fo = fullfile(MYDATAFOLDER , a , d , rec );
% info.fn = sprintf('m%s_%s_%s_CAM%s.ext,a,d,rec,camID)% for example m345_221215_001_CAM1.mat
% ext : file extensions could be .avi .mat .mp4 ...

%% Checking inputs

if numel(varargin)==2
    FilenameOrURL = varargin{1};
    MYDATAFOLDER = varargin{2};
elseif numel(varargin)==1
    if ischar( varargin{1})
        FilenameOrURL = varargin{1};
        MYDATAFOLDER = [];
    else
        gui = varargin{1};
        FilenameOrURL = mygcv(gui);
        MYDATAFOLDER = gui.mydatafolder;
    end
end

% find user defined data folder
if isempty(MYDATAFOLDER)
    [MYDATAFOLDER,fn0,ext] = fileparts(FilenameOrURL);
else
    [~,fn0,ext] = fileparts(FilenameOrURL);
end
if isempty(MYDATAFOLDER)
    MYDATAFOLDER = evalin('base','gui.mydatafolder;');
end

%% extract animal nb , data, and recording number of current experiment
temp = strsplit(fn0,'_');
if startsWith(fn0,'CAM')
    camID = temp{1}(4:end);
    a = temp{2};
    d = temp{3};
    rec = temp{4};
else
    a = temp{1};
    d = temp{2};
    rec = temp{3};
    camID = '1';%?

end
info.a =a;% 'm543'
info.d =d;%
info.rec =rec;
info.camID = camID;
anum = str2double(a(2:end));
dnum = str2double(d);
recnum = str2double(rec);
camIDnum = str2double(camID);

info.currentExpt = [anum dnum recnum];

%% define folders
EXPERIMENT_FOLDER = fullfile(MYDATAFOLDER,a,d,rec);
info.fo.annot = fullfile(EXPERIMENT_FOLDER,'annotations');
info.fo.features = fullfile(EXPERIMENT_FOLDER,'features');
info.fo.proc = fullfile(EXPERIMENT_FOLDER,'bento_proc');% processed data
info.fo.phys = EXPERIMENT_FOLDER;
info.fo.vid = EXPERIMENT_FOLDER;%winopen(info.fo.vid)
info.fo.tracking =EXPERIMENT_FOLDER;
info.fo.events = EXPERIMENT_FOLDER;

%% define extensions
info.ext.annot = '.annot';
info.ext.features = '__features.mat';
info.ext.phys = '.mat';
info.ext.vid = '*.mp4';
info.ext.tracking = '_tracking.mat';
info.ext.events = '.mat';


%% define filename without extension 
fn0 = sprintf('CAM%s_%s_%s_%s',camID,a,d,rec);
%fn0 =myfilename(mygcv(gui));

% if wildcard is used
wildcard_url = fullfile(info.fo.vid,[ fn0 info.ext.vid]);
if contains(wildcard_url,'*')
    fn0=my_find_wildcard_filename(wildcard_url);
    info.ext.vid = '.mp4';
end


info.fn0.annot = fn0;
info.fn0.features = fn0;% example: CAM1_m943_220413_003_motion_speed__features.mat
info.fn0.phys = sprintf('fluo_CAM%s',camID);
info.fn0.vid = fn0;
info.fn0.tracking = fn0;
info.fn0.events = sprintf('events_CAM%s',camID);



%% filename with extension 
info.fn.annot = sprintf('%s%s',info.fn0.annot,info.ext.annot);
info.fn.features = sprintf('%s%s',info.fn0.features,info.ext.features);
info.fn.phys = sprintf('%s%s',info.fn0.phys,info.ext.phys);
info.fn.vid = sprintf('%s%s',info.fn0.vid,info.ext.vid);
info.fn.tracking = sprintf('%s%s',info.fn0.tracking,info.ext.tracking);
info.fn.events = sprintf('%s%s',info.fn0.events,info.ext.events);



%% urls
info.url.annot = fullfile(info.fo.annot,info.fn.annot);
info.url.features = fullfile(info.fo.features,info.fn.features);
info.url.phys = fullfile(info.fo.phys,info.fn.phys);
info.url.vid = fullfile(info.fo.vid,info.fn.vid);
info.url.tracking =fullfile(info.fo.tracking,info.fn.tracking);
info.url.events =fullfile(info.fo.events,info.fn.events);

%% find the correct filename if the filename contains wildcard *
FIELDS = fieldnames(info.url);
nfi = numel(FIELDS);
for  ifi = 1:nfi
    fi = FIELDS{ifi};
    u = info.url.(fi);
    [~,oldfn,ext]=fileparts(u);
    oldfn = [oldfn ext];%#ok
    if contains(u,'*')
        files = dir(string(info.url.(fi)));%dir('Y:\Users\Karin\data\processed\aligned2vid\m943\220413\003\CAM1_mm943_220413_003*.mp4')
        newfn = files.name;
        info.fn.(fi) = replace(info.fn.(fi),oldfn,newfn);
        info.url.(fi)  = replace(u,oldfn,newfn);
    end

end


%% numbers
info.num.a = anum;
info.num.d = dnum;
info.num.rec = recnum;
info.num.camID = camIDnum;




