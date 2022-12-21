function [mouse,enabled,pth,hotkeys,params] = unpackExperiment(raw,skipvideos)
% parses the metadata in the excel sheet, then loads all of the listed
% files and formats their data.
%
% (C) Ann Kennedy, 2019
% California Institute of Technology
% Licensing: https://github.com/annkennedy/bento/blob/master/LICENSE.txt

% get rid of nans
for ise = 3:size(raw,1)
    inds = [4 5 9 10 15 16 17 19];
    inds(inds>size(raw,2))=[];
    mask = cellfun(@sum,cellfun(@isnan,raw(ise,inds),'uniformoutput',false));
    raw(ise,inds(find(mask))) = {''};
    inds = setdiff(1:17,inds);
    inds(inds>size(raw,2))=[];
    mask = cellfun(@sum,cellfun(@isnan,raw(ise,inds),'uniformoutput',false));
    raw(ise,inds(find(mask))) = {[]};
end

%
if(~exist('skipvideos','var'))
    skipvideos = 0;
end

%OS compatibility for file paths:
raw(cellfun(@isstr,raw)) = strrep(raw(cellfun(@isstr,raw)),'\',filesep);
raw(cellfun(@isstr,raw)) = strrep(raw(cellfun(@isstr,raw)),'/',filesep);
pth = raw{1,1};
if(~isempty(pth))
    if(pth(end)~=filesep)
        if(isstr(pth))
            pth = [pth filesep];
        else
            pth = [];
        end
    end
end

[data,match,~] = reconcileSheetFormats([],raw);
params = getParamsFromSheet(raw);
hotkeys = struct();

if(isnumeric(raw{1,3})&&~isnan(raw{1,3})) %if there's a common Ca framerate
    data(:,match.FR_Ca) = raw(1,3);
end
if(isnumeric(raw{1,5})&&~isnan(raw{1,5})) %if there's a common bhvr movie framerate
    data(:,match.FR_Anno) = raw(1,5);
end

%window visibility:
enabled.movie     = any(~cellfun(@isempty,data(:,match.Behavior_movie)))*[1 1];
enabled.annot     = [1 1];
enabled.legend    = [1 any(~cellfun(@isempty,data(:,match.Annotation_file)))];
enabled.traces    = any(~cellfun(@isempty,data(:,match.Calcium_imaging_file)))*[1 1];
enabled.tracker   = any(~cellfun(@isempty,data(:,match.Tracking)))*[1 1];
enabled.features  = any(~cellfun(@isempty,data(:,match.Tracking)))*[1 1];
enabled.audio     = any(~cellfun(@isempty,data(:,match.Audio_file)))*[1 1];

enabled.tsne      = [0 0];
enabled.scatter   = any(~cellfun(@isempty,data(:,match.Calcium_imaging_file)))*[1 0];
enabled.fineAnnot = any(~cellfun(@isempty,data(:,match.Annotation_file)))*[1 1];

%load the data:
mouse = struct();
prevCa = '';
rast   = [];
CaTime = [];
nse = size(data,1);
for ise=1:nse


    if(isempty(data{ise,match.Mouse}))
        continue;
    end


    d = cell2struct(data(ise,:)',fieldnames(match));%KM
    thisMouse = ['m' num2str(d.Mouse)];
    thisSession = num2str(d.Sessn);
    thisTrial = num2str(d.Trial,'%03.f');
    fprintf('\n> Uploaded experiment %g/%g: m%g_%s_%s',ise,nse,d.Mouse,thisSession,thisTrial);


    this_exp         = struct();
    this_exp.stim    = data{ise,match.Stim};
    this_exp.CaFR    = data{ise,match.FR_Ca};

    this_exp.annoFR  = data{ise,match.FR_Anno};


    offset          = data{ise,match.Offset};
    if(~isnumeric(offset))
        offset = str2num(offset);
    end
    hasOffset       = ~isempty(offset)&&~any(isnan(offset));

    %% load calcium traces----------------------------------------------------------
    if(~isempty(data{ise,match.Calcium_imaging_file}))
        Calcium_imaging_file = strip(strip(data{ise,match.Calcium_imaging_file},'left','.'),'left');
        fid     = fullfile(params.My_Data_Folder, Calcium_imaging_file);
        tstart  = data{ise,match.Start_Ca};
        tstop   = data{ise,match.Stop_Ca};
        CaFR    = data{ise,match.FR_Ca};

        if(~strcmpi(fid,prevCa)) %save some time by not re-loading files used previously
            [rast,CaTime,spikes,ROIs]    = unpackCaData(fid);
            if(size(rast,1)>size(rast,2)) %sometimes data is stored transposed. we usually have more timepoints than cells, so use this fact to detect this
                rast=rast';spikes=spikes';
            end
            prevCa  = fid;
            this_exp.rast    = rast;
            this_exp.spikes  = spikes;
            this_exp.ROIs    = ROIs;
        else
            this_exp.rast    = rast;
            this_exp.spikes  = spikes;
            this_exp.ROIs    = ROIs;
        end
        cutCa = ~any(isempty(tstart))&&~any(isnan(tstart));
        if(cutCa)
            if(~isnumeric(tstart))
                tstart = str2num(tstart); tstop = str2num(tstop);
            end
            try
                this_exp.rast    = rast(:,tstart:tstop);
            catch
                keyboard
            end
            if(~isempty(spikes))
                this_exp.spikes  = spikes(:,tstart:tstop);
            end
        end
        this_exp.sm_rast = this_exp.rast;
        drdt = [zeros(size(this_exp.rast(:,1),1),1) this_exp.rast(:,2:end)-this_exp.rast(:,1:end-1)];
        this_exp.ddt = smoothts(drdt,'g',50,10)*20;
        if(~isempty(CaTime))
            this_exp.CaTime   = CaTime;
            this_exp.CaFR     = 1/mean(CaTime(2:end)-CaTime(1:end-1)); % trust the timestamp over the user
        else
            this_exp.CaTime   = (1:size(this_exp.rast,2))/CaFR;
            this_exp.CaFR     = CaFR;
        end
        if(hasOffset)
            this_exp.CaTime  = this_exp.CaTime + offset;
        end

        % for 2d plots:
        this_exp.proj.d1 = [];
        this_exp.proj.d2 = [];
    else
        this_exp.CaTime=[];
        this_exp.rast=[];
        this_exp.proj.d1=[];
        this_exp.proj.d2=[];
    end

    % adds cross-day alignments if available. assumes for now that my data
    % format is being used.
    if(isfield(match,'Alignments')&&~isempty(data{ise,match.Alignments}))
        fid = [pth data{ise,match.Alignments}];
        sessn = ['session' num2str(data{ise,match.Sessn})];
        temp = load(fid);
        this_exp.rast_matched  = this_exp.rast(temp.alignedCells.(sessn),:);
        this_exp.match         = temp.alignedCells.(sessn);
        if(isfield(temp.bounds,['day' sessn(end)]))
            if(exist('ICs'))
                this_exp.units         = temp.ICs.(['day' sessn(end)]);
            else
                this_exp.units = [];
            end
            this_exp.bounds        = temp.bounds.(['day' sessn(end)]);
        elseif(isfield(temp.bounds,sessn))
            if(exist('ICs'))
                this_exp.units         = temp.ICs.(sessn);
            else
                this_exp.units = [];
            end
            this_exp.bounds        = temp.bounds.(sessn);
        else
            this_exp.units         = [];
        end
    else
        this_exp.rast_matched = [];
        this_exp.units = [];
    end

    %% ----------------------------------------------------------
    %% link movies----------------------------------------------------------
    if(~isempty(data{ise,match.Behavior_movie}))
        colList   = strsplit(data{ise,match.Behavior_movie},';;');

        for col = 1:length(colList)
            movieList = strsplit(colList{col},';');
            for j = 1:length(movieList)
                thisMovie = movieList{j};
                thisMovie_stripped = strip(strip(strtrim(thisMovie),'left','.'),'left',filesep);
                if contains(thisMovie,'.findme')
                    info = findmydata(thisMovie,params.My_Data_Folder);
                    this_exp.io.movie.fid{col,j} = ...
                        info.url.vid;
                else
                    this_exp.io.movie.fid{col,j} = ...
                        strtrim([pth strip(strip(strtrim(movieList{j}),'left','.'),'left',filesep)]);
                end

                %% SET FR for movie and annotations if commone FR has not been set% KM
                if isempty(this_exp.annoFR)% get it from the video
                    if myisfile(this_exp.io.movie.fid{col,j})
                        Reader = VideoReader(this_exp.io.movie.fid{col,j});%winopen()
                    else
                        keyboard
                    end
                    this_exp.annoFR = Reader.FrameRate;
                end

                switch(this_exp.io.movie.fid{col,j}(end-2:end))
                    case 'seq'
                        this_exp.io.movie.readertype{col,j} = 'seq';
                    otherwise
                        this_exp.io.movie.readertype{col,j} = 'vid';
                end

            end
        end


        % provisory
        this_exp.io.movie.FR  = this_exp.annoFR; %should change this to allow multiple FR's in the future

        if(~skipvideos)
            if(raw{1,9})
                this_exp.io.movie.tmin = data{ise,match.Start_Anno};
                this_exp.io.movie.tmax = data{ise,match.Stop_Anno};
            else
                this_exp.io.movie.tmin = 1;
                switch(this_exp.io.movie.fid{1}(end-2:end))
                    case 'seq'
                        tmax = inf;
                        for j = 1:length(this_exp.io.movie.fid)
                            info = seqIo(this_exp.io.movie.fid{j},'getInfo');
                            tmax = min([tmax info.numFrames]);
                        end
                    otherwise
                        tmax = inf;
                        for j = 1:length(this_exp.io.movie.fid)
                            try
                                info = VideoReader(this_exp.io.movie.fid{j});
                            catch
                                error(['I wasn''t able to find/load a video at: ' this_exp.io.movie.fid{j}]);
                                %mywinopen(strtemp.io.movie.fid{j})
                            end
                            timestamps = getVideoTimestamps(this_exp.io.movie.fid{j});
                            if timestamps % check for timestamp files, update movie data accordingly
                                tmax = min([tmax length(timestamps)]);
                                this_exp.io.movie.FR = 1/mean(timestamps(2:end)-timestamps(1:end-1));
                            else
                                tmax = min([tmax round(info.Duration*info.FrameRate)]);
                            end
                            if(~isempty(data{ise,match.Calcium_imaging_file}) && isempty(this_exp.CaTime)) % hack to check for accompanying Ca timestamps :[
                                timestamps = getVideoTimestamps(this_exp.io.movie.fid{j},'_Ca');
                                if timestamps
                                    this_exp.CaTime = timestamps(1:length(this_exp.rast))';
                                    this_exp.CaFR = 1/mean(timestamps(2:end)-timestamps(1:end-1));
                                end
                            end
                        end
                end
                this_exp.io.movie.tmax = tmax;
            end
        end

    else
        this_exp.io.movie = struct();
    end

    %% add tracking data----------------------------------------------------
    if(enabled.tracker(1))
        this_exp.trackTime = [];


        if(~isempty(data{ise,match.Tracking}))
            trackList = strsplit(data{ise,match.Tracking},';');


            for itrackFile = 1:length(trackList)
                cellContent = strip(strip(trackList{itrackFile},'left','.'),'left',filesep);
                 isAbsolutePath=contains(cellContent,':');
                if isAbsolutePath
                    fid = cellContent;
                else
                    fid = [pth cellContent];
                end
                [~,~,ext] = fileparts(fid);

                if(strcmpi(ext,'.mat'))

                    if myisfile(fid)
                        temp = load(fid); %virtual load would be faster/more memory friendly, but laggier

                        %                     f = fieldnames(temp);
                        %                     if(length(f)==2) %what is this for? i forget.
                        %                         temp=temp.(f{2});
                        %                     elseif(length(f)==1)
                        %                         temp = temp.(f{1});
                        %                     end
                        temp=myCheckFeaturesDataValidity(temp);%KM

                        this_exp.tracking.args{itrackFile} = temp;
                        this_exp.io.feat.fid{itrackFile} = fid;
                        this_exp.trackTime = timestamps;

                    else
                        warning('Did not find tracking/feature file: %s',fid);
                        mywinopen(fid);
                        keyboard
                    end

                elseif strcmpi(ext,'.findme')
                    if isempty(trackList{itrackFile})
                        info = findmydata(thisMovie,params.mydatafolder);
                    else
                        info = findmydata(trackList{itrackFile},params.mydatafolder);
                    end
                    if myisfile(info.url.features)
                        temp = load(info.url.features); %virtual load would be faster/more memory friendly, but laggier

                        %                     f = fieldnames(temp);
                        %                     if(length(f)==2) %what is this for? i forget.
                        %                         temp=temp.(f{2});
                        %                     elseif(length(f)==1)
                        %                         temp = temp.(f{1});
                        %                     end
                        temp=myCheckFeaturesDataValidity(temp);%KM

                        this_exp.tracking.args{itrackFile} = temp;
                        this_exp.io.feat.fid{itrackFile} = info.url.features;
                        this_exp.trackTime = timestamps;

                    else
                        warning('Did not find tracking/feature file: %s',fid);
                        mywinopen(fid);
                        keyboard
                    end


                elseif(strcmpi(ext,'.json'))
                    if(exist('jsondecode','builtin'))
                        disp('loading tracking data');
                        this_exp.tracking.args{itrackFile} = jsondecode(fileread(fid));
                    elseif(exist('loadjson','file'))
                        disp('Using loadjson.')
                        this_exp.tracking.args{itrackFile} = loadjson(fid);
                    else
                        disp('Please download jsonlab (https://github.com/fangq/jsonlab) or upgrade to Matlab 2016b or later.')
                        this_exp.tracking.args{itrackFile} = [];
                    end
                    this_exp.io.feat.fid{itrackFile} = fid;

                elseif(strcmpi(ext,'.h5')) %DeepLabCut or JAX output
                    trackType = promptTrackType({'JAX','DLC'},trackList{itrackFile});
                    this_exp.tracking.fun = trackType;
                    if contains(trackType,'DLC')
                        args = h5read(fid,'/df_with_missing/table');
                        this_exp.tracking.args{itrackFile} = args.values_block_0;
                    elseif contains(trackType,'JAX')
                        args = struct();
                        args.points = h5read(fid,'/poseest/points');
                        args.ids = h5read(fid,'/poseest/instance_track_id');
                        this_exp.tracking.args{itrackFile} = args;
                        this_exp.trackTime = (1:length(this_exp.tracking.args{1}.ids))/30;
                    end

                elseif(strcmpi(ext,'.csv')) %either DeepLabCut output or SimBA features
                    trackType = promptTrackType({'DLC','SimBA features'},trackList{itrackFile});
                    if contains(trackType,'DLC')
                        [numvals,textfields] = xlsread(fid);
                        args.data = numvals';
                        args.ids = string(textfields(:,2)) + '_' + string(textfields(:,3));
                    else
                        fh = fopen(fid);
                        headers = strsplit(fgetl(fh),',');
                        fclose(fh);
                        ncol = length(headers);
                        fh = fopen(fid);
                        C = textscan(fh,'%f','headerlines',1,'Delimiter',',');
                        fclose(fh);
                        C = reshape(C{1},ncol,[]);
                        args.data = C(2:end,:);
                        args.ids = char(headers(2:end));
                    end
                    this_exp.tracking.args{itrackFile} = args;
                else % I couldn't unpack the tracking data :[
                    this_exp.tracking.args={[]};
                end
            end

            % everything that follows is shameful hacks!
            % we need a better way to get timestamps for the tracking data...
            if(isfield(this_exp,'tracking'))

                if(isfield(this_exp.tracking.args{1},'keypoints') && isfield(this_exp.tracking.args{1},'fps'))
                    this_exp.trackTime = (1:length(this_exp.tracking.args{1}.keypoints))/double(this_exp.tracking.args{1}.fps);

                elseif(isfield(this_exp.tracking.args{1},'tMax')) %hacks for jellyfish
                    this_exp.trackTime = (1:this_exp.tracking.args{1}.tMax)/this_exp.CaFR;

                elseif isfield(this_exp.tracking.args{1},'fps')
                    if length(fieldnames(this_exp.tracking.args{1}))==2
                        datafield = setdiff(fieldnames(this_exp.tracking.args{1}),'fps');
                    elseif isfield(this_exp.tracking.args{1},'data')
                        datafield = 'data';
                    elseif isfield(this_exp.tracking.args{1},'data_smooth')
                        datafield = 'data_smooth';
                    else
                        f = fieldnames(this_exp.tracking.args{1});
                        ans = inputdlg('Which field holds the tracking data?');
                        if isfield(this_exp.tracking.args{1},ans{:})
                            datafield = ans{:};
                        else
                            datafield = [];
                        end
                    end
                    if datafield
                        this_exp.trackTime = (1:length(this_exp.tracking.args{1}.(datafield)))/double(this_exp.tracking.args{1}.fps);
                    end

                else
                    if(~isempty(data{ise,match.Behavior_movie}) && length(this_exp.io.movie.fid)==1)
                        if(~strcmpi(this_exp.io.movie.fid{1}(end-2:end),'seq'))
                            timestamps = getVideoTimestamps(this_exp.io.movie.fid{1});
                        else
                            temp         	= seqIo(data.io.movie.fid{1},'reader');
                            disp('getting timestamps...');
                            timestamps = getSeqTimestamps(data.io.movie.fid{1},temp);
                        end
                        if(timestamps)
                            this_exp.trackTime = timestamps;
                            continue;
                        end
                    end

                    % I give up, let's just ask the user
                    if 1%provisory
                        warning('KM provisory annoFR setting');
                        if ~isempty(this_exp.annoFR)
                            fps = this_exp.annoFR;
                        else
                            keyboard;
                        end
                    else
                        ans = inputdlg('What''s the framerate of the tracking data?','Frame rate',1,{'video'});
                        if strcmp(ans{1},'video')
                            fps = this_exp.annoFR;
                        else
                            fps = str2num(ans{:});
                        end
                    end
                    ARG1 = this_exp.tracking.args{1};
                    if isnumeric(ARG1)
                        this_exp.trackTime = (1:length(ARG1))/fps;
                    elseif length(fieldnames(ARG1))==1
                        f = fieldnames(ARG1);
                        this_exp.trackTime = (1:length(ARG1.(f{:})))/fps;
                    elseif length(fieldnames(ARG1))==2
                        if isfield(ARG1,'data')
                            this_exp.trackTime = (1:length(ARG1.data))/fps;
                        elseif isfield(ARG1,'data_smooth')
                            this_exp.trackTime = (1:length(ARG1.data_smooth))/fps;
                        elseif isfield(ARG1,'features')
                            this_exp.trackTime = (1:length(ARG1.features))/fps;
                        end
                    elseif length(fieldnames(ARG1))>2
                        f = fieldnames(ARG1);
                        this_exp.trackTime = (1:length(ARG1.(f{1})))/fps;
                    else
                        f = fieldnames(ARG1);
                        this_exp.trackTime = (1:length(ARG1.(f{:})))/fps;
                    end


                    if isempty(this_exp.trackTime)
                        f = fieldnames(ARG1);
                        this_exp.trackTime = (1:length(ARG1.(f{1})))/fps;
                        warning('KM: computed trackTime in features')
                    end

                end
            end
        else
            this_exp.tracking.args = [];
        end
        if isempty(this_exp.trackTime)
            keyboard
        end

    else
        this_exp.io.feat = [];
    end

    %% add audio data-------------------------------------------------------
    if(enabled.audio(1))
        if(~isempty(data{ise,match.Audio_file}))
            [~,~,ext] = fileparts(data{ise,match.Audio_file});
            strSpect  = strip(strip(strrep(data{ise,match.Audio_file},ext,'_spectrogram.mat'),'left','.'),'left',filesep);

            fid             = [pth strip(strip(data{ise,match.Audio_file},'left','.'),'left',filesep)];
            loadSpect = 0;
            if(~isempty(strfind(fid,'spectrogram.mat')))
                loadSpect   = 1;
            elseif(~isempty(ls([pth strSpect])))
                loadSpect   = 1;
                fid         = [pth strSpect];
            end

            if(loadSpect)
                disp('Loading spectrogram...');
                this_exp.audio   = matfile(fid,'Writable',true);
            else
                disp(['Processing file ' data{ise,match.Audio_file}]);
                disp('Reading audio...');
                fid             = [pth strip(strip(data{ise,match.Audio_file},'left','.'),'left',filesep)];
                [y,fs]          = audioread(fid);
                disp('Generating spectrogram...');
                win = hann(1024);
                [~,f,t,psd]     = spectrogram(y,win,[],[],fs,'yaxis');
                psd             = 10*log10(abs(double(psd)+eps));
                disp('Saving spectrogram for future use...');
                fid             = [pth strip(strip(strrep(data{ise,match.Audio_file},ext,'_spectrogram.mat'),'left','.'),'left',filesep)];
                save(fid,'-v7.3','f','t','psd','fs');
                this_exp.audio.f   = f;
                this_exp.audio.t   = t;
                this_exp.audio.psd = psd;
                this_exp.audio.fs  = fs;
            end
            disp('Done!');

            %             strtemp.audio.psd = imresize(strtemp.audio.psd,0.5);
            %             strtemp.audio.psd = strtemp.audio.psd(2:end-1,:);
            %             strtemp.audio.f   = strtemp.audio.f(3:2:end-1);
            %             strtemp.audio.t   = strtemp.audio.t(2:2:end);
            this_exp.audio.FR  = 1/(this_exp.audio.t(1,2)-this_exp.audio.t(1,1));

            %             if(hasOffset)
            %                 strtemp.audio.t  = strtemp.audio.t + offset;
            %             end
            this_exp.audio.tmin = this_exp.audio.t(1,1);
            this_exp.audio.tmax = this_exp.audio.t(1,end);
        else
            this_exp.audio = [];
        end
    end

    %% load annotations-----------------------------------------------------
    if(~isempty(data{ise,match.Annotation_file}))
        FILEEXISTS = true;
        if contains(data{ise,match.Annotation_file},'.findme')
            info=findmydata(data{ise,match.Annotation_file},params.My_Data_Folder);
            if myisfile(info.url.annot)
                [annot, tmin, tmax, allFR, this_exp.io.annot.fid, this_exp.io.annot.fidSave] = ...
                    unpackAnnotFromLoader(info.fo.annot, info.fn.annot, this_exp.annoFR, data{ise,match.Start_Anno},data{ise,match.Stop_Anno},raw{1,9});
            else
                warning('Annotation file does not exist yet. Skipping loading annotations.')
                FILEEXISTS = false;
            end
        else
            [annot, tmin, tmax, allFR, this_exp.io.annot.fid, this_exp.io.annot.fidSave] = ...
                unpackAnnotFromLoader(pth, data{ise,match.Annotation_file}, this_exp.annoFR, data{ise,match.Start_Anno},data{ise,match.Stop_Anno},raw{1,9});
        end

        if ~FILEEXISTS
            this_exp.annot           = struct();
            this_exp.annoTime        = (1:this_exp.io.movie.tmax)/this_exp.io.movie.FR;% KM???
            this_exp.io.annot        = struct();
            this_exp.io.annot.fid    = [];
            this_exp.io.annot.tmin   = 1;
            this_exp.io.annot.tmax   = length(this_exp.io.movie.tmax);
            this_exp.io.annot.FR     = this_exp.io.movie.FR;
            this_exp.annoFR          = this_exp.io.movie.FR;
            this_exp.annoFR_source   = this_exp.annoFR;

        else

            this_exp.annot           = orderfields(annot);
            tmin                    = min(tmin);
            tmax                    = max(tmax);
            [FR, this_exp.annot]     = setGlobalFR(allFR,this_exp.annot,this_exp.annoFR);
            this_exp.annoFR          = FR;
            this_exp.annoFR_source   = allFR; %for saving edited annotations back in their original FR

            if(isnan(tmax))
                tmin = this_exp.io.movie.tmin;
                tmax = this_exp.io.movie.tmax;
            end
            this_exp.io.annot.tmin   = tmin;
            this_exp.io.annot.tmax   = tmax;
            this_exp.io.annot.FR     = this_exp.annoFR;
            this_exp.annoTime        = (1:(tmax-tmin+1))/this_exp.annoFR;
            this_exp.io.movie.tmin   = tmin;
            this_exp.io.movie.tmax   = tmax;
        end

    elseif(~isempty(data{ise,match.Behavior_movie}))
        this_exp.io.annot        = struct();
        this_exp.io.annot.fid    = [];
        this_exp.io.annot.tmin   = this_exp.io.movie.tmin;
        this_exp.io.annot.FR     = this_exp.io.movie.FR;
        this_exp.annoFR          = this_exp.io.movie.FR; % change default annotation framerate to match the movie
        this_exp.io.annot.tmax   = ceil(this_exp.io.movie.tmax * this_exp.annoFR/this_exp.io.movie.FR);
        this_exp.annoFR_source = this_exp.annoFR;
        this_exp.annot           = struct();
        this_exp.annoTime        = (this_exp.io.annot.tmin:this_exp.io.annot.tmax)/this_exp.io.annot.FR;

    elseif(~isempty(data{ise,match.Audio_file}))
        this_exp.io.annot        = struct();
        this_exp.io.annot.fid    = [];
        this_exp.io.annot.tmin   = this_exp.audio.tmin;
        this_exp.io.annot.tmax   = ceil(this_exp.audio.tmax * this_exp.annoFR/this_exp.audio.FR);
        this_exp.io.annot.FR     = this_exp.audio.FR;
        this_exp.annoFR          = this_exp.io.movie.FR; % change default annotation framerate to match the audio
        this_exp.annoFR_source = this_exp.annoFR;
        this_exp.annoTime        = this_exp.io.annot.tmin:(1/this_exp.annoFR):this_exp.io.annot.tmax;
        this_exp.annot           = struct();

    elseif(~isempty(data{ise,match.Tracking}) && ~isempty(this_exp.trackTime))
        this_exp.io.annot        = struct();
        this_exp.io.annot.fid    = [];
        this_exp.io.annot.tmin   = 1;
        this_exp.io.annot.tmax   = length(this_exp.trackTime);
        this_exp.io.annot.FR     = 1/(this_exp.trackTime(2)-this_exp.trackTime(1));
        this_exp.annoFR          = 1/(this_exp.trackTime(2)-this_exp.trackTime(1));
        this_exp.annoFR_source = this_exp.annoFR;
        this_exp.annoTime        = this_exp.io.annot.tmin:(1/this_exp.annoFR):this_exp.io.annot.tmax;
        this_exp.annot           = struct();
    else
        this_exp.annot           = struct();
        this_exp.annoTime        = (1:this_exp.io.movie.tmax)/this_exp.io.movie.FR;% KM???
        this_exp.io.annot        = struct();
        this_exp.io.annot.fid    = [];
        this_exp.io.annot.tmin   = 1;
        this_exp.io.annot.tmax   = length(this_exp.io.movie.tmax);
        this_exp.io.annot.FR     = this_exp.io.movie.FR;
        this_exp.annoFR          = this_exp.io.movie.FR;
        this_exp.annoFR_source   = this_exp.annoFR;
    end

    try
        mouse(data{ise,match.Mouse}).(['session' num2str(data{ise,match.Sessn})])(data{ise,match.Trial}) = this_exp;
    catch
        keyboard
    end
end
