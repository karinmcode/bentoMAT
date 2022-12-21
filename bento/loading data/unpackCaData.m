function [rast,time,spikes,ROIs] = unpackCaData(URL)
%
% (C) Ann Kennedy, 2019
% California Institute of Technology
% Licensing: https://github.com/annkennedy/bento/blob/master/LICENSE.txt



[MyDataFolder,fname,ext] = fileparts(URL);
disp(['Loading Ca file ' fname '...']);

[rast,time,ROIs]=deal([]);
% add cases to this switch statement for other data types~?
switch ext
    case '.csv'
%         fid = fopen(pth);
%         temp = textscan(fid,'%f,%f,[%d %d],%d,%f');
%         neurons = unique(temp{1});
%         rast = [];
%         for i = neurons'
%             rast(i+1,:) = temp{6}(temp{1}==i);
%             rast(i+1,temp{5}(temp{1}==i)==0)=nan;
%         end
        temp = readtable(URL);
        if(size(temp,2)==5) % it's an Anderson lab miniscope file
            time = temp(:,1)';
            rast = temp(:,5)';
            % get rid of artifacts
            rast((1:500))    = nan;
            rast(end-119:end) = nan;
        else % it's probably an inscopix csv
            rast = table2array(temp)';
            time = (rast(1,:)+rast(1,2));
            rast = rast(2:end,:);
        end
        spikes=[];
    case {'.findme' '.me'}%KM  Karin's files
        % think of adding a info=findmydata(filename) function
%         [a , d, rec]=mygce(fname);% example : fname = 'CAM1_m523_211118_003_motion_speed'% gce
%         DATA_FOLDER = fullfile(MyDataFolder,a,d,rec);% example Y:\Users\Karin\data\processed\aligned2vid\m523\211118\003
%         camID = mycamID(fname);
%         CaFile = sprintf('fluo_CAM%g.mat',camID);
%         CaUrl = fullfile(DATA_FOLDER,CaFile);

        info=findmydata(URL);%URL should contains camID and DATA_FOLDER ideally
        load(info.url.phys,'fluo');%mywinopen(CaUrl)
        rast = double(fluo.values_common);
        time = fluo.time_common- fluo.time_common(1);% here the code assumes that the calcium trace and the video are the same duration
        %ROIs coming!
        spikes=[];
    case '.mat'
        temp = load(URL);
        f    = fieldnames(temp);
        if(length(f)==1 || isfield(temp,'neuron'))
            if(length(f)==1)
                use = f{1};
            else
                use = 'neuron';
            end
            S2d = strcmpi(class(temp.(use)),'Sources2D');
            if(S2d || (isstruct(temp.(use)) && isfield(temp.(use),'C_raw'))) %CNMFE data
                rast    = temp.(use).C_raw;
                if(S2d || (isfield(temp.(use),'options') && isfield(temp.(use),'A')))
                    d1      = temp.(use).options.d1;
                    d2      = temp.(use).options.d2;
                    ROIs = zeros(size(rast,1),d1,d2);
                    for i = 1:size(rast,1)
                        ROIs(i,:,:) = reshape(temp.(use).A(:,i),d1,d2);
                    end
                end
                if(S2d || isfield(temp.(use),'S'))
                    spikes  = temp.(use).S;
                else
                    spikes = [];
                end
            else
                rast = double(temp.(use)); %assume it's a matrix of traces
                spikes = [];
            end
            time = [];
        elseif(length(f)==2 && any(strcmpi(f,'time')))
            rast = temp.(f{~strcmpi(f,'time')});
            spikes = [];
            time = temp.(f{strcpmi(f,'time')});
        elseif any(strcmpi(f,'C')) && isstruct(temp.C) && isfield(temp.C,'data') % it's probably minian output
            rast = temp.C.data;
            spikes = [];
            time = []; %look for miniscope_timeStamps.csv in same folder, and video_timeStamps.csv for behavior video
        else
            disp(['unsure which variable to read in ' fname]);
            rast = []; spikes=[];
        end
        time = [];
    case '.flr'
        % Yatang Ca traces
        datapara = input_lyt(URL);
        rast = datapara.data';
        time=[];
        spikes=[];
end

if isempty(time)
    time = getVideoTimestamps(URL);
    if(time)
        if(size(time,2)==1)
            time=time';
        end
    end
end

disp('done');