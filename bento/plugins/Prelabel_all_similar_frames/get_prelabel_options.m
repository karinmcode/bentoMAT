function [opt,isloading_necessary,urlstep] = get_prelabel_options(gui,varargin)
%% Define steps

%% 1) Image processing steps
% - vidmotion indicates that the motion is already in the video
% - HOG
% - ROI
%% 2) Dimensionality reduction steps
% - PCA
% - tSNE
% - uMAP

%% 3) Clustering steps
% - uMAP output
% - dbscan
% - watershed
% - kmeans
% - fcm
% - gaussmix :gaussian mixture model clustering

if ~isempty(varargin)
    opt =  varargin{1};
else
    opt.steps = {'vidmotion' 'HOG','PCA','uMAP','fcm'};%'ROI','HOG','wavelet','PCA','tSNE','kmeans',dbscan,watershed,fcm
end

fprintf('\n\n   %s',get_params_string(opt));
nsteps = numel(opt.steps);
opt.DimRedName = 'none'; % initialize dimensionality reduction name variable
info = findmydata(gui);
proc_data_folder = info.fo.proc;
currFN0 = info.fn0.vid;

%% -----
%% PROCESSING STEPS PARAMETERS
for ist = 1:nsteps
    stepname = opt.steps{ist};
    switch stepname
        case 'ROI'

            opt.ROI.name = 'eye';
            FN0.(stepname) = [currFN0 '__ROI_' opt.ROI.name ]
            DoPlot.(stepname) = 1;
            currFN0 = FN0.(stepname);

            % make urls
            urlstep.(stepname) = fullfile(proc_data_folder,[currFN0,'.mat']);
            stepfile_exists(ist)=myisfile(urlstep.(stepname));
            continue;
        case 'motion'
            opt.motion.nfr = 2;
        case 'motioncat'
            opt.motioncat.nfr = 2;
        case 'motionabs'
            opt.motionabs.nfr = 2;
        case 'vidmotion'
            opt.vidmotion=struct();
        case 'HOG'
            opt.HOG.CellSize = [32 32];
            opt.HOG.CellSize = 20;% divide im by 20
            FN0.(stepname) = [currFN0 '__' stepname '_' replace(get_params_string(opt.(stepname)),'__','_')];
            DoPlot.(stepname) = 0;
       
        case 'PCA'
            opt.PCA.nPCs = 100;
            FN0.(stepname) = [currFN0 '__' stepname '_' replace(get_params_string(opt.(stepname)),'__','_')];
            DoPlot.(stepname) = 0;
        case 'tSNE'
            testClus = 0;

            opt.tSNE.NumPCAComponents = 3;% should not exceed the number of col in obs (nfeatures or 100 PCs, here)
            opt.tSNE.Exaggeration = 4;
            opt.tSNE.LearnRate = 250;
            opt.tSNE.Perplexity = 20;
            FN0.(stepname) = [currFN0 '__' stepname '_' replace(get_params_string(opt.(stepname)),'__','_')];
            DoPlot.(stepname) = 0;
            %% -   uMAP
        case 'uMAP'
            testClus = 0;
            opt.uMAP.min_dist = 0.06;%default 0.3, compaction 0<mindist<1
            opt.uMAP.n_neighbors = 199;%default 15, 2<n_neighbors<199
            opt.uMAP.template=0;

            FN0.(stepname) = [currFN0 '__' stepname '_' replace(get_params_string(opt.(stepname)),'__','_')];
            DoPlot.(stepname) = 1;
        case 'kmeans'
            opt.kmeans.nclusters = 16;
            FN0.(stepname) = [currFN0 '__' stepname '_' replace(get_params_string(opt.(stepname)),'__','_')];

        case 'dbscan'
            opt.dbscan.MinPoints = 5;% Minimum number of neighbors for a core point % set to nan if you want it computed by code
            opt.dbscan.Radius =0.25;% set to nan if you want it computed by code
            %opt.dbscan.nclusters = 16;
            % bestRadius=clusterDBSCAN.estimateEpsilon(xy,5,10);

            FN0.(stepname) = [currFN0 '__' stepname '_' replace(get_params_string(opt.(stepname)),'__','_')];
        case 'watershed'
            opt.watershed.nbins = [20 20];% 
        case 'fcm'
            opt.fcm.nclusters = 16;
        case 'gaussmix'
            opt.gaussmix.nclusters = 16;
    end


    FN0.(stepname) = [currFN0 '__' stepname '_' replace(get_params_string(opt.(stepname)),'__','_')];
    currFN0 = FN0.(stepname);

    %% make urls
    urlstep.(stepname) = fullfile(proc_data_folder,[currFN0,'.mat']);
    if numel(currFN0)>255
        keyboard;
    end
    stepfile_exists(ist)=myisfile(urlstep.(stepname));
end

for ist=1:nsteps-1
    stepname = opt.steps{ist};
    isloading_necessary.(stepname) =  stepfile_exists(ist+1)~=1;
end
isloading_necessary.(opt.steps{end})=true;
isloading_necessary.PCA = true;
