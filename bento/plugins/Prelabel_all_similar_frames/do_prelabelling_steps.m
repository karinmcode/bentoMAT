function [new_labels,new_annotations,data,opt] = do_prelabelling_steps(gui,opt,isloading_necessary,urlstep)
dbstop if error
%% Define Variables
Vin = gui.data.io.movie.reader{1}.reader;
w = Vin.Width;
h = Vin.Height;
motionStep =0;% irrelevant in this context
DoPlot.PCA = 0;
DoPlot.uMAP = 0;
%% GO THROUGH ALL STEPS
%
%% STEP ROI
if contains('ROI',opt.steps)

    % make ROI figure
    fClus=makegoodfig('fig_ROI','slide');
    add_analysis_params(gcf,opt);
    Vin.CurrentTime = 0;
    thisFrame = readFrame(Vin);
    I = thisFrame(:,:,1);
    image(thisFrame);%thisFrame(1,:,1)figure;image(thisFrame/mymax(thisFrame))
    goodax(gca,'axis',{'equal' 'tight'},'title',{sprintf('%s %s',Vin.Name,opt.ROI.name)},'box','off','colormap',gray(255));


    if isfile(urlstep.ROI)
        load(urlstep.ROI,'ROIpos')
        rect=drawrectangle('position',ROIpos,'label',opt.ROI.name);
    else
        % choose from existing ROI rand ROI name
        rect=drawrectangle();
        title('double-click on rectangle to continue');
        wait(rect);
        ROIpos = rect.Position;
        save(urlstep.ROI,'ROIpos');

    end
    imcropped = imcrop(I,ROIpos);
    h = size(imcropped,1);
    w = size(imcropped,2);

end


%% STEP : HOG
if contains('HOG',opt.steps)
    fprintf('\n==== HOG ====')
    % - Compute HOG
    if myisfile(urlstep.HOG)
        sz=fprintf('\nLoading HOG ... %s',getfilesize(urlstep.HOG));
        load(urlstep.HOG,'HOG')
    else

        NFR = Vin.NumFrames;
        Vin.CurrentTime = 0;
        frames_vec = 1:NFR;

        if numel(opt.HOG.CellSize)==2
            HOGCellSize = opt.HOG.CellSize;
        else
            if opt.HOG.CellSize<1
                HOGCellSize = round([h w]*opt.HOG.CellSize);
            else
                HOGCellSize = round([h w]/opt.HOG.CellSize);
            end
        end

        % find number of features to preallocate
        thisFrame = readFrame(Vin);
        I = thisFrame(:,:,1);

        if contains('ROI',opt.steps)
            I= imcrop(I,ROIpos);
        end

        if contains('motion',opt.steps)
            thisFrame = getframes(Vin,1:2);
            I = diff(thisFrame,1,3);
            dI = diff(thisFrame,1,3);
        elseif contains('vidmotion', opt.steps)% motion already in video
            thisFrame = readFrame(Vin);
            I = thisFrame(:,:,1)+thisFrame(:,:,2);
        else
            I = thisFrame;
        end

        if contains('ROI',opt.steps)
            I= imcrop(I,ROIpos);
        end

        [features,hogVisualization]=extractHOGFeatures(I,'CellSize',HOGCellSize);

        %figure
        makegoodfig('preproc2HOG','slide');
        add_analysis_params(gcf,opt);
        ax=axes();

        imh=image(thisFrame);%thisFrame(1,:,1)figure;image(thisFrame/mymax(thisFrame))
        goodax(ax,'axis',{'equal' 'tight'},'title',{sprintf('Frame%g',1)},'box','off','hold','on','colormap','gray','ydir','reverse');

        if motionStep
            AlphaData = double(abs(dI)+0.2);
            AlphaData(AlphaData>1)=1;
            imh.AlphaData = AlphaData;
        end

        pHOG=plot(hogVisualization);
        pHOG(1).Color = [1 1 1 ]*0.5;
        nfeatures = numel(features);
        HOG = nan(NFR,nfeatures);
        title({sprintf('nfeatures = %g (Opt.HOGsize=%s)',nfeatures,num2str(opt.HOG.CellSize))})
        nCells = floor(size(I,[1 2])./HOGCellSize);
        ylabel(ax,sprintf('%g row cells',nCells(1)))
        xlabel(ax,sprintf('%g column cells',nCells(2)))

        % go across all frames
        tic
        Vin.CurrentTime = 0;
        sz=0;
        for ifr = frames_vec%ifr=1
            cleanline(sz);
            sz=fprintf('Frame %g/%g   > %.0f%%',ifr,NFR,100*ifr/NFR);



            if contains('motion', opt.steps)
                if ifr~=NFR
                    thisFrame = getframes(Vin,ifr+(0:1));
                    I = diff(thisFrame,1,3);
                else
                    I = diff(thisFrame,1,3);
                end
            elseif contains('motioncat', opt.steps)
                if ifr~=NFR
                    I = getframes(Vin,ifr+(0:1));
                    I(:,:,3) = I(:,:,1);
                else
                    I = getframes(Vin,ifr-1+(0:1));
                    I(:,:,3) = I(:,:,1);
                end
            elseif contains('motionabs', opt.steps)
                if ifr~=NFR
                    thisFrame = getframes(Vin,ifr+(0:1));
                    I = abs(diff(thisFrame,1,3));
                else
                    I = abs(diff(thisFrame,1,3));
                end
            elseif contains('vidmotion', opt.steps)% motion already in video
                thisFrame = readFrame(Vin);
                I = thisFrame(:,:,1)+thisFrame(:,:,2);
            else
                thisFrame = readFrame(Vin);
                I = mat2gray(thisFrame);
            end

            if contains('ROI',opt.steps)
                I= imcrop(I,ROIpos);
            end

            features=extractHOGFeatures(I,'CellSize',HOGCellSize,'NumBins',9);%21*21*9*3
            %makegoodfig('preproc2HOG_temp','slide');plot(features)

            % store data
            HOG(ifr,:)=features;

        end


        frame_rate = Vin.FrameRate;
        nframes = NFR;
        img_size = size(I);
        mymkdir(urlstep.HOG);
        save(urlstep.HOG,'-v7.3','HOG','frame_rate','nframes','img_size','HOGCellSize');

    end

    obs = HOG;



    clear HOG;
else% upload all frames
    all_frames = getframes(Vin);
    obs = all_frames;
end

%% STEP: PCA
if contains('PCA',opt.steps)
    fprintf('\n==== PCA ====');
    if isloading_necessary.PCA
        tic
        if exist(urlstep.PCA,'file')==0
            sz=fprintf(' >> doing PCA ...') ;
            tic
            [coeff,score,latent,tsquared,explained,mu] = pca(obs,'NumComponents',opt.PCA.nPCs);
            cleanline(sz);

            % for explained plot
            %             obs_S = zscore(obs,0,1);
            %             [U,S,V] = svd(obs_S,'econ');%taakes too long

            fprintf(' >> did PCA in %.0f sec',toc) ;
            nPCs = opt.PCA.nPCs;
            save(urlstep.PCA,'-v7.3','score','nPCs','explained','coeff','explained','latent')
        else
            sz=fprintf('\nLoading PCA   ... %s',getfilesize(urlstep.PCA));
            load(urlstep.PCA,'score','nPCs','explained','coeff','explained','latent')
        end
        cleanline(sz);
        fprintf('\nGot %g first PCs scores in %.0f sec (%.0f min).',opt.PCA.nPCs,toc,toc/60)

        obs = score;
        if DoPlot.PCA
            % figure
            fPCA = makegoodfig('figPCA','slide');
            add_analysis_params(gcf,opt);

            ax = axes('position',[0.1 0.25 0.8 0.55]);
            NFR = Vin.NumFrames;
            fs = Vin.FrameRate;
            tcomp = linspace(0,(NFR/fs)-1/fs,NFR);
            [ax,ax2,pl]=myplotxxn(ax,1:NFR,tcomp,score(:,1:3),'xlabel',{'Frames' 'Time (s)'});
            set([ax ;ax2 ],'PositionConstraint','innerposition')
            legstr = mynum2str(1:3,'PC%g','cellstr');
            leg = legend(pl,legstr);
            goodax(ax,'ylabel','PCA scores','box','off');

%             savethisfig(fig_folder,fPCA,opt,fn0vid)

            % figure PCA explained
            GrpId = ones(NFR,1);
            varnames = 1:numel(explained);
            GrpNames = {'g1'};
            fPCAexpl= mypcaplots([],explained,GrpId,GrpNames,coeff,score,varnames);%f= mypcaplots(S,explained,ObsIds,coeff,score,varnames)
            add_analysis_params(gcf,opt);
            %% - PCA best var
            I = getframes(Vin,1);
            cell_size = opt.HOG.CellSize;
            nbOfBinsPerCell = 9;
            BlockSize = [2 2];
            if numel(opt.HOG.CellSize)==2
                HOGCellSize = opt.HOG.CellSize;
            else
                if opt.HOG.CellSize<1
                    HOGCellSize = round([h w]*opt.HOG.CellSize);
                else
                    HOGCellSize = round([h w]/opt.HOG.CellSize);
                end
            end
            BlockOverlap =1;
            nTopVar = 5;
            TopPCAVar = getTopPCAVar(coeff,nTopVar);
            nPCs_plot = 10;

            %figure
            fPCA_bestvar = makegoodfig('figPCA_bestvar','slide');
            ax = axes;
            imagesc(I)
            goodax(ax,'axis',{'equal' 'tight'} ,'colormap','gray','xyticklabel','','hold','on');

            CMblock = jet(nPCs_plot);

            for ipc = 1:nPCs_plot
                idxx = TopPCAVar(:,ipc);
                for idx = idxx(:)';
                    pCell = 6;%percentExplained
                    cellpos=HOGind2framecell(idx,I,HOGCellSize,nbOfBinsPerCell,BlockSize,BlockOverlap);
                    rect=drawrectangle(ax,'Position',cellpos);
                    set(rect,'label',['PC' num2str(ipc) '-' num2str(idx)],'color',CMblock(ipc,:),'linewidth',1,'MarkerSize',1,'FaceAlpha',0.3,'LabelAlpha',0.3)
                end
            end
            pExplained = cumsum(explained(1:nPCs));
            title(sprintf('showing %g best features for %g first PCs (%%explained = %.0f)',nTopVar, nPCs_plot,pExplained(nPCs_plot)));
            p = ax.Position;
            cb=mycolorbar(ax,CMblock,mynum2str(1:nPCs_plot,'PC%g','cellstr'));
            cb.Position(1)=p(1)+p(3)*0.9;
            add_analysis_params(gcf,opt);
            %% - pca fft
            fPCA_fft = makegoodfig('figPCA_fft','slide');
            nrow = 4;
            subplot(nrow,1,1,'replace');
            normscores = score(:,1:nPCs)';
            normscores = normscores./mymax(abs(normscores),2);
            imagesc(normscores)
            goodax(gca,'xlabel','frames','ylabel',sprintf('%g PCs',opt.PCA.nPCs));
            caxis([-1 1]);
            cb=colorbar(gca);
            ylabel(cb,'norm PCA scores');

            %
            for iPC = 1:3;
                subplot(nrow,1,iPC+1,'replace','colormap',turbo(100));
                signal = score(:,iPC);
                pspectrum(signal,fs,'spectrogram', ...
                    'TimeResolution',1,'FrequencyLimits',[0.01 5]);
                colormap turbo;
                title(sprintf('PC%g   ',iPC))
            end
            add_analysis_params(gcf,opt);
        end
    end
end


%% STEP DIMENSIONALITY REDUCTION: tSNE
if contains('tSNE',opt.steps)
    opt.DimRedName = 'tSNE';
    if isloading_necessary.tSNE
        % tsne params
        rng default % for reproducibility
        testClus=0;
        if testClus
            clu_params=test_tSNEinputs(obs,size(obs,1));
            keyboard;
        else
            clu_params = keepfield(opt.tSNE,{'NumPCAComponents','Exaggeration','LearnRate','Perplexity'});
        end

        % tsne
        tic
        rng default % for reproducibility
        XYclusters = tsne(obs,...
            'Algorithm','barneshut','NumPCAComponents',clu_params.NumPCAComponents,'Exaggeration',clu_params.Exaggeration,'LearnRate',clu_params.LearnRate,'perplexity',clu_params.Perplexity);
        fprintf(' >> did t-SNE in %.0f sec',toc) ;
    end
end

%% STEP DIMENSIONALITY REDUCTION: uMAP
if contains('uMAP',opt.steps)
    addpath(genpath('G:\My Drive\code\downloaded'))
    opt.DimRedName = 'uMAP';
    opt.ClusteringMethod = 'uMAP';
    testClus =0;
    if isloading_necessary.uMAP
        tic;

        if testClus==1
            fig_test_uMAP=test_uMAP(obs);
        end

        if isfile(urlstep.uMAP)==0
            nvar = size(obs,2);
            uMAP_var = mynum2str(1:nvar,'v%g','cellstr');

            min_dist=opt.uMAP.min_dist;
            n_neighbors=opt.uMAP.n_neighbors;

            if opt.uMAP.template==0
                [XYclusters,umap,clusterIdentifiers,extras]=run_umap(obs,'parameter_names',uMAP_var,'min_dist',min_dist,'n_neighbors',n_neighbors,'randomize',false);
                close(extras.fig);
            else
                url_template_uMAP = fullfile('Y:\Users\Karin\data\processed\videos\batch_2p_211108_all_cells\get_beh_states_220608\m522\211112\002','m522_211112_002_CAM2__HOG_CellSize_20__PCA_nPCs_100__uMAP_min_dist_0o06_n_neighbors_42.mat');
                [XYclusters,umap,clusterIdentifiers,extras]=run_umap(obs,'parameter_names',uMAP_var,'min_dist',min_dist,'n_neighbors',n_neighbors,'randomize',false,...
                    'template_file',url_template_uMAP,'see_training',true);
                close(extras.fig);
            end
            try
                save(urlstep.uMAP,'-v7.3','XYclusters','umap','clusterIdentifiers');
            catch err
                warning('did not save');
                warning(err.message);
            end
        else
            load(urlstep.uMAP,'XYclusters','umap','clusterIdentifiers');
        end
        clusterIdentifiers = clusterIdentifiers(:)+1;
        p_clusters = unique(clusterIdentifiers);
        nclusters = numel(p_clusters);
        if nclusters==1
            keyboard
        end
        Centroids=mycentroids(XYclusters(:,1),XYclusters(:,2),clusterIdentifiers);



        fprintf(' >> did uMAP in %.0f sec (%.0f min)',toc,toc/60) ;
        if DoPlot.uMAP
            figUMAP = makegoodfig('uMAP','slide_half_width');
            add_analysis_params(gcf,opt);
            nrow = 5;
            ax=subplot(nrow,1,2:3,'replace');
            NFR = size(XYclusters,1);
            gs=scatter(XYclusters(:,1),XYclusters(:,2),ones(NFR,1)*6,clusterIdentifiers,'filled');
            goodax(ax,'colormap','turbo','colorbar',{'ylabel','cluster ids'},'axis',{'equal' 'tight'},'title',{{'uMAP',sprintf('min dist = %g , n neighbors = %g , n found clus = %g',opt.uMAP.min_dist,opt.uMAP.n_neighbors,nclusters)}});


            ax=subplot(nrow,1,4:5,'replace');
            gs=scatter(XYclusters(:,1),XYclusters(:,2),ones(NFR,1)*6,1:NFR,'filled');
            goodax(ax,'colormap','turbo','axis',{'equal', 'tight'},'colorbar',{'ylabel','frame idx'});
        end


    end
end

%% ------------------------------------------------
%% STEP : CLUSTER IDENTIFICATION kmeans
if contains('kmeans',opt.steps)
    opt.ClusteringMethod = 'kmeans';

    rng default % for reproducibility
    dim_red_position = find(ismember(opt.steps,opt.DimRedName));
    clusterOnPCA = false;
    if ~isempty(dim_red_position)
        clusterOnPCA = find(strcmp(opt.steps,'kmeans'))<dim_red_position;
    end

    if clusterOnPCA
        [clusterIdentifiers,Centroids]=kmeans(XYclusters,opt.kmeans.nclusters,"Distance","cityblock",'start','cluster','MaxIter',500);
    else
        [clusterIdentifiers,Centroids]=kmeans(obs,opt.kmeans.nclusters,"Distance","cityblock",'start','cluster','MaxIter',500);
    end
    nclu = opt.kmeans.nclusters;
end

%% STEP : CLUSTER IDENTIFICATION dbscan

if contains('dbscan',opt.steps)
    opt.ClusteringMethod = 'dbscan';

    rng default % for reproducibility
    MinPoints = opt.dbscan.MinPoints;% Minimum number of neighbors for a core point ndim+1
    Radius =opt.dbscan.Radius;

    if isnan(MinPoints)
        MinPoints = 5;
    end

    if isnan(Radius)
        thresholdDiff = 0.1;% vary between 0.01 - 0.1 ; low means more clusters
        %Radius = nearDist(find(diff(nearDist)/max(diff(nearDist))>0.05,1));
        Radius=clusterDBSCAN.estimateEpsilon(XYclusters,MinPoints,MinPoints+5);
        clusterDBSCAN.estimateEpsilon(XYclusters,MinPoints,MinPoints+5);

        % code i wrote
%     % deifne radius
%     allSmallestDistances = pdist2(XYclusters,XYclusters,'euc','Smallest',MinPoints);% look for the beggining of the knee
%     nearDist = sort(allSmallestDistances(end,:));
%     makegoodfig('dbscan_knee');
%     subplot(2,1,1);
%     plot(nearDist);
%     hold on;
%     plot(xlim,[1 1]*Radius,':r')
%     axis square;
%     title('k-distance graph')
%     xlabel('Points sorted with 50th nearest distances')
%     ylabel('50th nearest distances')
%     add_analysis_params(gcf,opt);
%     grid;

    end
    opt.dbscan.Radius = round(Radius,2);



    % dbscan
    if ~isfield(opt.dbscan,'nclusters')
        clusterIdentifiers = dbscan(XYclusters,Radius,MinPoints);% (Radius,minpoints)
        if ismember(-1,clusterIdentifiers)
            clusterIdentifiers(clusterIdentifiers==-1)=mymax(clusterIdentifiers)+1;
        end
        pclu = unique(clusterIdentifiers);
        nclu = numel(pclu);
    else
        nclu = opt.dbscan.nclusters;
        %H=my_test_dbscan(XYclusters)
        [clusterIdentifiers,Radius,MinPoints,nclu] = mydbscan(XYclusters,nclu);
       
    end
    Centroids = mycentroids(XYclusters(:,1),XYclusters(:,2),clusterIdentifiers);

    subplot(2,1,2);

    gs=gscatter(XYclusters(:,1),XYclusters(:,2),clusterIdentifiers);
    legend(gs,'location','bestoutside')
end
%% STEP : CLUSTER IDENTIFICATION ---> watershed

if contains('watershed',opt.steps)
    fprintf('\n==== watershed ====')
    opt.ClusteringMethod = 'watershed';
    [XYdensity, centers] = hist3(XYclusters,'nbins',opt.watershed.nbins);%Opt.watershed.nbins = [21 21]

    % figure;
    fwat=makegoodfig('watershed','slide');
    add_analysis_params(gcf,opt);
    ncol = 3;
    ax=subplot(1,ncol,1,'replace')
    imagesc(centers{:}, XYdensity.')
    if contains(opt.steps,'uMAP')
        pclu = unique(clusterIdentifiers);
        nclu = numel(pclu);
        goodax(ax,'axis',{'equal' 'tight'},'ydir','normal','title',sprintf('nclu dimred= %g',nclu),'colormap','hot','colorbar',{'ylabel','density'});
    else
        goodax(ax,'axis',{'equal' 'tight'},'ydir','normal','colormap','hot','colorbar',{'ylabel','density'});
    end
    % watershed
    opt.watershed.pixelConnectivity = 4;
    RidgeLinesImg = watershed(XYdensity,opt.watershed.pixelConnectivity);% (Radius,minpoints)
    hold on

    %ax2 watershed image
    ax=subplot(1,ncol,2,'replace');
    imw = imagesc(centers{1},centers{2},RidgeLinesImg);
    imw.AlphaData = RidgeLinesImg~=0;
    pRidges = unique(RidgeLinesImg);
    nRi = numel(pRidges);
    goodax(ax,'title',sprintf('nclu watershed = %g',nRi),'colormap','turbo','ydir','normal','axis',{'equal' 'tight'});

    % ax with points colored by cluid
    ax=subplot(1,ncol,3,'replace');
    xclu = XYclusters(:,1);
    yclu = XYclusters(:,2);
    clusterIdentifiers=mydiscretize2(xclu,yclu,centers{1},centers{2},RidgeLinesImg);
    Centroids = mycentroids(xclu,yclu,clusterIdentifiers);
    CMclu = turbo(nRi);
    sc=scatter(xclu,yclu,yclu*0+6,clusterIdentifiers,'filled','marker','o');
    goodax(ax,'colormap',CMclu,'caxis',myminmax(pRidges),'colorbar',{'ylabel','watershed segments ID'},'axis',{'equal','tight'});

end

%% STEP : CLUSTER IDENTIFICATION ---> fcm fuzzy c-means clustering.

if contains('fcm',opt.steps)%
    %{
URLxy = 'G:\My Drive\code\GUIs\bentoMAT-master\data\ExptKM\proc_data\CAM1_m943_220413_003_motion_speed__vidmotion___HOG_CellSize_20__PCA_nPCs_100__uMAP_min_dist_0o06_n_neighbors_199_template_0.mat'
load(URLxy,'XYclusters','clusterIdentifiers');
opt.fcm.nclusters = 16
    %}

    [Centroids,U] = fcm(XYclusters,opt.fcm.nclusters,[2 100 1e-5 false]);

    % U , Fuzzy partition matrix, returned as a matrix with Ncenters rows and Nd columns.
    % Element U(i,j) indicates the degree of membership of the jth data point in the ith cluster.
    % For a given data point, the sum of the membership values for all clusters is one
    maxU = max(U);
    nclu = opt.fcm.nclusters;
    nobs = size(XYclusters,1);
    clusterIdentifiers = nan(nobs,1);
    for iclu = 1:nclu
        index = find(U(iclu,:) == maxU);
        clusterIdentifiers(index) = iclu;

    end
    %{
 %plotting
makegoodfig('fuzzy partition matrix');
CM = jet(nclu);
for iclu = 1:nclu
i4clu = clusterIdentifiers==iclu;
p(1)=plot(XYclusters(i4clu,1),XYclusters(i4clu,2),'.b');
hold on

p(2)=plot(Centroids(iclu,1),Centroids(iclu,2),'xb','MarkerSize',15,'LineWidth',3);
set(p,'color',CM(iclu,:));
end
hold off
    %}
    opt.ClusteringMethod = 'fcm';

end
%% CREATE DATA OUTPUT
data.pclu = unique(clusterIdentifiers);
data.nclu = numel(data.pclu);
data.cluids = clusterIdentifiers;
data.clu_xy = XYclusters;
new_labels = cellfun(@(x) sprintf('prelabel_%g',x) ,num2cell(data.pclu),'uniformoutput',false);
new_annotations = cellfun(@(x) sprintf('prelabel_%g',x) ,num2cell(clusterIdentifiers),'uniformoutput',false);

%% -------------------------------
% --------------------------------------
% --------------------------------------
% --------------------------------------
% --------------------------------------
% --------------------------------------
% --------------------------------------
% --------------------------------------
% --------------------------------------
% --------------------------------------
% --------------------------------------
% --------------------------------------

