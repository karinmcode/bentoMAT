function figs = do_prelabelling_figs(gui,opt,data)

Vin = gui.data.io.movie.reader{1}.reader;
NFR = Vin.NumFrames;
url_vid = fullfile(Vin.Path,Vin.Name);
% get possible cluster ids
pclu = data.pclu;
nclu = numel(pclu);
cluids = data.cluids;
cluxy =data.clu_xy;
ROIpos = [0 0 0 0];
if nclu==1
    disp('only 1 cluster')
    keyboard
end
figs=gobjects(0,1);

%% plot clusters
if any(pclu==-1)
    GrpCM = [0 0 0;turbo(nclu-1)];
else
    GrpCM = turbo(nclu);
end




%% FIGURE clustersmontage
plotFig_clustersmontage = 1;
if plotFig_clustersmontage
fig=makegoodfig('clustersmontage','slide');
figs = vertcat(figs,fig);
add_analysis_params(gcf,opt);

AX(1)=subplot(2,2,[1 3],'replace');
sc=scatter(cluxy(:,1),cluxy(:,2),cluxy(:,2)*0+6,cluids,'filled');
goodax(AX(1),'colormap',GrpCM,'axis',{'equal' 'tight'},'title',{{opt.DimRedName,opt.ClusteringMethod,'clustering results',sprintf('nclu = %g',nclu)}},'color',[1 1 1]*rem(now,1),'colorbar',{'ylabel','clusterIDs'});

% plot cluster montage

h=Vin.Height;
w = Vin.Width;
clu_examples= nan(Vin.Height,Vin.Width,3,nclu);
tic
Margin = round(h/30);
fprintf(' >> \nFinding example images for each cluster...') ;
    motionStep = 0;

for iclu = 1:nclu

    i4clu = find(cluids==iclu);
    if isempty(i4clu)
        continue;
    end

    if motionStep
        i4clu(i4clu==NFR)=[];
        i4clu = randsample(i4clu,1);
        Iexamples=getframes(Vin,i4clu);
        Iexamples2=getframes(Vin,i4clu+1);
        Iexamples = cat(4,Iexamples,Iexamples2,Iexamples);
        Iexamples = permute(Iexamples,[1 2 4 3]);
        Iexamples = mean(Iexamples,4);
        dI = abs(diff(Iexamples(:,:,1:2),1,3));
        dI3 = repmat(dI,1,1,3);
        i4dI3 = dI3<0.1;
        Iexamples(i4dI3)= Iexamples(i4dI3)+dI3(i4dI3);
        Iexamples(Iexamples>1)=1;
    else
        try
        i4clu = randsample(i4clu,3);
        end
        Iexamples=getframes(Vin,i4clu,1);
        Iexamples = mymean(Iexamples,4);
        Iexamples = Iexamples/mymax(Iexamples(:));% normalize
    end

    if contains('ROI',opt.steps)
        Iexamples= imcrop(Iexamples,ROIpos);
    end
    cm_clu = permute(GrpCM(iclu,:),[1 3 2]);
   
    Iexamples(1:Margin,:,:)= repmat(cm_clu,Margin,w);
    Iexamples(:,1:Margin,:)= repmat(cm_clu,h,Margin);
    Iexamples(end-(Margin-1):end,:,:)= repmat(cm_clu,Margin,w);
    Iexamples(:,end-(Margin-1):end,:)= repmat(cm_clu,h,Margin);

    clu_examples(:,:,1:3,iclu)= Iexamples;
end

clu_examples= (clu_examples-mymin(clu_examples(:)))/(mymax(clu_examples(:))-mymin(clu_examples(:)));
fprintf(' >> found example images in %.0f sec',toc) ;

% plot example clusters
AX(2)=subplot(2,2,2,'replace');
ncol = 4;
nrow = round(nclu/ncol);
while ncol*1.5<nrow
    ncol = ncol+1;
    nrow = round(nclu/ncol);
end
montage(clu_examples,'Size',[nrow ncol]);
goodax(AX(2),'axis',{'equal'},'title',Vin.Name);

AX(3)=subplot(2,2,4,'replace');
montage(clu_examples,'Size',[nrow ncol]);
goodax(AX(3),'axis',{'equal'});
axis([0 w 0 h]*0.6);

end

%% FIGURE with frames displayed
if 0
fClus=makegoodfig('figClustering','slide');
figs = vertcat(figs,fClus);
add_analysis_params(gcf,opt);

NFR = Vin.NumFrames;
dots = gobjects(NFR,1);

% ax1
AX(1) = subplot(1,2,1,'replace');
hold(AX(1),'on');

for ifr= 1:NFR
    this_fr_cluid = cluids(ifr);
    i4clu = pclu==this_fr_cluid;
    cluco = GrpCM(i4clu,:);

    dots(ifr)=plot(cluxy(ifr,1),cluxy(ifr,2),'.','color',cluco);
    menuu=uicontextmenu(fClus);
    uimenu(menuu,'Label',sprintf('show frame %g',ifr),'Callback',@FIG_showFrame,'UserData',ifr)
    uimenu(menuu,'Label','show all cluster frames','Callback',@FIG_showFrame,'UserData',ifr)
    uimenu(menuu,'Label','label cluster frames','Callback',@FIG_showFrame,'UserData',ifr)
    uimenu(menuu,'Label','label this frame as','Callback',@FIG_showFrame,'UserData',ifr)

    set(dots(ifr),'ContextMenu',menuu,'UserData',ifr,'ButtonDownFcn',@FIG_showFrame)

end
goodax(AX(1),'axis',{'equal'},'colormap',GrpCM,'colorbar',{'ylabel','cluster IDs'},'caxis',[0.5 nclu+0.5],'title',sprintf('nclu =  %g',nclu));

%% prep ax to show frame
AX(2) = subplot(1,2,2,'replace');
I = getframes(Vin,1);
if contains('ROI',opt.steps)
    I= imcrop(I,ROIpos);
end
imagesc(I);
goodax(AX(2),'axis',{'equal' 'tight'},'tick',[],'colormap','gray','color','none','box','on');
% if opt.ROI
%     rect= drawrectangle(AX(2),'position',ROIpos,'label',opt.ROI.name,'FaceAlpha',0);
% end
%% ----
%% FIND LABELS
addSpeedToLabeledImg= true;
if ismember('motion',opt.steps) && addSpeedToLabeledImg
    label_folder = fullfile(fig_folder,'frame labels',sprintf('motion_%s_%s_%s_speed',GCE.a,GCE.d,GCE.exp.vid));
elseif ismember('motion',opt.steps)
    label_folder = fullfile(fig_folder,'frame labels',sprintf('motion_%s_%s_%s',GCE.a,GCE.d,GCE.exp.vid));
else
label_folder = 'hello';
    %label_folder = fullfile(fig_folder,'frame labels',sprintf('pose_%s_%s_%s',GCE.a,GCE.d,GCE.exp.vid));
end

labels = cell(NFR,2);
labels(:,1)= cellfun(@(this_ifr) sprintf('%s_frame%g.png',myfilename(Vin.Name),this_ifr),num2cell(1:NFR),'uniformoutput',false);
if exist(label_folder,'dir')
    imgst = imageDatastore(label_folder, 'LabelSource', 'foldernames','IncludeSubfolders',true);

    if numel(imgst.Files)>0
        [filepaths,filenames,ext]= fileparts(imgst.Files);
        labeled_files = cellfun(@(x,y) [x y],filenames,ext,'uniformoutput',false);
        frame_ind = nan(numel(labeled_files),1);
        for ifr = 1:numel(labeled_files)
            labeled_filename = labeled_files{ifr};
            temp = strsplit(labeled_filename,{'frame' '.'});
            frame_ind(ifr)= str2num(temp{2});
        end
        labels(frame_ind,2) = cellstr(imgst.Labels);


        selected_dots = dots(frame_ind);
        set(selected_dots,'marker','x')


    end


    clear imgst;
end


%% --------
%% USER DATA
GCE.url_vid = url_vid;
cam_id = mycamID(url_vid);
fClus.UserData=var2struct({'GCE' 'Vin' 'AX' 'dots' 'pclu' 'nclu' 'cluids'...
    'GrpCM' 'cluxy' 'cam_id' 'opt'...
    'ROIpos'  'NFR' 'motionStep'  'labels','label_folder'});
fClus.KeyPressFcn = @FIG_KEYPRESS;
menuu2 = menuu.Children;
FIG_showFrame(menuu2,0);
end
fprintf('\nDONE !!')
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% FUNCTIONS

%% showFrame
function FIG_showFrame(SRC,eee)


CLASS = class(SRC);
Opt_ShowThisFrame = false;
Opt_ShowAllFrames = false;
Opt_LabelAllFrames = false;
Opt_LabelThisFrame = false;
switch CLASS
    case 'matlab.ui.container.Menu'
        mmenu = SRC;
        ifr = mmenu.UserData;
        f = ancestor(mmenu,'figure');
        try
            f = f{1};
        end
        if ~isa(eee,'double')
            Opt_ShowThisFrame = any(contains('show frame',{mmenu.Text}));
            Opt_ShowAllFrames = any(strcmp('show all cluster frames',{mmenu.Text}));
            Opt_LabelAllFrames = any(strcmp('label cluster frames',{mmenu.Text}));
            Opt_LabelThisFrame = any(strcmp('label this frame as',{mmenu.Text}));
        end

    case  'matlab.graphics.chart.primitive.Line'
        dot=SRC;
        ifr = dot.UserData;
        f = ancestor(dot,'figure');
    otherwise
        keyboard;
end
ANY_MENU_CLICK = any([Opt_ShowAllFrames Opt_LabelThisFrame Opt_LabelAllFrames Opt_ShowThisFrame]);
fud = f.UserData;
struct2var(fud,'fud');

if isempty(ifr)

    keyboard
end
dot = dots(ifr);
clu_id = cluids(ifr);
if Opt_ShowAllFrames || Opt_LabelAllFrames
    ifr = find(cluids==clu_id)';
    disp(['   >>found ' num2str(numel(ifr)) ' frames in cluster'])
end
url_vid = GCE.url_vid;
if Opt_LabelAllFrames  || Opt_LabelThisFrame
    label_list = { 'running' 'grooming' 'whisking' 'sniffing' 'reajusting' 'immobile' 'walking' 'other' 'keyboard' };
    label_name_idx = listdlg('ListString',label_list,'SelectionMode','single');
    if isempty(label_name_idx)
        return;
    end
    label_name = label_list{label_name_idx};


end
w = Vin.Width;
h = Vin.Height;

% make fig stop
fstop=make_stop_figure(f);
figure(f);

% change selected to different marker
if ANY_MENU_CLICK
    selected_dots = dots(ifr);
    if Opt_LabelAllFrames || Opt_LabelThisFrame
        set(selected_dots,'marker','x')
    elseif Opt_ShowAllFrames || Opt_ShowThisFrame
        set(selected_dots,'marker','+')
    end
end

% show frame on second ax
AX2=AX(2);
im=findobj(AX2,'type','image');
tic;
sz=0;
N=numel(ifr);
i=0;
nfr = Vin.NumFrames;

for this_ifr = ifr(:)'
    i = i+1;
    cleanline(sz);
    sz=fprintf('cluster %g :frame %g/%g (%.f%%)',clu_id,i,N,100*i/N);

    if Opt_LabelAllFrames || Opt_LabelThisFrame
        f.UserData.labels{this_ifr,2}=label_name;
        continue;

    end

    if motionStep
        if this_ifr~=NFR
            I = getframes(Vin,this_ifr+(0:1));
        else
            I = getframes(Vin,1:2);
        end
        I(:,:,3) = I(:,:,1);
    else
        I = getframes(Vin,this_ifr,1);
    end
    if contains('ROI',opt.steps)
        I= imcrop(I,fud.ROIpos);
    end

    im.CData = I;
    if numel(size(I))==3
    im.CDataMapping = 'direct';
    else
    im.CDataMapping = 'scaled';
    end

    if motionStep
        dI = I(:,:,1)-I(:,:,2);
        im.AlphaData=abs(dI)+0.3;

    end


    goodax(AX2,'xycolor',dot.Color,'linewidth',2);

    title(AX2,sprintf('cluster %g :frame %g/%g',clu_id,this_ifr,nfr),'color',dot.Color)
    pause(0.1);

    if toc>2
        drawnow;
        tic;
    end
    if ~isvalid(fstop)
        break;
    end
end

    nlabeled = sum(~cellfun('isempty',f.UserData.labels(:,2)));
    fprintf('\n ==> %.0f %% frames labeled (%g frames)',100*nlabeled/Vin.NumFrames,nlabeled);

end

%% play_frames
function play_frames(ax,Vin,ifr,Colors,varargin)
Vinfs = Vin.FrameRate;
im=findobj(ax,'type','image');%im=findobj(gca,'type','image')
im.CDataMapping = 'direct';
nfr = Vin.NumFrames;
if ~isempty(varargin)
    opt.ROI =1
    ROIpos = varargin{1};
else
    opt.ROI =0
end
if size(Colors,1)==1
    Colors = repmat(Colors,numel(ifr),1);
end
for this_ifr = ifr
    Vin.CurrentTime = (this_ifr-0.5)/Vin.FrameRate;
    thisFrame = readFrame(Vin);
    I = thisFrame;
    if opt.ROI
        I= imcrop(I,ROIpos);
    end
    im.CData = I;
    thisColor = Colors(this_ifr,:);
    goodax(ax,'xycolor',Colors(this_ifr,:),'linewidth',2);
    nfr = Vin.NumFrames;
    title(ax,sprintf('cluster %g: frame %g/%g',clu_id,this_ifr,nfr),'color',thisColor)
    pause(0.1);
end

end

%% montage_frames(ax,Vin,ifr,Colors)
function  montage_frames(ax,Vin,ifr,Colors,varargin)
nfr_montage = numel(ifr);
Vinfs = Vin.FrameRate;
im=findobj(ax,'type','image');
nfr = Vin.NumFrames;

% options
Opt_motion = 0;
opt.ROI =0;

if ~isempty(varargin)
    ROIpos = varargin{1};
    if ~isempty(ROIpos)
        opt.ROI =1;
    end
    if numel(varargin)>1
        Opt_motion = 1;
    end

end

if opt.ROI==0
    h = Vin.Height;
    w = Vin.Width;
else
    Vin.CurrentTime=0;
    I= Vin.readFrame();
    I= imcrop(I,ROIpos);
    h = size(I,1);
    w = size(I,2);
end

if size(Colors,1)==1;
    Colors = repmat(Colors,numel(ifr),1);
end

Iexamples= nan(h,w,3,nfr_montage);
tic

Margin = round(h/30);
Iexamples = nan(h,w,3,nfr_montage);
pos =0;
tic;
f = findobj('type','figure','name','figClustering');

fstop = make_stop_figure(f);


for this_ifr = ifr(:)'
    pos = pos+1;


    if Opt_motion==1
        if this_ifr==nfr
            thisFrame = getframes(Vin,this_ifr*[1 1 1]);
        else
            thisFrame = getframes(Vin,this_ifr+(0:1));
            thisFrame(:,:,3) = thisFrame(:,:,1);
            dI = abs(diff(thisFrame(:,:,1:2),1,3))<0.1;
            [irows,icols]=find(dI);
            i3 = repmat(1:3,numel(icols),1);
            dI3 = sub2ind([h w 3],repmat(irows,3,1),repmat(icols,3,1),i3(:)) ;
            thisFrame(dI3) = thisFrame(dI3)+0.5;
            thisFrame(thisFrame>1)=1;
        end
    else
        thisFrame = repmat(getframes(Vin,this_ifr),1,1,3);

        if pos>1
            thisFrame(:,:,1) = Iexamples(:,:,1,1);
        end
    end
    if max(thisFrame(:))<=1
        isRGB=1;
    else
        isRGB=0;
    end
    if opt.ROI
        thisFrame= imcrop(thisFrame,ROIpos);
    end

    thisColor = Colors(pos,:);
    if isRGB
        thisColor2 = permute(thisColor,[1 3 2]);
    else
        thisColor2 = permute(thisColor*255,[1 3 2]);
    end

    % make margins
    thisFrame(1:Margin,:,:)= repmat(thisColor2,Margin,w);
    thisFrame(:,1:Margin,:,:)= repmat(thisColor2,h,Margin);
    thisFrame(end-(Margin-1):end,:,:)= repmat(thisColor2,Margin,w);
    thisFrame(:,end-(Margin-1):end,:)= repmat(thisColor2,h,Margin);

    Iexamples(:,:,1:3,pos)= thisFrame;%figure;imagesc(thisFrame)

    if toc>2
        drawnow;
    end

    if ~isvalid(fstop)
        break;
    end

end

Iexamples= (Iexamples-mymin(Iexamples(:)))/(mymax(Iexamples(:))-mymin(Iexamples(:)));



axes(ax);
cla(ax);
ncol = 4;
nrow = nfr_montage/ncol;
montage(Iexamples,'Size',[nrow ncol]);
goodax(ax,'axis',{'equal'});

end
%% FIG_KEYPRESS
function FIG_KEYPRESS(f,EV)
U = f.UserData;
struct2var(U,'U');


switch EV.Key
    case {'c' 'q' 'r' 'd'}
        if ~isfield(f.UserData,'fluo')
            [FLUO,~,~,~,~ ]=evalin('base','T.call_proc_s2p_resp_cells_x1rec();');%[FLUO,FOV,EVENTS,PSTHs,LOC]=evalin('base','T.call_proc_s2p_resp_cells_x1rec();');
            f.UserData.fluo = FLUO;
            fluo = FLUO;

        end

        if EV.Key=='r' && ~isfield(f.UserData,'EVENTS')
            [~,~,EVENTS,~,~]=evalin('base','T.call_proc_s2p_resp_cells_x1rec();');
            f.UserData.EVENTS = EVENTS;
        end

        nce = fluo.cnt;
        tphy = fluo.time;%t(end)


        if ~isfield(f.UserData,'cam')
            continuous=evalin('base',"T.call_get_proc_data('continuous');");
            cam = continuous.cam;
            f.UserData.cam = cam;
        end
        cam = cam(cam_id);
        i4rec = strcmp(cam.filename,Vin.Name);
        tcam = cam.time(i4rec);

end

switch EV.Key



    case 'c'% cell
        disp '"c"'

        fstop = make_stop_figure(f);
        figure(f);
        tic;
        sz=fprintf('\n ');
        for cell_ind = 1:nce
            cleanline(sz);
            sz=fprintf('\n cell %g/%g ');

            cell_id = fluo.id(cell_ind,:);
            cellid_str = join(num2cellstr(cell_id),'_');
            phy = fluo.data.norm(cell_ind,:);
            phy_int = interp1(tphy,phy,tcam);%figure;plot(tphy,phy,'-b');hold on;plot(tcam,phy_int,'-r')


            % make figure
            fig(cell_ind)=makegoodfig(sprintf('cell_activity_xbehStates_m%s',cellid_str{1}),'slide');

            inonan = ~isnan(phy_int);
            X = cluxy(inonan,1);
            Y = cluxy(inonan,2);
            C = phy_int(inonan);
            nfr = sum(inonan);

            AXce(1) = subplot(1,2,1,'replace');
            scatter(X,Y,ones(nfr,1)*6,C,'filled');
            colormap('turbo');
            CLIM = [0 25];
            CLIM = myminmax(phy_int);
            goodax(gca,'axis',{'equal' 'tight'},'title',{sprintf('Cell %s',num2str(cell_id))},'colorbar',{'limits',CLIM,'caxis',CLIM,'ylabel','norm fluo'},'xlabel','t-SNE dim 1 (a.u.)','ylabel','t-SNE dim 2 (a.u.)');
            hold(gca,'on');

            % add outline of clusters
            outclu = gobjects(nclu,1);
            n4clu = nan(nclu,1);
            isvalidout = false(nclu,1);
            txtgrp = gobjects(nclu,1);
            for iclu = 1:nclu
                cluid = pclu(iclu);
                i4clu = cluids==cluid;
                n4clu(iclu)=sum(i4clu);

                if n4clu(iclu)==0
                    continue
                end
                xclu = cluxy(i4clu,1);
                yclu = cluxy(i4clu,2);
                i4out = boundary(cluxy(i4clu,:),0.1);
                if ~isempty(i4out)
                    I4clu = pclu==cluid;
                    xout = xclu(i4out);
                    yout = yclu(i4out);
                    outclu(iclu)=plot(xout,yout,'-','color',GrpCM(I4clu,:),'LineWidth',1);
                    isvalidout(iclu)=true;
                    [xtxt,i4y] = max(xout);
                    ytxt = yout(i4y);
                    txtgrp(iclu)=text(xtxt,ytxt,num2str(cluid),'color',GrpCM(I4clu,:),'HorizontalAlignment','left');
                end
            end
            leglabelid = num2cellstr(pclu(isvalidout));
            leglabeln = num2cellstr(n4clu(isvalidout));
            leglabel = cellfun(@(x,y) [x ' (' y ')'],leglabelid,leglabeln,'UniformOutput',false);
            legend(outclu(isvalidout),leglabel,'location','southoutside','NumColumns',6);

            AXce(2) = subplot(1,2,2,'replace');


            % sort by best clusters and select highest activtiy frames
            clu_ids = cluids(inonan);
            avg_fluo_xClu = groupsummary(table(C,clu_ids),'clu_ids','median');
            [Csort,best_clu_order]=sort(avg_fluo_xClu.median_C ,'descend');
            [Csort,isort]=sort(C ,'descend');
            selectedFrames = isort(1:24);
            selectedFrames_clu = clu_ids(selectedFrames);

            % sort by best cluster
            selectedFrames_sort = [];
            for ibe = 1:nclu
                clufr = selectedFrames_clu==best_clu_order(ibe);
                selectedFrames_sort = vertcat(selectedFrames_sort,selectedFrames(clufr));
            end

            selectedFrames = selectedFrames_sort;
            selectedFrames_clu = clu_ids(selectedFrames);

            % get colors sorted
            [~,selectedFrames_clu_pos] = ismember(selectedFrames_clu,pclu);
            Colors = GrpCM(selectedFrames_clu_pos,:);
            if motionStep
                if contains('ROI',opt.steps)
                    montage_frames(AXce(2),Vin,selectedFrames,Colors,ROIpos,'motion');
                else
                    montage_frames(AXce(2),Vin,selectedFrames,Colors,[],'motion');
                end
            else
                if contains('ROI',opt.steps)
                    montage_frames(AXce(2),Vin,selectedFrames,Colors,ROIpos);
                else
                    montage_frames(AXce(2),Vin,selectedFrames,Colors);
                end
            end
            goodax(AXce(2),'title',{{'Top 24 frames with cell highest activity';'sorted by descending activity cluster'}});
            cb=mycolorbar(AXce(2),GrpCM,pclu);
            ylabel(cb,'Cluster IDs');

            % do stats
            AXstats = axes();
            cla;
            p =AXce(2).Position;
            set(fig(cell_ind),'units','normalized');
            y0 =0.05;
            yspace = p(2)-y0;
            set(AXstats,'position',[p(1) y0 p(3) yspace*1.4],'colormap',GrpCM,'PositionConstraint','innerposition')
            % sort groups by descencing order of neuron's activity
            [~,descendOrder] = sort(avg_fluo_xClu.median_C,'descend');
            sortedGrp = avg_fluo_xClu.clu_ids(descendOrder);
            [isort,bx_id,bx_vals]=sortby(sortedGrp,clu_ids,C);
            bx_CM = GrpCM(descendOrder,:);
            bo=boxplot(bx_vals,num2cellstr(bx_id),'boxstyle','filled','plotstyle','compact','colors',bx_CM,'OutlierSize',0.5,'FactorDirection','data');

            Stats=KMStats(AXstats,1:nclu,C,clu_ids,'groupOrder',sortedGrp,'plotbars',false);
            goodax(AXstats,'box','off','xlabel',{'cluster id'},'ylabel','avg activity');
            set(AXstats,'position',[p(1) y0 p(3) yspace*1.5],'colormap',GrpCM,'PositionConstraint','innerposition')
            add_analysis_params(gcf,Opt);

            if toc>2
                drawnow;
            end

            if ~isvalid(fstop)
                break;
            end

            savethisfig(fig_folder,fig(cell_ind),Opt,fn0vid)



        end
    case 'q'
        %% --- "q" quantify number of clusters with neuronal activity different from others
        disp '"q"'

        tic;

        nce = fluo.cnt;

        % compute inonan
        PHY = fluo.data.norm;
        phy_int = interp1(tphy,PHY',tcam(tcam>0))';%figure;plot(tphy,phy,'-b');hold on;plot(tcam,phy_int,'-r')
        inonan = ~isnan(phy_int(1,:));
        nfr = sum(inonan);
        phy_int = phy_int(:,1:nfr);


        clu_ids = cluids(tcam>0);
        clu_ids = clu_ids(inonan);

        % compute number of frames in cluster
        n4clu = nan(nclu,1);
        i4clu = cell(nclu,1);
        for iclu = 1:nclu
            cluid = pclu(iclu);
            clu_idx = pclu==cluid;
            i4clu{iclu} = clu_ids==cluid;
            n4clu(iclu)=sum(i4clu{iclu});


        end

        % exclude small clusters
        i4smallgrp=find(n4clu<10);
        cluid2excl = pclu(i4smallgrp);
        i4clu2excl = ismember(clu_ids,cluid2excl);

        phy_int(:,i4clu2excl)=nan;

        %prealloacte
        clu_issig = false(nce,nclu);
        clu_vals = nan(nce,nclu);
        fprintf('\n Computing mean values across clusters for each cell ...')
        sz=0;
        for cell_ind = 1:nce
            cleanline(sz);
            sz=fprintf('\ncell %g/%g  %.0f%%',cell_ind,nce,100*cell_ind/nce);
            phy = phy_int(cell_ind,:);

            % compute stats anova
            Stats=KMStats([],1:nclu,phy(:),clu_ids);
            c = Stats.sigdata.comparisons;
            issig_xclu = false(nclu,1);

            for iclu=1:nclu

                % find if cluster is significantly different from others
                i4gr = any(c(:,1:2)==iclu,2);
                issig_xclu = c(i4gr,end);
                if any(issig_xclu<Stats.Alphas(1))
                    clu_issig(cell_ind,iclu)=true;
                end

                % compute mean fluo value for cluster
                clu_vals(cell_ind,iclu)=mymean(phy(i4clu{iclu}));

            end

        end%cell_ind

        % value should be above 2 STD
        nSTD = 2;
        thSTD = nSTD*mystd(clu_vals,2)+mymean(clu_vals,2);
        i4val4SigClu = clu_vals>=thSTD;
        i4diffClu = clu_issig & i4val4SigClu;

        n4SigClu = sum(i4diffClu,2);
        p4SigClu = n4SigClu/nclu;

        % plot quantification histogram across cells
        figNH = makegoodfig('figNeuHist','slide_half_width');
        his_edges = [0 1 2 3 4 nclu+1];
        his=histcounts(n4SigClu,his_edges,'Normalization','probability');

        ba = bar(1:numel(his),his);
        ax=gca;
        goodax(ax,'xtick',1:numel(his),'xticklabel',mynum2str(his_edges(1:end-1),[],'cellstr'),...
            'xlabel',sprintf('Nb of clusters with diff. neuronal activity\n(pval<%.3f,threshold = %g STDs, %g clusters)',Stats.Alphas(1),nSTD,nclu),...
            'ylabel','probability','axis','square','title',sprintf('N = %g cells',nce));
        ax.XTickLabel(end,:) = {'>4'};
        add_analysis_params(figNH,Opt);
        savethisfig(fig_folder,figNH,Opt,fn0vid)
    case "r"
        %% --- "r" quantify neuronal activity during sound vs cluster

        PHY = fluo.data.norm;
        nce = size(PHY,1);
        if nce==0
            fprintf('\n NO CELLS')
        end
        figRespxClu=fig_sound_resp_xClu_x1rec(tphy,PHY,tcam,cluids,EVENTS,Opt,cluxy);

        %% - save fig
        savethisfig(fig_folder,figRespxClu,Opt,fn0vid)
    case "d"
        %% --- "d" clusterless analysis dist vs activity diff
        disp(' >>> d : clusterless analysis dist vs activity diff');
        tic;
        th_prct = 90;

        nce = fluo.cnt;

        % compute inonan
        PHY = double(fluo.data.norm');
        inonan = tcam>0 & ~isnan(interp1(tphy,PHY(:,1),tcam));
        phy_int = interp1(tphy,PHY,tcam(inonan));%figure;plot(tphy,phy,'-b');hold on;plot(tcam,phy_int,'-r')
        clear PHY;

        %prealloacte

        fprintf('\n Computing distances and activity differences  ...')
        sz=0;
        PC_scores = load(urlstep.PCA,'score');
        PC_scores=PC_scores.score(inonan,:);
        dist_dim=pdist2(PC_scores,PC_scores,"squaredeuclidean");
        dist_dim = dist_dim./max(dist_dim(:));
        nfr = size(PC_scores,1);
        RHO=nan(nce,1);
        PVAL = nan(nce,1);
        DIFF_DIST =   nan(nce,1);
        PVAL_ttest =  nan(nce,1);
        nrow = 2;
        ncol = 1;
        XCLU = cluxy(inonan,1);
        YCLU = cluxy(inonan,2);
        ticPARFOR=tic;
        parfor ice=1:nce
            f4ce = phy_int(:,ice);%figure;plot(f4ce)
            fpairs = sum(cat(3,repmat(f4ce',nfr,1),repmat(f4ce,1,nfr)),3);
            %             all_d = all_d/mymax(all_d(:));
            df =fpairs(:);
            d = dist_dim(:);
            tic
            [RHO(ice),PVAL(ice)]=corr(discretize(d,100),discretize(df,100));
            toc

            % test if high activity frames have closer distances
            th_df = prctile(df,th_prct);
            md= mymean(d);
            i4hi = df>=th_df;
            i4lo = i4hi ==0;
            md_lo = mymean(d(i4lo));
            md_hi = mymean(d(i4hi));
            DIFF_DIST(ice) = md_lo-md_hi;
            tic
            [~,PVAL_ttest(ice)]=ttest2(d(i4lo)',d(i4hi)');
            toc
            %             tic
            %             if 0
            %                 ftemp=makegoodfig('cell_dist_example');
            %                 rightMargin = 0.7;
            %                 h=binScatterPlot(tall(d),tall(df),10);
            %                 ax2 = mysubplot(ancestor(h,'axes'),nrow,ncol,2,1,'rightMargin',rightMargin);
            %                 goodax(ax2,'xlabel','Distance (norm uMAP scores)','ylabel','Summed activity (norm)','title',sprintf('rho = %.2f, p = %.3f, diff_dist = %.2f',RHO(ice),PVAL(ice),DIFF_DIST(ice)));
            %                 hold(ax2,'on');
            %                 fi=my_linear_fit(ax2,d,df,1,0.05);
            %
            %
            %                 ax1 = mysubplot([],nrow,ncol,1,1,'rightMargin',rightMargin-0.1);
            %                 sc= scatter(XCLU,YCLU,ones(nfr,1)*6,phy_int(:,ice),'filled');
            %                 goodax(ax1,'ylabel',sprintf('%s 1 (a.u.)',opt.DimRedName),'xlabel',sprintf('%s 2 (a.u.)',opt.DimRedName),...
            %                     'colormap','turbo','colorbar',{'title','fluo'},'axis','equal',...
            %                     'title',sprintf('cell %g/%g',ice,nce));
            %                 pause();
            %             end
            %             toc
        end
        TOCparfor= toc(ticPARFOR)

        % plot quantification histogram across cells
        figNH = makegoodfig('figNeuDistHist','slide');
        nrow =1;
        ncol =2;
        AX = mysubplot([],nrow,ncol,[],[],'bottomMargin',0.2,'topMargin',0.8);

        % sub 1
        axes(AX(1));%cla;
        Alpha = 0.001/nfr^2;
        his_edges = -1:0.01:1;
        his=histogram(RHO,his_edges);
        his_sig=histogram(RHO(PVAL<=Alpha),his_edges);

        set(his,'FaceColor',  [0 0 1],'EdgeAlpha',0);
        set(his_sig,'FaceColor',  [1 0 0],'EdgeAlpha',0.5,'FaceAlpha',0.5,'EdgeColor',[0 0 0]);
        goodax(AX(1),'xlim',[-1 1]*mymax(abs(RHO)),'xlabel',{{'high activity clustered                       low activity clustered';'correlation coeff.'}},...
            'title',{{'Cells with correlation between' 'distance and summed activity ' mymeanstd(RHO,3)}});

        % sub 2
        axes(AX(2));%cla;

        his=histogram(DIFF_DIST,20);
        gid=discretize(DIFF_DIST,his.BinEdges);
        his_sig=histogram(DIFF_DIST(PVAL_ttest<=Alpha),his.BinEdges);
        set(his,'FaceColor',  [0 0 1],'EdgeAlpha',0);
        set(his_sig,'FaceColor',  [1 0 0],'EdgeAlpha',0.5,'FaceAlpha',0.5,'EdgeColor',[0 0 0]);
        goodax(AX(2),'xlim',1.05*[-1 1]*mymax(abs(DIFF_DIST)),'title',{{'Comparison high vs low activity distances' 'for each cell'   sprintf('th_Prctile = %g ; %s',th_prct,mymeanstd(DIFF_DIST,3))}},'xlabel',{{'low activity clustered                 high activity clustered';'dist(low activity)-dist(high activity)'}});

        add_analysis_params(figNH,Opt);
        savethisfig(fig_folder,figNH,Opt,fn0vid)

        %% ------  figure with all cells

        [fig_all,AXlist,i4row,i4col,xpos,ypos,ANN,i4fig,AXrcf] = preplotXfigs(10,9,nce);
        [~,cell_vec] = sort(DIFF_DIST,'descend');
        iax = 0;
        for ice=cell_vec(:)'
            iax = iax+1;
            ax = AXlist(iax);
            f4ce = phy_int(:,ice);%figure;plot(f4ce)

            szm = 4;
            sc=scatter(ax,XCLU,YCLU,ones(nfr,1)*6,f4ce,'marker','.');
            set(sc,'SizeData',ones(nfr,1)*szm,'AlphaData',1)
            h=goodax(ax,'colormap','turbo','axis','equal','xycolor','none',...
                'title',{{sprintf('cell %g/%g',ice,nce), sprintf('rho=%.2f,ddist= %.2f',RHO(ice),DIFF_DIST(ice))},'fontsize',8});;

        end
        CLIM = [0 10];
        for iax=1:nce
            ax = AXlist(cell_vec(iax));
            caxis(ax,CLIM)
        end
    case "l"% label frames
        %% --- "l" label frames and save them to folder
        speed = getSpeed(f);
        saveLabeledFrames(labels,label_folder,Vin,Opt,speed)

    case "n"% train neural network to classify images
        %% --- "n" train artificial neural network to classify labelled images'
        train_neural_network_to_classify_images


    case 'leftarrow'


    case 'rightarrow'


    otherwise
        disp 'HELP:'
        disp 'e: examples of frames for each cluster'
        disp 'c: show activity of cells for each frame'
        disp 'q: quantify number of clusters with neuronal activity different from others'
        disp 'r: quantify sound response during sound vs beh state cluster'
        disp 'l: label frames.'
        disp 'n: train artificial neural network to classify labelled images'
        disp 'left/arightarrow: show previous/next frames'



end

end

%% FUNCTION S

%% get_timestamps_frame_rate
function [t,fs]=get_timestamps_frame_rate(cam_id,Vin)
continuous=evalin('base',"T.call_get_proc_data('continuous');");
cam = continuous.cam;
cam = cam(cam_id);

fs = cam.frame_rate_comp(1);
i4rec = strcmp(cam.filename,Vin.Name);

t = cam.time(i4rec);
end

%% savethisfig
function savethisfig(fig_folder,fig,Opt,fn0vid)
%savethisfig(fig_folder,fig,Opt,fn0vid)
if evalin('base','T.el.option_fig.Value')

    OptFN = rmfield(Opt,'steps');
    params_str = get_params_string(OptFN);
    params_str = replace(params_str,'___','_');
    params_str = replace(params_str,'__','_');

    FN = [fig.Name '_' fn0vid '_' params_str];
    if numel(FN)>250
        keyboard;
    end

    % short name for illustrator
    %     stepnames = join(opt.steps,'_');
    %     urlfig = fullfile(fig_folder,[fig.Name '_' stepnames{1} '.emf']);
    %     mysaveas(fig(cell_ind),urlfig);
    %     %                 urlfig = fullfile(FO,[FN '.emf']);
    %                 mysaveas(fig(cell_ind),urlfig);
    urlfig = fullfile(fig_folder,[FN '.png']);
    mysaveas(fig,urlfig);

    disp(['*' FN '.png*'])% for ppt



end

end

%%  fstop = make_stop_figure()
function fstop = make_stop_figure(f)

fstop = makegoodfig('stop');
p = f.Position;
set(fstop,'Position' , [p(1)+p(3) p(2) 150 150],'color','k');
ann=annotation(fstop,'textbox');
set(ann,'color','r','fontweight','bold','String','CLOSE TO STOP','position',[0 0 1 1]);
end


%% saveLabeledFrames(labels,label_folder,Opt)
function saveLabeledFrames(labels,label_folder,Vin,Opt,all_speeds)

% make sure all label sub folders exist
plabels = unique(labels(~cellfun('isempty',labels(:,2)),2));
for ila = 1:numel(plabels)
    label = plabels{ila};
    mymkdir(fullfile(label_folder,label))
end


% save all labeled images that have not been saved already
NFR = size(labels,1);

speed_bins = -0.1:0.01:0.4;
all_speeds(all_speeds<min(speed_bins))=min(speed_bins);
all_speeds(all_speeds>max(speed_bins))=max(speed_bins);
speed_bins = 0:0.01:0.4;
speed_groups = discretize(all_speeds,speed_bins);
ngr = numel(speed_bins)-1;
nneg = sum(speed_bins<0.01);
negCM = cool(nneg*2);
CM = [negCM(nneg+1:end,:);jet(ngr-nneg)];
nrow = Vin.Height;
ncol = Vin.Width;
irows = round(0.97*nrow):nrow;
icols = 1:ncol;
i4speed = myindices(nrow ,ncol,irows,icols);

sz = 0;

for ifr = 1:NFR
    cleanline(sz);
    sz=fprintf('frame %g/%g (%.f%%)',ifr,NFR,100*ifr/NFR);

    label = labels{ifr,2};
    isLabeled = ~isempty(label);
    filename = labels{ifr,1};

    fo_img = fullfile(label_folder,label);
    mymkdir(fo_img);
    url_img = fullfile(label_folder,label,filename);
    isSaved = exist(url_img,'file')==2;
    
    if isLabeled% && ~isSaved%provisory
        IMG = getFrameImg(Vin,ifr,Opt);

        % add speed to image
        speed_gr = speed_groups(ifr);
        if isnan(speed_gr)
            co= [0.1 0.1 0.1];
        else
            co= CM(speed_gr,:);
        end
        IMG=mycolorpixels(IMG,i4speed,co);
       
        % write image
        imwrite(IMG,url_img);%winopen(url_img)
    end

end
nlabeled = sum(~cellfun('isempty',labels(:,2)));
fprintf('\n %.0f %% frames labeled',100*nlabeled/NFR);


end
%%   I = getFrameImg(Vin,ifr,Opt);
function         IMG = getFrameImg(Vin,this_ifr,Opt)

NFR = Vin.NumFrames;
motionStep = contains('motion',opt.steps);
if contains('ROI',opt.steps)
    f = findobj('figure','figClustering');
    fud = f.UserData;
    ROIpos = fud.ROIpos;
end

if motionStep
    if this_ifr~=NFR
        I = getframes(Vin,this_ifr+(0:1));
    else
        I = getframes(Vin,1:2);
    end
    I(:,:,3) = I(:,:,1);
else
    I = getframes(Vin,this_ifr);
end

if contains('ROI',opt.steps)
    I= imcrop(I,ROIpos);
end
IMG=I;
if motionStep
    IMG=I;

    dI = abs(I(:,:,1)-I(:,:,2));
    dI3 = repmat(dI,1,1,3);

    enhancing = 1.3;
    IMG(:,:,1)=I(:,:,1).*(1+enhancing*dI);
    IMG(:,:,2)=I(:,:,2).*(1+enhancing*dI);
    IMG=1-5*dI3.*(1-IMG);% transparency
    %             IMG = IMG-0.7*IMG(:,:,3);
    IMG(IMG>1)=1;

    %figure;
    %image(IMG);
end


end

%% getSpeed
function  speed = getSpeed(f);

if isfield(f.UserData,'speed')
    speed = f.UserData.speed;
else
    ud = f.UserData;
    GCE= ud.GCE;
    time_src = 'time';

    % get encoder 
    continuous=evalin('base',"T.call_get_proc_data('continuous');");

    enc = continuous.encoder;
    tenc = enc.(time_src)(:);
    enc_speed = enc.speed_ms(:);
    tbegs = unique(enc.tbeg);
    % get cam data
    cam = continuous.cam;
    icam= f.UserData.cam_id;
    tcam = cam(icam).(time_src)(:);

    % get cam data only for recorded session
    filenames = cam(icam).filename;
    filenames(cellfun('isempty',filenames))={'NA'};
    pfilenames = unique(filenames);
    
    [~,~,ext] = fileparts(filenames{1});
    filenames = cellfun(@(x) replace(x,ext,''),filenames,'UniformOutput',false);
    this_filename = sprintf('CAM%g_%s',icam,replace(GCE.fn0.vid,'CAM1_',''));
    i4filename =strcmp(filenames,this_filename);

    % get tcam for filenames
    tcam = tcam(i4filename);

    % check that number of frames matches timestamps
    ntcam = numel(tcam);
    NFR = ud.NFR;
    if ntcam>NFR
    tcam = tcam(1:NFR);
    else ntcam<NFR% more frames in video than timestamps
        nmoreframes = NFR-ntcam;
        dt = 1/cam(icam).frame_rate_comp(i4filename(1));
        tcam = [tcam; (tcam(end)+(1:1:nmoreframes)*dt)'];
    end

    % temporarily get rid of cam times without speed info
    i4cam = tcam>=tenc(1) & tcam<=tenc(end);
    tcam_temp = tcam(i4cam);
    tenc_temp = tenc;

    % interpolate speed data to match frames
    speed_temp = interp1(tenc_temp,enc_speed,tcam_temp,"makima");

    % pad with nans when speed unknown
    speed = nan(NFR,1);
    speed(ismember(tcam,tcam_temp))=speed_temp;

    % make sure speed does not have more frames than 


    % figure
     fspeed = makegoodfig('speed','slide_half_height');
     plot(tenc,enc_speed,'-k');
     hold on;
     plot(tcam,speed(:),'-r');

     f.UserData.speed = speed;
end


end


