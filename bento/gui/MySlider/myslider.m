classdef myslider<handle

    properties
        panel;
        Value;
        Min;
        Max;
        Scale;
        SliderStep;
        Pos;
        arrowW;
        textW;
        timer;
        LeftArrow;
        RightArrow;
        Ctr;
        bg;
        Marker;
        text;
        units;
        period;
    end

    methods

        %% s=myslider(gui)
        function s=myslider(gui,row,scale,nRows)

            s.init_slider(gui,row,scale,nRows);

        end

        %% init_slider(s)
        function s=init_slider(s,gui,row,scale,nRows)


            s.panel = uipanel('parent',gui.ctrl.panel,...
                'position',[0.01 (row-.5)/(nRows+1) 0.98 scale/(nRows+1)],'bordertype','none');

            s.Value      = 1;
            s.Min        = 1;
            s.Max        = 2;
            s.Scale      = 1-.015; %extra term corrects for marker width
            s.SliderStep = [1 1];
            s.Pos        = gui.ctrl.panel.Position .* s.panel.Position;
            s.arrowW     = 0.025;
            s.textW      = 0.125;
            s.timer      = tic;
            try
                reader = s.getActiveVideoReader(gui);
                s.period     = 1/reader.FrameRate;% time between frames when video is played, can be modified with f7 and f9
            catch
                s.period     = 0.1;
            end
            aW = s.arrowW;
            tW = s.textW;

            %make the slider
            s.LeftArrow = uipanel('Position',[0 0.05 aW 0.95],...
                'parent',s.panel,...
                'backgroundcolor',[.92 .92 .92],...
                'bordertype','line',...
                'HighlightColor',[.5 .5 .5],...
                'tag','leftArrow',...
                'buttondownfcn',{@s.sliderCheck,gui});
            makeButton(s.LeftArrow,char(9664));


            s.RightArrow = uipanel('Position',[1-aW-tW 0.05 aW 0.95],...
                'parent',s.panel,...
                'backgroundcolor',[.92 .92 .92],...
                'bordertype','line',...
                'HighlightColor',[.5 .5 .5],...
                'tag','rightArrow',...
                'buttondownfcn',{@s.sliderCheck,gui});%set(s.Ctr,'buttondownfcn',{@s.sliderCheck,gui,evt})
            makeButton(s.RightArrow,char(9656));

            %% bar Ctr
            s.Ctr        = uipanel('position',[aW 0.05 1-2*aW-tW 0.95],...
                'parent',s.panel,...
                'bordertype','line',...
                'HighlightColor',[.5 .5 .5],...
                'tag','bar',...
                'ButtonDownFcn',{@s.sliderCheck,gui});%set(s.Ctr,'buttonupfcn',{@s.sliderCheck,gui})
            axes('parent',s.Ctr,'position',[.015/2 0 1-.015 1]);
            s.bg    = imagesc([],'hittest','off'); axis tight; axis off;

            s.addContextMenu(s.Ctr,gui);


            %% Marker
            s.Marker     = uipanel('position',[s.Value*s.Scale 0.05 .015 0.9],...
                'Parent',s.Ctr,...
                'bordertype','line',...
                'HighlightColor',.7*ones(1,3),...
                'tag','marker',...
                'buttondownfcn',{@s.sliderCheck,gui});
            makeButton(s.Marker,'');


            %make the text box
            s.text   = uicontrol('Style','edit',...
                'String','',...
                'parent',s.panel,...
                'units','normalized',...
                'position',[1-tW 0.05 .95*tW*.8 0.95],...
                'Tag','timeBox',...
                'Callback',@updatePlot);
            s.units  = uicontrol('Style','pushbutton',...
                'string',char(8644),...
                'parent',s.panel',...
                'units','normalized',...
                'position',[1-tW*.2 0.075 .95*tW*.2 0.9],...
                'Tag','toggleTime',...
                'TooltipString','Switch to frame # (behavior video)',...
                'Callback',@changeTimeUnits);

        end

        %% gui=sliderCheck(s,gui,evt)
        function gui=sliderCheck(s,src,evt,gui)
            % gui=sliderCheck(s,gui,evt)

            s.timer = tic;
            switch evt.Source.Tag
                case 'leftArrow'
                    gui.Action = [-s.SliderStep(1) 1];   %left arrow
                case 'rightArrow'
                    gui.Action = [s.SliderStep(1) 1];    %right arrow
                case 'bar'
                    pClick  = get(0,'pointerlocation');

                    tracker = getpixelposition(s.Marker,true);
                    tracker(1) = gui.h0.Position(1) + tracker(1);

                    if(pClick(1) < tracker(1))
                        gui.Action = [-s.SliderStep(2) 1];   %track left
                    elseif(pClick(1) > tracker(1)+tracker(3))
                        gui.Action = [s.SliderStep(2) 1];    %track right
                    end
                    gui=s.incSliderVal(gui);

                case 'marker'
                    gui.Action = 'dragslider';
            end

        end

        %% gui=sliderDownCheck(s,gui,evt)
        function gui=sliderDownCheck(s,src,evt,gui)
            % gui=sliderDownCheck(s,gui,evt)

            s.timerDown = tic;
            switch evt.Source.Tag
                case 'leftArrow'
                    gui.Action = [-s.SliderStep(1) 1];   %left arrow
                case 'rightArrow'
                    gui.Action = [s.SliderStep(1) 1];    %right arrow
                case 'bar'
                    pClick  = get(0,'pointerlocation');

                    tracker = getpixelposition(s.Marker,true);
                    tracker(1) = gui.h0.Position(1) + tracker(1);

                    if(pClick(1) < tracker(1))
                        gui.Action = [-s.SliderStep(2) 1];   %track left
                    elseif(pClick(1) > tracker(1)+tracker(3))
                        gui.Action = [s.SliderStep(2) 1];    %track right
                    end
                    gui=s.incSliderVal(gui);

                case 'marker'
                    gui.Action = 'dragslider';
            end

        end

        %% gui=incSliderVal(s,gui)
        function gui=incSliderVal(s,gui)

            newValue=s.Value + gui.Action(1);

            % check if out-of-range
            if(((s.Value + gui.Action(1))>s.Max)||((s.Value + gui.Action(1))<s.Min))
                return;
            end

            % check for overshoot
            p       = get(0,'PointerLocation');         %click location
            marker  = getpixelposition(s.Marker,true);
            marker  = gui.h0.Position(1) + marker(1) + marker(3)/2; %marker coordinates
            try
                if(gui.Action(2) && (gui.h0==gcf) && (sign(p(1) - marker) ~= sign(gui.Action(1))))
                    return;
                end
            end

            s.Value = newValue;

            pos = (s.Value - s.Min)/(s.Max-s.Min)*s.Scale;
            if(isempty(pos))
                pos = 0;
            end
            s.Marker.Position(1) = pos;
            gui = gui;
            gui.ctrl.slider = s;
            eventdata.Source.Tag='slider';
            updatePlot(gui.h0,eventdata);
        end

        %% gui=setSliderVal(s,gui)
        function gui=setSliderVal(s,gui)
            p       = get(0,'PointerLocation');         % absolute location of the click
            figPos  = getpixelposition(gui.h0);         % absolute location of the figure

            
            parent = gui.Action(5:end);
            if isempty(parent)
                parent = 'slider';
            end
            if(strcmpi(parent,'slider'))
                obj     = gui.ctrl.slider.Ctr;
                scale   = (gui.ctrl.slider.Max - gui.ctrl.slider.Min);
                offset  = gui.ctrl.slider.Min - gui.ctrl.slider.Marker.Position(3)/2;
            else
                if(~strcmpi(parent,'features'))
                    obj = gui.(parent).axes;
                else
                    obj = gui.features.feat(1).axes;
                end
                scale   = -(obj.XLim(2) - obj.XLim(1));
                offset  = gui.(parent).clickPt(1);
            end

            dragTo    = getpixelposition(obj,true);     % location of the relevant box within the figure
            dragTo(1) = dragTo(1) + figPos(1);
            dragTo    = (p(1) - dragTo(1))/dragTo(3);

            s.Value = dragTo*scale + offset;
            s.Value = min(max(s.Value,s.Min),s.Max);

        end


        %% updateSliderAnnot(s,gui)
        function updateSliderAnnot(s,gui)
            % (C) Ann Kennedy, 2019
            % California Institute of Technology
            % Licensing: https://github.com/annkennedy/bento/blob/master/LICENSE.txt



            if(gui.enabled.annot(2))
                if(~isempty(gui.data.annoTime))
                    time 	= floor((gui.data.annoTime(end)-gui.data.annoTime(1))*gui.data.annoFR);
                elseif(all(gui.enabled.traces))
                    time = round(gui.data.CaTime(end) - gui.data.CaTime(1)*gui.data.CaFR);
                elseif(isfield(gui.data,'trackTime'))
                    time = floor(gui.data.trackTime(end) - gui.data.trackTime(1)/(gui.data.trackTime(2)-gui.data.trackTime(1)));
                else
                    time = 0;%KM
                end

                try
                    offset = gui.data.io.annot.tmin/gui.data.annoFR;
                catch
                    offset = 0;%KM
                end
                inds = find((gui.data.annoTime>(s.Min-offset)) & (gui.data.annoTime<=(s.Max-offset)));

                img     = makeBhvImage(gui.annot.bhv,gui.annot.cmap,inds,time,gui.annot.show);

                bgsmall = displayImage(img,s.Pos(3)*gui.h0.Position(3)*5,1);
                set(s.bg,'cdata',bgsmall);
            else
                set(s.bg,'cdata',[]);
            end
        end


        %%  gui = updateSlider(s,gui)
        function updateSlider(s)
            % updateSlider(s)
            pos = (s.Value - s.Min)/(s.Max-s.Min)*s.Scale;
            if(isempty(pos))
                pos = 0;
            end
            s.Marker.Position(1) = pos;



        end


        %% gui=applySliderUpdates(s,gui,type,info)
        function gui=applySliderUpdates(s,gui,type,info)
            %
            % (C) Ann Kennedy, 2019
            % California Institute of Technology
            % Licensing: https://github.com/annkennedy/bento/blob/master/LICENSE.txt



            % remove existing annotation data:

            set(s.bg,'CData',[]);

            switch type
                case 'movie'
                    if(strcmpi(info.readertype{1,1},'seq'))
                        Fr = info.FR; % use the experimenter-set value
                        tMax = info.reader{1,1}.numFrames/Fr;
                    else

                        try
                            Fr = info.FR;
                            tMax = info.tmax/Fr;
                        catch%KM
                            if numel(info.reader)==1
                            Fr = info.reader{1,1}.reader.FrameRate;%provisory
                            info.FR = Fr;
                            tMax = info.reader{1,1}.reader.NFrames;%provisory
                            else
                                keyboard;
                            end
                        end
                    end
                    maxVal = min(tMax,info.tmax/info.FR);
                    minVal = max(1/Fr,info.tmin/info.FR);
                case 'Ca'
                    minVal  = 0;
                    maxVal  = info.annoTime(end);
                    Fr      = 1/(info.annoTime(2)-info.annoTime(1));
                case 'audio'
                    minVal  = 0;
                    maxVal  = info.audio.t(1,end);
                    Fr      = 1/(info.audio.t(1,2)-info.audio.t(1,1));
                case 'tracker'
                    minVal = 0;
                    maxVal = info.trackTime(end);
                    Fr     = 1/(info.trackTime(2)-info.trackTime(1));
                otherwise % annotations really shouldn't be the last choice- it should
                    % override all other fields for setting the slider bounds,
                    % since annotations are what's displayed on the slider.
                    Fr      = info.annoFR;
                    minVal  = 0;
                    maxVal  = info.annoTime(end) - info.annoTime(1);
            end

            % change the slider limits and resets it to start
            s.Max         = maxVal;
            s.Min         = minVal;
            s.SliderStep  = [1/Fr 1];
            s.Value       = s.Min;
            s.updateSlider();

            gui.ctrl.slider = s;

        end

        %% playVideo(s,gui)
        function   playVideo(s,gui)
            ftemp = findobj('name','vidplayer','type','figure');
            ISPLAYING = isvalid(ftemp);
            if ISPLAYING
                close(ftemp)
                return;
            end

            ftemp = makegoodfig('vidplayer');
            set(ftemp,'color','c','MenuBar','none','position',[gui.h0.Position(1)+gui.h0.Position(3) gui.h0.Position(2) 144 65])
            drawnow;

            V=s.getActiveVideoReader(gui);%help VideoReader frame idx 1 inlcudes time == Period 
           
            Period = 1/V.FrameRate;
            gui.Action(1) = Period;
            % set slider value to currentTime
            s.Value= gui.data.annoTime(find(gui.data.annoTime<V.CurrentTime ,1,'last'));
            try
                if ftemp.Color~='g'
                    set(ftemp,'color','g');
                end
            end
            figure(gui.h0);

            while isvalid(ftemp)
                t0 = tic();
                gui = s.incSliderVal(gui);

                while toc(t0)<=s.period
                    drawnow;
                end

            end
            try
            close(ftemp);
            end


        end

        %%  playSelectBhv(s,src,evt,gui)
        function playSelectBhv(s,gui,varargin)
            if isempty(varargin)
                selectBhv=  {gui.annot.activeBeh};
            else
                selectBhv = varargin{1};
            end
            fprintf('\nSelected behaviors : %s' ,strjoin(selectBhv,', '))
          

            annoTime = gui.data.annoTime;
            if numel(selectBhv)==1
                selectBhv = selectBhv{1};
                BehTimes = annoTime(gui.annot.bhv.(selectBhv));
            else
                ns = numel(selectBhv);
                rast = zeros(size(annoTime));
                for i = 1:ns
                    rast = rast | gui.annot.bhv.(selectBhv{i});
                end
                BehTimes = annoTime(rast);
            end

            if isempty(BehTimes)
                fprintf('  >> No frames labelled.')
                return;
            end

            % make figure
            ftemp = findobj('name','vidplayer','type','figure');
            if isvalid(ftemp)
                close(ftemp)
            elseif isempty(ftemp)
                ftemp = makegoodfig('vidplayer');set(ftemp,'color','c','MenuBar','none','position',[2393 269 144 65]);
            else
                ftemp = makegoodfig('vidplayer');set(ftemp,'color','c','MenuBar','none','position',[2393 269 144 65]);
            end
            drawnow;

            % show frames of select behavior
            cnt = 0 ;
            Nfr = numel(BehTimes);
            while isvalid(ftemp) && cnt<Nfr
                cnt =cnt+1;
                Time = BehTimes(cnt);
                s.MoveToTime(gui,Time);
                try
                    if ftemp.Color~='g'
                        set(ftemp,'color','g');
                    end
                end
                pause(s.period);
                drawnow;
                figure(gui.h0);
            end
            try
            close(ftemp);
            end

        end

        %% addContextMenu
        function addContextMenu(s,hobject,gui)
            menuu=get(s.Ctr,'UIContextMenu');
            if isempty(menuu)
                if strcmp(class(hobject),'matlab.ui.Figure')
                    menuu = uicontextmenu(hobject);
                else
                    hfig = ancestor(hobject,'figure');

                    if strcmp(class(hfig),'matlab.ui.Figure')
                        menuu = uicontextmenu(hfig);
                    else
                        keyboard
                    end
                end

            end

            % display when right click
            uimenu(menuu,'Label','move here','Callback',{@s.MoveHere,gui},'UserData','m');
            % add right click menu
            set(hobject,'UIContextMenu',menuu)
        end

        %%  MoveHere(s,src,evt,gui)
        function MoveHere(s,src,evt,gui)
            % get mouse click position
            gui = guidata(gui.h0);
            gui= s.setSliderVal(gui);
            s.updateSlider();
            eventdata.Source.Tag = 'dummy';
            updatePlot(gui.h0,eventdata);
           
        end

        function MoveToTime(s,gui,Time)

            s.Value = Time;
            s.Value = min(max(s.Value,s.Min),s.Max);

            s.updateSlider();
            eventdata.Source.Tag = 'dummy';
            updatePlot(gui.h0,eventdata);
        end

        function reader = getActiveVideoReader(s,gui)
            if ~isfield(gui,'data')
                gui = guidata(gui.h0);
            end
            
            if numel(gui.data.io.movie.reader)==1
                reader =gui.data.io.movie.reader{1, 1}.reader;
            else
                keyboard;
            end
        end
    end





end