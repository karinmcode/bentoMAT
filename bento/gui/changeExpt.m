function changeExpt(source,~)
%
% (C) Ann Kennedy, 2019
% California Institute of Technology
% Licensing: https://github.com/annkennedy/bento/blob/master/LICENSE.txt


useSource = source.Parent.Parent.Parent;
gui = guidata(useSource);

mList    = get(gui.ctrl.expt.mouse,'String');
sessList = get(gui.ctrl.expt.session,'String');
trList   = get(gui.ctrl.expt.trial,'String');

m       = str2num(mList{get(gui.ctrl.expt.mouse,'Value')});
sess    = str2num(sessList{get(gui.ctrl.expt.session,'Value')});
tr      = str2num(trList{get(gui.ctrl.expt.trial,'Value')});

% make sure we have a valid session for this mouse
matches = gui.allPopulated(gui.allPopulated(:,1)==m,:);
if(~any(matches(:,2)==sess))
    sess = matches(1,2);
    set(gui.ctrl.expt.session,'Value',1);
end

% make sure we have a valid trial for this session
matches = matches(matches(:,2)==sess,:);
if(~any(matches(:,3)==tr))
    tr = matches(1,3);
    set(gui.ctrl.expt.trial,'Value',1);
end

sess = ['session' num2str(sess)];
    

% determine whether we have to load a new seq movie (loading can take a while)
if(~gui.enabled.movie(1)) % movies aren't enabled
    newMovie = false;
    
elseif(~isfield(gui.allData(m).(sess)(tr).io.movie,'fid')) % new trial doesn't have a movie
    newMovie = false;
    set(gui.movie.img,'cdata',32*ones(1,1,3)); %blank screen
    
else % new trial has a movie- now check if it's the same
    if(~strcmpi(source.Tag,'trial')) %diffent mouse or session -> always assume a new movie
        newMovie = true;
    else
        len = [length(gui.data.io.movie.fid) length(gui.allData(m).(sess)(tr).io.movie.fid)];
        newMovie = len(1)~=len(2);  %check if the number of movies is the same
        for i = 1:min(len)          %check if the identity of each movie is the same
            newMovie = newMovie | ~strcmpi(gui.data.io.movie.fid{i}, ...
                                           gui.allData(m).(sess)(tr).io.movie.fid{i});
        end
    end
end
gui = setActiveMouse(gui,m,sess,tr,newMovie);
gui = redrawPanels(gui);

guidata(useSource,gui);
dummy.Source.Tag = 'slider';
gui.ctrl.slider.updateSliderAnnot(gui);

updatePlot(gui.h0,dummy);