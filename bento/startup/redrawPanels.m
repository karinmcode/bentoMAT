function gui = redrawPanels(gui)
% this gets called when the user makes a change to the data being viewed-
% eg, opens a movie, adds tracking data, etc. It is terrible.
%
% (C) Ann Kennedy, 2019
% California Institute of Technology
% Licensing: https://github.com/annkennedy/bento/blob/master/LICENSE.txt


enabled     = fieldnames(gui.enabled);
config      = gui.config;

ctrlSize = ones(size(gui.config.ctrl));
for i = 1:length(gui.config.ctrl)
    if(isfield(gui.config.ctrlSc,gui.config.ctrl{i}))
        ctrlSize(i) = gui.config.ctrlSc.(gui.config.ctrl{i});
    end
end
ctrlSize = sum(ctrlSize)*gui.config.rowscale;

%% figure out which panels are enabled, and set their visibility
for i = 1:length(enabled)
    if(gui.enabled.(enabled{i})(2))
        gui.(enabled{i}).panel.Visible = 'on';
    else
        gui.(enabled{i}).panel.Visible = 'off';
    end
end
if(gui.enabled.tracker(2)) %one exception- we still have the movie panel on if just tracking is enabled
    gui.movie.panel.Visible = 'on';
end

%% determine whether our display is one column or two
leftOn  = gui.enabled.movie(2) || gui.enabled.tracker(2);
rightOn = gui.enabled.traces(2) || gui.enabled.features(2) || gui.enabled.scatter(2) ||...
            (~gui.enabled.movie(2) && (gui.enabled.audio(2) || gui.enabled.fineAnnot(2)));
numOn   = (gui.enabled.traces(2)||gui.enabled.scatter(2)) + ...
           gui.enabled.features(2) + ...
          (gui.enabled.audio(2)||gui.enabled.fineAnnot(2));

if(leftOn && rightOn)
    leftWidth = config.midline;
elseif(rightOn)
    leftWidth = 0;
else
    leftWidth = 1;
end
rightStart = leftWidth;
rightWidth = 1-leftWidth;
legendSize = .125;

%time for the world's shittiest panel layout code!
bump = 0;
if(leftOn)
    
    if(gui.enabled.ctrl(2)) % control panel at the bottom
        gui.ctrl.panel.Position = [0 0 leftWidth ctrlSize];
        bump = ctrlSize;
    end
    
    if(gui.enabled.audio(2)) %then audio/fineAnnot
        gui.audio.panel.Position = [0 bump leftWidth .15];
        gui.fineAnnot.panel.Visible='off'; % merge fineAnnot with audio if both are turned on
        bump=bump+.15;
    elseif(gui.enabled.fineAnnot(2))
        gui.fineAnnot.panel.Position = [0 bump leftWidth .15];
        bump=bump+.15;
    end
    
    %remaining vertical space goes to the movie
    if(gui.enabled.legend(2)) %show a legend of annotations
        
        gui.legend.panel.Position = [0 bump legendSize 1-bump];
        gui.movie.panel.Position  = [legendSize bump leftWidth-legendSize 1-bump];
    else % just show movie
        gui.movie.panel.Position  = [0 bump leftWidth 1-bump];
    end
    
end

bump = 0;
if(rightOn)
    % if left column wasn't displayed, need to place ctrl/audio/fineAnnot
    if(~leftOn && (gui.enabled.ctrl(2)||gui.enabled.features(2)||gui.enabled.audio(2)||gui.enabled.fineAnnot(2)))
        
        if(gui.enabled.ctrl(2)) %control panel first
            gui.ctrl.panel.Position = [0 0 1 ctrlSize];
            bump = ctrlSize;
        end
        
        lBump = 0; %squeeze in the legend panel if it hasn't been placed yet
        if(gui.enabled.legend(2))
            lBump = 1/4;
            gui.legend.panel.Position = [0 bump 1/4 1-bump];
        end
        
        if(gui.enabled.audio(2)) %then audio/fineAnnot if applicable
            gui.audio.panel.Position = [lBump bump 1-lBump (.95-bump)/numOn];
            gui.fineAnnot.panel.Visible='off';
            bump = bump + (.95-bump)/numOn;
        elseif(gui.enabled.fineAnnot(2))
            gui.fineAnnot.panel.Position = [lBump bump 1-lBump .3];
            bump = bump+.3;
        end
        
        if(gui.enabled.scatter(2) || gui.enabled.features(2) || gui.enabled.traces(2)) %then divide remaining space between features + traces
            if(gui.enabled.scatter(2)), use = 'scatter'; %don't allow both traces and scatterplot to be on, for now
            elseif(gui.enabled.traces(2)), use = 'traces';
            else use = [];
            end
            if(gui.enabled.features(2))
                gui.features.panel.Position = [lBump bump 1-lBump (1-bump)/(1+~isempty(use))];
                if(~isempty(use))
                    gui.(use).panel.Position   = [lBump bump+(1-bump)/2 1-lBump (1-bump)/2];
                end
            else
                gui.(use).panel.Position   = [lBump bump 1-lBump 1-bump];
            end
        else
            gui.features.panel.Position     = [lBump bump 1-lBump 1-bump];
        end
        
    else % don't have to worry about ctrl/audio/fineAnnot placement (because they're in the left column)
        
        if(gui.enabled.scatter(2) || gui.enabled.traces(2))
            if(gui.enabled.scatter(2)), use = 'scatter';
            else, use = 'traces'; end
            if(gui.enabled.features(2))
                gui.(use).panel.Position = [rightStart .5 rightWidth .5];
                gui.features.panel.Position = [rightStart 0 rightWidth .5];
            else
                gui.(use).panel.Position = [rightStart 0 rightWidth 1];
            end
        else
            gui.features.panel.Position = [rightStart 0 rightWidth 1];
        end
        
    end
end
