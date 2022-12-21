function img = makeAllChannelBhvImage(gui,data,cmap,inds,idx_max,showAnnot)
%
% (C) Ann Kennedy, 2019
% California Institute of Technology
% Licensing: https://github.com/annkennedy/bento/blob/master/LICENSE.txt



chList = fieldnames(data);

if(~isempty(inds))
    numframe_win = length(inds);
    mask            = (inds<1)|(inds>idx_max);
    inds(inds<1)    = 1;
    inds(inds>idx_max) = idx_max;
else
    numframe_win = idx_max;
end
img = ones(length(chList),numframe_win,3);
for ich = 1:length(chList)
    ch = chList{ich};
    bhvList = fieldnames(data.(ch));
    
    for ibe = 1:length(bhvList)
        be = bhvList{ibe};

        show = (~isfield(showAnnot,be)) || (isfield(showAnnot,be) && showAnnot.(be));
        if(show && ~strcmpi(be,'other'))
            
            if(strcmpi(ch,gui.annot.activeCh))
                hits = gui.annot.bhv.(be)(inds);
                img(ich,hits~=0,:)   = ones(sum(hits),1)*cmap.(be);
            elseif(~isempty(data.(ch).(be)))
                bouts = data.(ch).(be);
                % make convertToRast faster
                if(~isempty(inds))
                    bouts(bouts(:,2)<inds(1),:) = [];
                    bouts(bouts(:,1)>inds(end),:) = [];
                    %bouts = bouts - inds(1);%original code
                    bouts = bouts - inds(1)+1;%KM fixed shift that was specific to not active channel
                end

                if~isempty(bouts)
                    hits                = convertToRast(bouts,numframe_win);
                    img(ich,hits~=0,:)   = ones(sum(hits),1)*cmap.(be);
                end
            end
        end
    end
    
end

if(~isempty(inds))
    img(:,mask,:)   = permute(ones(sum(mask),3),[3 1 2]);
    img = flip(img,1);
end
