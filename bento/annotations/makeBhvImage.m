function img = makeBhvImage(bhv,cmap,inds,max_index,showAnnot)
%
% (C) Ann Kennedy, 2019
% California Institute of Technology
% Licensing: https://github.com/annkennedy/bento/blob/master/LICENSE.txt



bhvList = fieldnames(bhv);

if(~isempty(inds))
    L = length(inds);
    mask            = (inds<1)|(inds>max_index);
    inds(inds<1)    = 1;
    inds(inds>max_index) = max_index;
else
    L = max_index;
end
img = ones(1,L,3);
for i = 1:length(bhvList)
    b = bhvList{i};

    if((isempty(showAnnot) || ~isfield(showAnnot,b) || showAnnot.(b)) && ~isempty(bhv.(b)) && ~strcmpi(b,'other'))

        if(min(size(bhv.(b)))==2 || length(bhv.(b))==2)
            if(~isempty(inds))
                rast = bhv.(b);
                rast(rast(:,2)<inds(1),:) = [];
                rast(rast(:,1)>inds(end),:) = [];
                rast = rast - inds(1);
                hits = convertToRast(rast,L);
            else
                hits = convertToRast(bhv.(b),max_index);
            end
        else
            if(~isempty(inds))
                hits = bhv.(b)(inds);
            else
                hits = bhv.(b);
            end
        end
        img(1,hits~=0,:)    = ones(sum(hits),1)*cmap.(b);

    end
end

if(~isempty(inds)) %crop image for traces/tracker browsers
    img(1,mask,:)   = permute(ones(sum(mask),3),[3 1 2]);
end