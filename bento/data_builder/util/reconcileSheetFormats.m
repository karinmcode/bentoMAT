function [M,matches,fields] = reconcileSheetFormats(gui,raw)
%
% (C) Ann Kennedy, 2019
% California Institute of Technology
% Licensing: https://github.com/annkennedy/bento/blob/master/LICENSE.txt

% get rid of blank columns at the end of the columns %KM
params = raw(1,:); % list of data fields in the excel sheet
headers = raw(2,:); % list of data fields in the excel sheet
isBlankCol = false(numel(headers),1);
for i=1:numel(headers)
    try
        isBlankCol(i)=isnan(headers{i}) ;
    end
end
raw = raw(:,~isBlankCol);

fieldset = raw(2,:); % list of data fields in the excel sheet
fieldset =fieldset(~cellfun('isempty',fieldset)); % list of data fields in the excel sheet

fields = struct();
for i=1:length(fieldset)
    str = strrep(fieldset{i},' ','_');
    str = strrep(str,'_#','');
    fields.(str) = i;
end
fset = fieldnames(fields);

if(~isempty(gui))
    matchset = get(gui.t,'columnname');
else %hacks~
    matchset = {'Mouse','Sessn','Trial','Stim','Calcium imaging file','Start Ca',...
                'Stop Ca','FR Ca','Alignments','Annotation file','Start Anno','Stop Anno',...
                'FR Anno','Offset','Behavior movie','Tracking','Audio file','tSNE'};
end
matchset = strrep(matchset,' ','_');
matchset = strrep(matchset,'_#','');

matches = struct();
for i=1:length(matchset)
    str = strrep(matchset{i},' ','_');
    str = strrep(str,'_#','');
    matches.(str) = i;
end

nodata = find(isnan([raw{3:end,1}]))+2;
raw(nodata,:) = [];

[startFlag,stopFlag,FRFlag]=deal(0);
for f = 1:length(fset)
    ind = find(~cellfun(@isempty,strfind(lower(matchset),lower(fset{f}))));
    if(length(ind)>1)
        switch fset{f}
            case 'Start'
                if(startFlag)
                    ind = find(~cellfun(@isempty,strfind(matchset,'Start_Anno')));
                else
                    ind = find(~cellfun(@isempty,strfind(matchset,'Start_Ca')));
                    startFlag = 1;
                end
            case 'Stop'
                if(stopFlag)
                    ind = find(~cellfun(@isempty,strfind(matchset,'Stop_Anno')));
                else
                    ind = find(~cellfun(@isempty,strfind(matchset,'Stop_Ca')));
                    stopFlag = 1;
                end
            case 'FR'
                if(FRFlag)
                    ind = find(~cellfun(@isempty,strfind(matchset,'FR_Anno')));
                else
                    ind = find(~cellfun(@isempty,strfind(matchset,'FR_Ca')));
                    FRFlag = 1;
                end
        end
    end
    M(:,ind(1)) = raw(3:end,f);
end

%append blank columns for unreported fields
for i = 1:(length(matchset)-size(M,2))
    M(:,end+1)={[]};
end
    