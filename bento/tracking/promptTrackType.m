function trackType = promptTrackType(custom_list, filename)
%
% (C) Ann Kennedy, 2019
% California Institute of Technology
% Licensing: https://github.com/annkennedy/bento/blob/master/LICENSE.txt



pth         = fileparts(which('promptTrackType'));
trackerList = dir([pth filesep 'tracking unpackers' filesep '*.m']);
trackerList = {trackerList.name};
trackerList = strrep(trackerList,'.m','');

if exist('custom_list','var')
    useList = custom_list;
else
    useList = trackerList;
end

if ~exist('filename','var')
    filename = '';
end

if numel(useList)==1
    % no need for a list dlg if only one format available
    s = 1;
    v = 1;
else
    [s,v] = listdlg('PromptString',{'Specify format of file', filename},...
        'SelectionMode','single',...
        'ListSize',[160 100],...
        'ListString',useList);
end
                
if(~v)
    warndlg('Tracking will be disabled until a tracking format has been specified');
    trackType=[];
    return;
end
    
trackType = useList{s};