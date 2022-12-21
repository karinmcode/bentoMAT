function rast = convertToRast(bhvr,max_index)
%% rast = convertToRast(bhvr,max_index)
% (C) Ann Kennedy, 2019
% California Institute of Technology
% Licensing: https://github.com/annkennedy/bento/blob/master/LICENSE.txt



indexes = 1:max_index;
rast = false(size(indexes));
nrows = size(bhvr,1);
for i = 1:nrows
    inds = (indexes>=bhvr(i,1))&(indexes<=bhvr(i,2))&(indexes<=max_index);
    rast(inds)=true;
end