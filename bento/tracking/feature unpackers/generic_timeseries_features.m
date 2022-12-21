function [features, names] = generic_timeseries_features(args)
% features : rows are number of features and columns time/samples
f = fieldnames(args);
if length(f)>1
% 	disp('not sure which variables to unpack, features will not be loaded.');
%     questdlg()
% 	features = [];
% 	names = '';
% 	return
end
num_feats = numel(f);
num_sample_xfeat = nan(num_feats,1);
for ifi = 1:num_feats
    fi = f{ifi};
    num_sample_xfeat(ifi) = numel(args.(fi));
    sz = size(args.(fi));
end
nsamp = num_sample_xfeat(1);
if ~all(num_sample_xfeat==nsamp)
	disp('not sure which variables to unpack, features will not be loaded.');
    disp(num_sample_xfeat)

    ANS=questdlg('Fix data to make them all same length?');
    if ANS(1)=='Y'
        n = min(num_sample_xfeat);
        for ifi = 1:num_feats
            fi = f{ifi};
            args.(fi) = args.(fi)(1:n);
        end
    else
	features = [];
	names = '';
	return;
    end
end

features = nan(num_feats,nsamp);
for ifi = 1:num_feats
    fi = f{ifi};
    features(ifi,:) = args.(fi);
end

names = f;