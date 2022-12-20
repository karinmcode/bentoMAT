function     [newFeatures,Nmax]=myCheckFeaturesDataValidity(newFeatures)
%% temp=myCheckFeaturesDataValidity(temp)

% make sure all data is vector of N sample
FI =  fieldnames(newFeatures);
nfi = numel(FI);
N = cell(nfi,1);
Nmax = nan(nfi,1);
for ifi=1:nfi
    fi = FI{ifi};
    N{ifi}=size(newFeatures.(fi));
    Nmax(ifi) = max(N{ifi});
end

if all(Nmax==Nmax(1))
    for ifi=1:nfi
        fi = FI{ifi};
        newFeatures.(fi) = newFeatures.(fi)(:)';

        Nmax(ifi) = numel(newFeatures.(fi));
    end
end

        