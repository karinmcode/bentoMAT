function myUnfoldAnnotStruct(gui)
if ~isfield(gui,'data')
    gui = guidata(gui.h0);
end
gui.annot.activeCh;
pch = gui.annot.channels;
nch = numel(pch);
for ich = 1:nch
    thisch = pch{ich};
    fprintf('\n-----------')
    if strcmp(thisch,gui.annot.activeCh)
        fprintf('CHANNEL %s (active)',thisch);
    else
        fprintf('CHANNEL %s ',thisch);
    end

    if isfield(gui.data.annot,thisch)
        ann = gui.data.annot.(thisch);
        labels = fieldnames(ann);
        nla = numel(labels);
        fprintf('\n')

        for ila = 1:nla
            lab = labels{ila};

            fprintf(' %s : %g ,',lab,size(ann.(lab),1))
        end
    else
        fprintf('     >CHANNEL %s  IS NOT FIELD',thisch);
    end
end
end