function [train_inds,test_inds] = balanced_partition(y,test_portion)

lints           = unique(y);
nvals           = length(lints);

train_inds      = [];
test_inds       = [];
% Go over the distinct label values.
% For each label partition the relevant examples to train set and test set
% according to the desired portion:
for li = 1:nvals
    lint        = lints(li);
    inds        = find(y == lint);
    ni          = length(inds);
    inds        = inds(randperm(ni));
    nts         = ceil(ni*test_portion);
    
    test_inds   = [test_inds;inds(1:nts)];
    train_inds  = [train_inds;inds((nts+1):end)];
end

end