function [accuracy,mean_diff] = get_regression_scores(gt,pred,diff_thresh)

N           = length(gt);
N2          = length(pred);
if N ~= N2
    error('gt and pred must both be vectors of the same dimension');
end

diffs       = pred - gt;
mean_diff   = mean(diffs);
correct     = (abs(diffs)<diff_thresh);
accuracy    = mean(correct);

end