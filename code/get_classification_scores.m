function [accuracy,confusion_mat,gt_counts] = get_classification_scores(gt,pred,label_ints)

N           = length(gt);
N2          = length(pred);
if N ~= N2
    error('gt and pred must both be vectors of the same dimension');
end

accuracy    = mean(gt == pred);
n_labels    = length(label_ints);
gt_counts   = zeros(1,n_labels);

confusion_mat   = zeros(n_labels,n_labels);
for li = 1:n_labels
    label_int1      = label_ints(li);
    
    gt_mark         = gt == label_int1;
    n_gt            = sum(gt_mark);
    gt_counts(li)   = n_gt;
    for lj = 1:n_labels
        label_int2  = label_ints(lj);
        pred_mark   = pred == label_int2;
        hits        = gt_mark & pred_mark;
        n_hits      = sum(hits);
        confusion_mat(li,lj)    = n_hits/n_gt;
    end
end

end