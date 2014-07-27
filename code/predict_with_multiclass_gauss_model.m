function [y_hat] = predict_with_multiclass_gauss_model(X,model)

N               = size(X,1);
n_labels        = length(model.label_ints);
label_sq_dists     = Inf(N,n_labels);
not_label_sq_dists = Inf(N,n_labels);

% For every label, calculate the distance of all the points from the
% label's center and from the not-in-label's center:
for li = 1:n_labels
    if isnan(model.label_gauss.means(li,1))
        continue;
    end
    label_sq_dists(:,li)       = calc_distances_from_center(X,model.label_gauss,li,model.is_diag);
    not_label_sq_dists(:,li)   = calc_distances_from_center(X,model.not_label_gauss,li,model.is_diag);
end

% Calculate for every point the relative distances (the difference between
% the distance to the label's center and the distance to the other labels'
% center):
relative_sq_dists  = label_sq_dists - not_label_sq_dists;

% For every point, select the label for which the relative distance is
% smallest:
[vals,inds]     = min(relative_sq_dists,[],2);
y_hat           = reshape(model.label_ints(inds),N,1);

end

function [sq_dists] = calc_distances_from_center(X,gausses_model,li,is_diag)

mean_vec        = gausses_model.means(li,:);

onevec          = ones(size(X,1),1);
differences     = X - (onevec*mean_vec);
if is_diag
    inv_var     = gausses_model.inv_vars(li,:);
    sq_dists    = sum((differences.^2) .* (onevec*inv_var),2);
else
    inv_cov     = gausses_model.inv_covars{li};
    sq_dists    = diag(differences * inv_cov * differences');
end

end