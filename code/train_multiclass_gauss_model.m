function [model] = train_multiclass_gauss_model(y,X,params)

model               = struct();
model.is_diag       = params.is_diag;
model.label_ints    = params.label_ints;
model.label_gauss   = struct();
model.not_label_gauss   = struct();

dim                 = size(X,2);
n_labels            = length(params.label_ints);
model.label_gauss.means     = NaN(n_labels,dim);
model.not_label_gauss.means = NaN(n_labels,dim);
if params.is_diag
    model.label_gauss.inv_vars      = zeros(n_labels,dim);
    model.not_label_gauss.inv_vars  = zeros(n_labels,dim);
else
    model.label_gauss.inv_covars    = cell(n_labels,1);
    model.not_label.gauss.inv_covars= cell(n_labels,1);
end

for li = 1:length(params.label_ints)
    label           = params.label_ints(li);
    pos             = (y == label);
    if (sum(pos) <= 0)
        continue;
    end
    X_i             = X(pos,:);
    X_not_i         = X(~pos,:);
    [i_mean,i_scaling]          = calc_gaussian_model(X_i,params.is_diag);
    [not_i_mean,not_i_scaling]  = calc_gaussian_model(X_not_i,params.is_diag);
    
    model.label_gauss.means(li,:)       = i_mean;
    model.not_label_gauss.means(li,:)   = not_i_mean;
    if params.is_diag
        model.label_gauss.inv_vars(li,:)    = i_scaling;
        model.not_label_gauss.inv_vars(li,:)= not_i_scaling;
    else
        model.label_gauss.inv_covars{li}    = i_scaling;
        model.not_label_gauss.inv_covars{li}= not_i_scaling;
    end
end

end

function [g_mean,g_scaling] = calc_gaussian_model(X,is_diag)

epsilon                 = 10e-6;
g_mean                  = mean(X,1);
if is_diag
    var_vec             = var(X,1) + epsilon;
    g_scaling           = ones(size(g_mean)); %%1./var_vec;
else
    cov_mat             = cov(X);
    cov_mat             = cov_mat + epsilon*eye(size(X,2));
    g_scaling           = inv(cov_mat);
end

end