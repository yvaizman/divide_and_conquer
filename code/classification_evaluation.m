function [eval] = classification_evaluation(features,labels,params)

addpath(genpath('C:\Users\Yonatan\Documents\ucsd\tools\liblinear-1.94'));

N                       = length(labels);
[d,N2]                  = size(features);
if N2~=N
    error('Second dimension of features must be equal to the dimension of labels');
end

n_labels                = length(params.int2label);
params.label_ints       = 1:n_labels;

% Set the labels:
switch params.model_type
    case 'linear'
        y               = zeros(N,1);
        for ii = 1:N
            y(ii)       = params.label2int.(labels{ii});
        end
    case 'linear_regression'
        y               = labels;
end
% Set the features:
X                       = features';

n_trials                = 5;
test_portion            = 0.2;
eval                    = struct();
eval.accuracy_per_trial = zeros(1,n_trials);
eval.details_per_trial  = cell(1,n_trials);
eval.model_per_trial    = cell(1,n_trials);

switch params.model_type
    case 'linear'
        params.ll_flags = '-s 2 -B 1 -q';
        X               = sparse(X);
        eval.conf_mat_per_trial = zeros(n_labels,n_labels,n_trials);
        eval.avr_ac_per_trial   = zeros(1,n_trials);
    case 'linear_regression'
        params.ll_flags = '-s 11 -B 1 -q';
        X               = sparse(X);
        eval.avr_diff_per_trial = zeros(1,n_trials);
    case 'gauss'
%        params.is_diag  = true;
end


for ti = 1:n_trials
    fprintf('=== trial %d:\n',ti);
    [train_inds,test_inds]  = balanced_partition(y,test_portion);
    eval.details_per_trial{ti}.test_inds    = test_inds;
    X_tr                    = X(train_inds,:);
    y_tr                    = y(train_inds);
    
    % Train a model:
    model                   = train_model(y_tr,X_tr,params);
    eval.model_per_trial{ti}        = model;
    
    % Predict on the test set:
    X_ts                    = X(test_inds,:);
    y_ts                    = y(test_inds,:);
    [y_hat]                 = predict_with_model(X_ts,params.model_type,model);

    switch params.model_type
        case 'linear'
            [accuracy,confusion_mat]    = ...
                get_classification_scores(y_ts,y_hat,params.label_ints);
            inst_ac_vals            = diag(confusion_mat);
            inst_ac_vals            = inst_ac_vals(~isnan(inst_ac_vals));
            avr_inst_ac             = mean(inst_ac_vals);

            fprintf('accuracy: %f. Average instrument accuracy: %f.\n',accuracy,avr_inst_ac);

            % Store the performance for this trial:
            eval.accuracy_per_trial(ti)     = accuracy;
            eval.conf_mat_per_trial(:,:,ti) = confusion_mat;
            eval.avr_ac_per_trial(ti)       = avr_inst_ac;
        case 'linear_regression'
            [accuracy,mean_diff]            = ...
                get_regression_scores(y_ts,y_hat,diff_thresh);
            fprintf('accuracy: %f. Average value difference: %f.\n',accuracy,mean_diff);

            % Store the performance for this trial:
            eval.accuracy_per_trial(ti)     = accuracy;
            eval.avr_diff_per_trial(ti)     = mean_diff;
    end
end

% Summarize performance:
eval.accuracy               = mean(eval.accuracy_per_trial);
eval.accuracy_std           = std(eval.accuracy_per_trial);
switch params.model_type
    case 'linear'
        eval.confusion_mat          = mean(eval.conf_mat_per_trial,3);
        eval.confusion_mat_std      = std(eval.conf_mat_per_trial,[],3);
        eval.avr_label_accuracy     = mean(eval.avr_ac_per_trial);
        eval.avr_label_accuracy_std = std(eval.avr_ac_per_trial);
    case 'linear_regression'
        eval.avr_diff               = mean(eval.avr_diff_per_trial);
        eval.avr_diff_std           = std(eval.avr_diff_per_trial);
end

end

function [model] = train_model(y_tr,X_tr,params)

switch params.model_type
    case 'linear'
        model   = liblinear_train(y_tr,X_tr,params.ll_flags);
    case 'linear_regression'
        model   = liblinear_train(y_tr,X_tr,params.ll_flags);
    case 'gauss'
        model   = train_multiclass_gauss_model(y_tr,X_tr,params);
end

end

function [y_hat] = predict_with_model(X_ts,model_type,model)

switch model_type
    case 'linear'
        y_hat   = liblinear_predict(zeros(size(X_ts,1),1),X_ts,model);
    case 'linear_regression'
        y_hat   = liblinear_predict(zeros(size(X_ts,1),1),X_ts,model);
    case 'gauss'
        y_hat   = predict_with_multiclass_gauss_model(X_ts,model);
end

end