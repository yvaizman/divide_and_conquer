function [eval] = classification_evaluation(features,labels,label2int,int2label)

addpath(genpath('C:\Users\Yonatan\Documents\ucsd\tools\liblinear-1.94'));

N                       = length(labels);
[d,N2]                  = size(features);
if N2~=N
    error('Second dimension of features must be equal to the dimension of labels');
end

n_labels                = length(int2label);
label_ints              = 1:n_labels;

% Set the labels:
y                       = zeros(N,1);
for ii = 1:N
    y(ii)               = label2int.(labels{ii});
end
% Set the features:
X                       = sparse(features');

n_trials                = 5;
test_portion            = 0.2;
ll_flags                = '-s 2 -B 1 -q';

eval                    = struct();
eval.accuracy_per_trial = zeros(1,n_trials);
eval.avr_ac_per_trial   = zeros(1,n_trials);
eval.conf_mat_per_trial = zeros(n_labels,n_labels,n_trials);
eval.details_per_trial  = cell(1,n_trials);
eval.model_per_trial    = cell(1,n_trials);

for ti = 1:n_trials
    fprintf('=== trial %d:\n',ti);
    [train_inds,test_inds]  = balanced_partition(y,test_portion);
    eval.details_per_trial{ti}.test_inds    = test_inds;
    X_tr                    = X(train_inds,:);
    y_tr                    = y(train_inds);
    
    % Train a model:
    model                   = liblinear_train(y_tr,X_tr,ll_flags);
    eval.model_per_trial{ti}        = model;
    
    % Predict on the test set:
    X_ts                    = X(test_inds,:);
    y_ts                    = y(test_inds,:);
    [y_hat]                 = liblinear_predict(y_ts,X_ts,model);
    [accuracy,confusion_mat]    = ...
        get_classification_scores(y_ts,y_hat,label_ints);
    inst_ac_vals            = diag(confusion_mat);
    inst_ac_vals            = inst_ac_vals(~isnan(inst_ac_vals));
    avr_inst_ac             = mean(inst_ac_vals);
    
    fprintf('accuracy: %f. Average instrument accuracy: %f.\n',accuracy,avr_inst_ac);
    
    % Store the performance for this trial:
    eval.accuracy_per_trial(ti)     = accuracy;
    eval.conf_mat_per_trial(:,:,ti) = confusion_mat;
    eval.avr_ac_per_trial(ti)       = avr_inst_ac;
end

% Summarize performance:
eval.accuracy               = mean(eval.accuracy_per_trial);
eval.accuracy_std           = std(eval.accuracy_per_trial);
eval.confusion_mat          = mean(eval.conf_mat_per_trial,3);
eval.confusion_mat_std      = std(eval.conf_mat_per_trial,[],3);
eval.avr_inst_accuracy      = mean(eval.avr_ac_per_trial);
eval.avr_inst_accuracy_std  = std(eval.avr_ac_per_trial);

end
