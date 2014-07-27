function [] = h__classification()

addpath(genpath('C:\Users\Yonatan\Documents\ucsd\tools\liblinear-1.94'));

data_supdir     = 'C:\Users\Yonatan\Documents\ucsd\music_data\UIowa\';
proc_supdir     = 'C:\Users\Yonatan\Documents\ucsd\music_data\UIowa_processing';
segments_file   = [data_supdir,filesep,'segmentations.mat'];
load(segments_file);

listfile        = [data_supdir,filesep,'uiowa.csv'];

[metadata,colinds]      = read_list_file(listfile);
fieldname               = 'type';
all_labels              = metadata(:,colinds.(fieldname));
for li = 1:length(all_labels)
    if (strcmp('2012',all_labels{li}(end-3:end)))
        all_labels{li}  = all_labels{li}(1:end-4);
    end
end

% Keep only part of the instances:
% is_piano                = strcmp('piano',all_labels);
% is_guitar               = strcmp('guitar',all_labels);
% keep                    = is_piano | is_guitar;
% metadata                = metadata(keep,:);
% all_labels              = all_labels(keep);

[label2int,int2label]   = create_label2int_map(all_labels);

M                               = 14;

mfcc_params.fs                  = 22050;
mfcc_params.fft_size            = 2048;
mfcc_params.n_mel_bins          = 128;
mfcc_params.num_ceps_coeffs     = 40;
mfcc_params.use_first_coeff     = false;
[mel_mat,dct_mat,c] = get_mel_scale_matrix_and_DCT_matrix(mfcc_params);
spec_size           = size(mel_mat,2);

mfcc_params.n_mel_bins          = 256;
[mel_mat256] = get_mel_scale_matrix_and_DCT_matrix(mfcc_params);
mfcc_params.n_mel_bins          = 34;
[mel_mat34] = get_mel_scale_matrix_and_DCT_matrix(mfcc_params);


temporal_file           = [proc_supdir,filesep,'temporal_gainfb.mat'];
[temp_features,temp_labels]       = get_features_and_labels(temporal_file,all_labels,segmentations);
pitch_file           = [proc_supdir,filesep,'pitch.mat'];
[pitch_features,pitch_labels]       = get_features_and_labels(pitch_file,all_labels,segmentations);
f0s                     = pitch_features(7,:);
% [peak_mfs,peak_spec] = pitch2peak_spectra(f0s,c.fft_freq,mel_mat);

lpc1_filters_file       = [proc_supdir,filesep,'representatives.As.lpcl1M' num2str(M) '.mat'];
lpc2_filters_file       = [proc_supdir,filesep,'representatives.As.lpcl2M' num2str(M) '.mat'];
lpc2np_filters_file     = [proc_supdir,filesep,'representatives.As.lpcl2M' num2str(M) '_nopreemph.mat'];
nolpc_specs_file        = [proc_supdir,filesep,'representatives.Ss.nolpc.mat'];

params.label2int        = label2int;
params.int2label        = int2label;
params.model_type       = 'linear';
params.is_diag          = false;

reloct_params               = params;
reloct_params.model_type    = 'linear_regression';
reloct_params.diff_thresh   = 1/24; % In octaves (equivalent of quarter tone).
% [lpc1_As,lpc1_labels,lpc1_Ss,lpc1_mfs] = ...
%     prepare_lpc_features(lpc1_filters_file,all_labels,spec_size,mel_mat,dct_mat);
% lpc1_mfs_eval           = classification_evaluation(lpc1_mfs,lpc1_labels,label2int,int2label);

[lpc2_As,lpc2_labels,lpc2_Ss,lpc2_mfs,lpc2_mfcc] = ...
    prepare_lpc_features(lpc2_filters_file,all_labels,spec_size,mel_mat,dct_mat);
%lpc2_mfs_eval           = classification_evaluation(lpc2_mfs,lpc2_labels,params);


[nolpc_Ss,nolpc_labels,nolpc_reloct,nolpc_file_inds,nolpc_seg_inds] = get_features_and_labels(nolpc_specs_file,all_labels,segmentations);
[nolpc_mfcc,nolpc_mfs]  = spectra2mfcc(nolpc_Ss,mel_mat,dct_mat);
%nolpc_mfs_eval          = classification_evaluation(nolpc_mfs,nolpc_labels,params);


[filter_mfs,excite_mfs,filter_power,excite_power,filter_log_spec,excite_log_spec] = ...
    separate_excitation_filter_with_dft_smoothing(nolpc_Ss,150,22050,mel_mat);
reloct_eval             = classification_evaluation(excite_power(1:500,:),440*2.^nolpc_reloct',reloct_params);

cb_inds = nolpc_mfs_eval.details_per_trial{1}.test_inds(1:3:end);
ex_inds = sort(setdiff(1:length(nolpc_labels),cb_inds))';
cb_lpc2_Ss = lpc2_Ss(:,cb_inds);
inv_cb = 1./cb_lpc2_Ss;

norms2 = ((abs(inv_cb').^2) * (abs(nolpc_Ss(:,ex_inds)).^2)).^0.5;
norms1 = (abs(inv_cb')) * (abs(nolpc_Ss(:,ex_inds)));
norms0_5 = ((abs(inv_cb').^0.5) * (abs(nolpc_Ss(:,ex_inds)).^0.5)).^2;

norm_eval = classification_evaluation([log(norms1);log(norms2);nolpc_mfcc(:,ex_inds);temp_features(:,ex_inds)],nolpc_labels(ex_inds),params);

end



function [peak_mfs,peak_spec] = pitch2peak_spectra(f0s,freq_bins,mel_mat)

N               = length(f0s);
spec_size       = length(freq_bins);
peak_spec       = ones(spec_size,N);

for ii = 1:N
    f0          = f0s(ii);
    harmonics   = f0:f0:freq_bins(end);
    for hi=1:length(harmonics)
        [val,bini]  = min(abs(freq_bins - harmonics(hi)));
        peak_spec(bini,ii)   = 10;
    end
end

peak_mel        = mel_mat*peak_spec;
peak_mfs        = 10*log10(peak_mel);

end

function [contrib] = estimate_feature_contribution(features,eval)

mean_val    = mean(abs(features'));
mean_coef   = mean(abs(eval.model_per_trial{1}.w));
contrib     = [mean_val,1].*mean_coef;

end

function [lpc_As,lpc_labels,lpc_Ss,lpc_mfs,lpc_mfcc] = prepare_lpc_features(lpc_filters_file,all_labels,spec_size,mel_mat,dct_mat)

load(lpc_filters_file);
lpc_As          = collection;
lpc_labels      = all_labels(file_inds);
%lpc_Ss          = filters2spectra(lpc_As,ones(1,size(lpc_As,2)),gains,spec_size);
lpc_Ss          = filters2spectra(lpc_As,ones(1,size(lpc_As,2)),ones(1,size(lpc_As,2)),spec_size);
[lpc_mfcc,lpc_mfs]  = spectra2mfcc(lpc_Ss,mel_mat,dct_mat);

end

function [features,labels,relocts,file_inds,seg_inds] = get_features_and_labels(features_file,all_labels,segmentations)

load(features_file);
features        = collection;
labels          = all_labels(file_inds);

relocts         = zeros(size(file_inds));
for ii = 1:length(file_inds)
    if isempty(segmentations{file_inds(ii)}.relocts)
        relocts(ii) = NaN;
    else
        relocts(ii) = segmentations{file_inds(ii)}.relocts(seg_inds(ii));
    end
end

end