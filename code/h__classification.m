function [] = h__classification()

addpath(genpath('C:\Users\Yonatan\Documents\ucsd\tools\liblinear-1.94'));

data_supdir     = 'C:\Users\Yonatan\Documents\ucsd\music_data\UIowa\';
proc_supdir     = 'C:\Users\Yonatan\Documents\ucsd\music_data\UIowa_processing';

listfile        = [data_supdir,filesep,'uiowa.csv'];

[metadata,colinds]      = read_list_file(listfile);
fieldname               = 'type';
all_labels              = metadata(:,colinds.(fieldname));
for li = 1:length(all_labels)
    if (strcmp('2012',all_labels{li}(end-3:end)))
        all_labels{li}  = all_labels{li}(1:end-4);
    end
end
[label2int,int2label]   = create_label2int_map(all_labels);

M                               = 14;

mfcc_params.fs                  = 22050;
mfcc_params.fft_size            = 2048;
mfcc_params.n_mel_bins          = 128;
mfcc_params.num_ceps_coeffs     = 14;
mfcc_params.use_first_coeff     = false;
[mel_mat,dct_mat,c] = get_mel_scale_matrix_and_DCT_matrix(mfcc_params);
spec_size           = size(mel_mat,2);

lpc1_filters_file       = [proc_supdir,filesep,'representatives.As.lpcl1M' num2str(M) '.mat'];
lpc2_filters_file       = [proc_supdir,filesep,'representatives.As.lpcl2M' num2str(M) '.mat'];
lpc2np_filters_file     = [proc_supdir,filesep,'representatives.As.lpcl2M' num2str(M) '_nopreemph.mat'];
nolpc_specs_file        = [proc_supdir,filesep,'representatives.Ss.nolpc.mat'];

% [lpc1_As,lpc1_labels,lpc1_Ss,lpc1_mfs] = ...
%     prepare_lpc_features(lpc1_filters_file,all_labels,spec_size,mel_mat,dct_mat);
% lpc1_mfs_eval           = classification_evaluation(lpc1_mfs,lpc1_labels,label2int,int2label);

[lpc2_As,lpc2_labels,lpc2_Ss,lpc2_mfs] = ...
    prepare_lpc_features(lpc2_filters_file,all_labels,spec_size,mel_mat,dct_mat);
lpc2_mfs_eval           = classification_evaluation(lpc2_mfs,lpc2_labels,label2int,int2label);

[lpc2np_As,lpc2np_labels,lpc2np_Ss,lpc2np_mfs] = ...
    prepare_lpc_features(lpc2np_filters_file,all_labels,spec_size,mel_mat,dct_mat);
lpc2np_mfs_eval           = classification_evaluation(lpc2np_mfs,lpc2np_labels,label2int,int2label);

[nolpc_Ss,nolpc_labels] = get_features_and_labels(nolpc_specs_file,all_labels);
[nolpc_mfcc,nolpc_mfs]  = spectra2mfcc(nolpc_Ss,mel_mat,dct_mat);
nolpc_mfs_eval          = classification_evaluation(nolpc_mfs,nolpc_labels,label2int,int2label);

end

function [lpc_As,lpc_labels,lpc_Ss,lpc_mfs] = prepare_lpc_features(lpc_filters_file,all_labels,spec_size,mel_mat,dct_mat)

load(lpc_filters_file);
lpc_As          = collection;
lpc_labels      = all_labels(file_inds);
%lpc_Ss          = filters2spectra(lpc_As,ones(1,size(lpc_As,2)),gains,spec_size);
lpc_Ss          = filters2spectra(lpc_As,ones(1,size(lpc_As,2)),ones(1,size(lpc_As,2)),spec_size);
[lpc_mfcc,lpc_mfs]  = spectra2mfcc(lpc_Ss,mel_mat,dct_mat);

end

function [features,labels] = get_features_and_labels(features_file,all_labels)

load(features_file);
features        = collection;
labels          = all_labels(file_inds);

end