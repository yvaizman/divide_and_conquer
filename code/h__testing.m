function [] = h__testing()

data_supdir     = 'C:\Users\Yonatan\Documents\ucsd\music_data\UIowa\wav';

%wav_file        = [data_supdir,filesep,'piano\piano\Piano.ff.Eb4.wav'];
wav_file        = [data_supdir,filesep,'strings\violin\Violin.arco.ff.sulE.C6Gb7.wav'];

params.sr       = 22050;
params.M        = 20;
params.hoplen   = 1024;
params.winlen   = 2*params.hoplen;
params.do_preemph = true;

[w,sr_orig]     = wavread(wav_file);
w               = mean(w,2);
w               = resample(w,params.sr,sr_orig);

seg             = w(50001:55000);

params.criterion = 'L2';
[A_l2,g_l2,e_l2,out_params_l2,S_l2] = lpc_and_spec(w,params);

params.criterion = 'acL1';
params.lags     = 0:70;
[A_acl1,g_acl1,e_acl1,out_params_acl1,S_acl1] = lpc_and_spec(w,params);

end

function [A,g,e,out_params,S] = lpc_and_spec(w,params)

[A,g,e,out_params] = lpc_analysis_by_frames(w,params);
S               = filters2spectra(A,ones(1,size(A,2)),ones(1,size(A,2)),512);

end