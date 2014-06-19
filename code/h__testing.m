function [] = h__testing()

data_supdir     = 'C:\Users\Yonatan\Documents\ucsd\music_data\UIowa\wav';

%wav_file        = [data_supdir,filesep,'piano\piano\Piano.ff.Eb4.wav'];
%wav_file        = [data_supdir,filesep,'strings\violin\Violin.arco.ff.sulE.C6Gb7.wav'];
%wav_file        = [data_supdir,filesep,'woodwind\EbAltoSaxophone\AltoSax.NoVib.ff.C5Ab5.wav'];
wav_file        = [data_supdir,filesep,'woodwind\BbSopranoSaxophone\SopSax.NoVib.ff.C5B5.wav'];

params.sr       = 22050;
params.M        = 20;
params.hoplen   = 1024;
params.winlen   = 2*params.hoplen;
params.do_preemph = true;
params.lags     = 0:40;
params.version  = 'nonsymmetric';

[w,sr_orig]     = wavread(wav_file);
w               = mean(w,2);
w               = resample(w,params.sr,sr_orig);

%[S,powers,mean_psd,mean_ac,mac_al2] = estimate_using_bulk_examples(w,params);
seg             = w(50001:55000);
wseg            = hamming(length(seg)).*seg;

yinr_xn = estimate_pitch_track(w,params);
S = spectrogram(w,params.winlen,params.hoplen,params.winlen);

params.criterion    = 'L2';
params.use_window   = true;
[res_l2w,S_l2w,Se_l2w] = lpc_and_spec(w,params);

params.criterion    = 'L1';
params.use_window   = true;
[A_l1w,g_l1w,e_l1w,out_params_l1w,S_l1w,Se_l1w] = lpc_and_spec(w,params);

params.criterion    = 'L1';
params.use_window   = false;
[A_l1,g_l1,e_l1,out_params_l1,S_l1,Se_l1,yinr_l1] = lpc_and_spec(w,params);

params.criterion    = 'acL1';
params.version      = 'nonsymmetric';
[A_acl1,g_acl1,e_acl1,out_params_acl1,S_acl1] = lpc_and_spec(w,params);

params.version      = 'convex';
[A_acl1c,g_acl1c,e_acl1c,out_params_acl1c,S_acl1c] = lpc_and_spec(w,params);


end

function [res,S,Se] = lpc_and_spec(w,params)

[res]               = lpc_analysis_by_frames(w,params);
S                   = filters2spectra(res.A,ones(1,size(res.A,2)),ones(1,size(res.A,2)),512);
Se                  = spectrogram(res.e,params.winlen,params.hoplen,params.winlen);
%yinr                = estimate_pitch_track(res.e,params);

end

function [r] = estimate_pitch_track(w,params)

addpath(genpath('C:\Users\yonatan\Documents\ucsd\tools\yin_changed\'));
yin_params.sr       = params.sr;
%yin_params.hop      = 1;
yin_params.wsize    = params.winlen - 10;

r = yin(w,yin_params);

end

function [S,powers,mean_psd,mean_ac,mac_al2] = estimate_using_bulk_examples(w,params)

S       = spectrogram(w,params.winlen,params.hoplen,params.winlen);
PSD     = abs(S).^2;
powers  = sum(PSD);

max_pow = max(powers);
thresh  = 10e-5*max_pow;
inds    = find(powers > thresh);

S       = S(:,inds);
powers  = powers(inds);
gains   = sqrt(powers);

S       = S ./ (ones(size(S,1),1)*gains);

mean_psd    = mean(abs(S).^2,2);
sym_psd     = [mean_psd;mean_psd(end:-1:2)];
mean_ac     = ifft(sym_psd);

mac_al2     = levinson(mean_ac,params.M);


end

