function [] = h__testing()

yin_path        ='C:\Users\yonatan\Documents\ucsd\tools\yin_changed\';
data_supdir     = 'C:\Users\Yonatan\Documents\ucsd\music_data\UIowa\wav';

%wav_file        = [data_supdir,filesep,'piano\piano\Piano.ff.Eb4.wav'];
%wav_file        = [data_supdir,filesep,'strings\violin\Violin.arco.ff.sulE.C6Gb7.wav'];
%wav_file        = [data_supdir,filesep,'woodwind\EbAltoSaxophone\AltoSax.NoVib.ff.C5Ab5.wav'];
%wav_file        = [data_supdir,filesep,'woodwind\BbSopranoSaxophone\SopSax.NoVib.ff.C5B5.wav'];
wav_file        = [data_supdir,filesep,'woodwind\BbClarinet\BbClar.ff.C4B4.wav'];
wav_file2       = [data_supdir,filesep,'strings\cello\Cello.arco.pp.sulA.C4Bb4.wav'];

params.sr       = 22050;
params.M        = 20;
params.hoplen   = 1024;
params.winlen   = 2*params.hoplen;
params.do_preemph = true;
params.yin_path = yin_path;

params.lags     = 0:40;
params.version  = 'nonsymmetric';

[w,sr_orig]     = wavread(wav_file);
w               = mean(w,2);
w               = resample(w,params.sr,sr_orig);

[w2,sr_orig]    = wavread(wav_file2);
w2              = mean(w2,2);
w2              = resample(w2,params.sr,sr_orig);

seg             = w(50001:55000);
wseg            = hamming(length(seg)).*seg;

seg2            = w(343500:446800);
gain_sr         = 200; %Hz.
gain_filt_dur   = 1.5; % seconds
gain_filt_len   = gain_filt_dur*gain_sr;
gain_filterbank = get_gamma_filterbank(gain_filt_len,gain_sr);
features        = calc_temporal_features(seg2,params.sr,gain_sr,gain_filterbank);

test_dwt_lpc(wseg,params);

S1              = spectrogram(w,params.winlen,params.hoplen,params.winlen);
S2              = spectrogram(w2,params.winlen,params.hoplen,params.winlen);

params.criterion    = 'L1';
params.use_window   = false;
[res_l1w1,S_l1w1,Se_l1w1] = lpc_and_spec(w,params);
[res_l1w2,S_l1w2,Se_l1w2] = lpc_and_spec(w2,params);

mfcc_params.fs  = params.sr;
mfcc_params.fft_size    = params.winlen;
mfcc_params.n_mel_bins  = 128;
mfcc_params.num_ceps_coeffs     = 13;
mfcc_params.use_first_coeff     = false;
[mel_mat,dct_mat,c] = get_mel_scale_matrix_and_DCT_matrix(mfcc_params);
p.mel_mat = mel_mat;
p.dct_mat = dct_mat;
[mfcc1,mfs1] = spectra2mfcc(S1,p);
[mfcc2,mfs2] = spectra2mfcc(S2,p);

[mfcc1_l1,mfs1_l1] = spectra2mfcc(S_l1w1,p);
[mfcc2_l1,mfs2_l1] = spectra2mfcc(S_l1w2,p);



%[S,powers,mean_psd,mean_ac,mac_al2] = estimate_using_bulk_examples(w,params);

%yinr_xn = estimate_pitch_track(w,params);
%S = spectrogram(w,params.winlen,params.hoplen,params.winlen);

params.criterion    = 'L2';
params.use_window   = true;
%[res_l2w,S_l2w,Se_l2w] = lpc_and_spec(w,params);

params.criterion    = 'L1';
params.use_window   = true;
[res_l1w,S_l1w,Se_l1w] = lpc_and_spec(w,params);

params.criterion    = 'L1';
params.use_window   = false;
[res_l1,S_l1,Se_l1] = lpc_and_spec(w,params);

params.criterion    = 'acL1';
params.version      = 'nonsymmetric';
[A_acl1,g_acl1,e_acl1,out_params_acl1,S_acl1] = lpc_and_spec(w,params);

params.version      = 'convex';
[A_acl1c,g_acl1c,e_acl1c,out_params_acl1c,S_acl1c] = lpc_and_spec(w,params);


end

function [] = test_dwt_lpc(wseg,params)

params.independent_scales = false;
params.poles_per_level = [10,16];
[a,poles_per_level,filt_per_level,origfilt_per_level,gain_per_level] = dwt_lpc(wseg,params);
y = fft(wseg);spec_size=length(wseg)/2+1;ha=freqz(1,a,spec_size);
m = length(a)-1;
alpc = lpc(wseg,m);
hlpc = freqz(1,alpc,spec_size);

N = length(poles_per_level);
hs = zeros(spec_size,N);
for ii = 1:N
    hs(:,ii) = freqz(gain_per_level(ii),filt_per_level{ii},spec_size);
end

hs_flat = reshape(hs(:,end:-1:1),N*spec_size,1);

sumh = sum(hs,2);
sumabsh = sum(abs(hs),2);
sumhdb = sum(20*log10(abs(hs)),2);

figure;plot(20*log10(abs(y(1:spec_size))),'b');hold on;plot(20*log10(abs(ha)),'m');
hold on;plot(20*log10(abs(hlpc)),'g');title('resulted filter');

figure;plot(20*log10(abs(y(1:spec_size))),'b');hold on;plot(20*log10(abs(hs)));title('hs');
% figure;plot(20*log10(abs(y(1:spec_size))),'b');hold on;plot(20*log10(abs(sumh)),'g');title('sum hs');
% figure;plot(20*log10(abs(y(1:spec_size))),'b');hold on;plot(20*log10(sumabsh),'g');title('sum abs hs');
% figure;plot(20*log10(abs(y(1:spec_size))),'b');hold on;plot(sumhdb,'g');hold on;plot(20*log10(abs(ha)),'m');title('sum abs-hs (dB) and resulted filter in magenta');

en = filter(a,1,wseg);
elpc = filter(alpc,1,wseg);
yen = fft(en);
yelpc = fft(elpc);
figure;plot(20*log10(abs(yen(1:spec_size))));title('resulted error');
figure;plot(20*log10(abs(yelpc(1:spec_size))));title('lpc error');

[Lo_D,Hi_D,Lo_R,Hi_R] = wfilters('db45');
hlor = freqz(Lo_R,1,spec_size);
hhir = freqz(Hi_R,1,spec_size);

% Original filters:
ho = zeros(spec_size,N);
for ii = 1:N
    ho(:,ii) = freqz(gain_per_level(ii),origfilt_per_level{ii},spec_size);
end

ho_flat = reshape(ho(:,end:-1:1),N*spec_size,1);
end

function [res,S,Se] = lpc_and_spec(w,params)

[res]               = lpc_analysis_by_frames(w,params);
S                   = filters2spectra(res.A,ones(1,size(res.A,2)),ones(1,size(res.A,2)),1025);
Se                  = spectrogram(res.e,params.winlen,params.hoplen,params.winlen);
%yinr                = estimate_pitch_track(res.e,params);

end

function [r] = estimate_pitch_track(w,params)

%addpath(genpath('C:\Users\yonatan\Documents\ucsd\tools\yin_changed\'));
yin_params.sr       = params.sr;
yin_params.hop      = 64;
yin_params.wsize    = 0.75*params.winlen;

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

