function [mel_mat,c,dct_mat] = get_mel_scale_matrix_and_DCT_matrix(p)
%%
%% This is copied from Slaney's ma_mel just for the purpose of getting the
%% mel scale filters matrix and the DCT matrix
%%
%%
%% based on function "mfcc" in Auditory Toolbox by Malcolm Slaney
%% (http://www.slaney.org/malcom) see Auditory Toolbox for details
%%
%% INPUT
%%  p (struct) parameters:
%%  p.fs: sampling rate of time signal (Hz)
%%  p.fft_size: the size of DFT over which we will apply Mel scaling
%%  p.n_mel_bins: number of Mel-scaled frequency bins
%%  p.min_freq: minimal frequency to cover with Mel bins
%%  p.max_freq: maximal frequency to cover with Mel bins
%%  p.num_ceps_coeffs: number of cepstral coefficients (output of DCT)
%%  p.use_first_coeff: should we include the first cepstral coefficient?
%% defaults
if ~isfield(p,'fs'),    p.fs = 22050; end
if ~isfield(p,'fft_size'),  p.fft_size = 2048; end
if ~isfield(p,'n_mel_bins'),    p.n_mel_bins = 128; end
if ~isfield(p,'min_freq'),  p.min_freq = 0; end
if ~isfield(p,'max_freq'),  p.max_freq = p.fs/2; end
if ~isfield(p,'num_ceps_coeffs'), p.num_ceps_coeffs = 40; end
if ~isfield(p,'use_first_coeff'), p.use_first_coeff = 1; end

% The DFT frequencies:
c.fft_freq = linspace(0,p.fs/2,p.fft_size/2+1);

c.num_filt  = p.n_mel_bins;

% Calculate the Mel scale:
tmp.f       = p.min_freq:p.max_freq;
tmp.mel     = log(1+tmp.f/700)*1127.01048;
tmp.m_idx = linspace(1,max(tmp.mel),c.num_filt+2);

tmp.f_idx = (1:c.num_filt+2)*0;
for i=1:c.num_filt+2,
    [tmp.dummy tmp.f_idx(i)] = min(abs(tmp.mel - tmp.m_idx(i)));
end
c.freqs = tmp.f(tmp.f_idx);

% Handle cases in which the freq resolution is so dense, there are two bins
% with the same freq:
sameAsNext = find(c.freqs(1:end-1) == c.freqs(2:end));
c.freqs(sameAsNext) = [];
c.num_filt = c.num_filt - length(sameAsNext);

% Calculate the Mel filters:
c.lower  = c.freqs(1:c.num_filt);
c.center = c.freqs(2:c.num_filt+1);
c.upper  = c.freqs(3:c.num_filt+2);

%% ignore filters outside of spectrum
[tmp.dummy, tmp.idx] = max(find(c.center <= p.fs/2));
c.num_filt = min(tmp.idx,c.num_filt);

mel_mat = zeros(c.num_filt,p.fft_size/2+1);
c.triangleHeight = 2./(c.upper-c.lower);

for i=1:c.num_filt,
    mel_mat(i,:) = ...
        (c.fft_freq > c.lower(i) & c.fft_freq <= c.center(i)).* ...
        c.triangleHeight(i).*(c.fft_freq-c.lower(i))/(c.center(i)-c.lower(i)) + ...
        (c.fft_freq > c.center(i) & c.fft_freq < c.upper(i)).* ...
        c.triangleHeight(i).*(c.upper(i)-c.fft_freq)/(c.upper(i)-c.center(i));
end




dct_mat = 1/sqrt(c.num_filt/2)*cos((1-p.use_first_coeff:(p.num_ceps_coeffs-1))' * (2*(0:(c.num_filt-1))+1) * pi/2/c.num_filt);
if p.use_first_coeff,
    dct_mat(1,:) = dct_mat(1,:) * sqrt(2)/2;
end

%% We only need up to here
return;