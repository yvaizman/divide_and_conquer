% function [A,g,e,preemph] = lpc_analysis_by_frames(wav,sr,n_poles,do_preemph,hoplen,winlen)
%
% Perform linear prediction analysis of a wave signal, by independent
% frames. For each frame an optimal all-pole filter is fitted to describe
% the generation of the frame's wave as an excitation signal passed through
% the all-pole filter.
%
% Input:
% -----
% wav: Either a string (name of a wav file to load) or a wave vector (L x 1).
% sr: Sampling rate (in Hz). If a wave vector was given, this will be
%     regarded as the sampling rate of the vector. If a filename was given,
%     the loaded wave will be resampled to this desired sampling rate.
% M: Positive integer. Order of the filters (number of poles).
% do_preemph: Boolean. Should we do a pre-emphasis? 
%             If true: a high-pass FIR 1-zero filter is applied to the given signal
%             before the frame-by-frame analysis.
% hoplen: Positive integer. Number of samples to hop from frame to frame.
% winlen: Positive integer. Size of each frame.
%
% Output:
% ------
% A: ((M+1) x n_frames) matrix. The fitted all-pole filters. Each column
%    has the coefficients of a frame's fitted filter (direct form).
% g: (1 x n_frames) vector. The gain of every frame.
% e: (L x 1) vector. The prediction-error/excitation/residual signal.
% preemph: vector. The coefficients of the pre-emphasis FIR filter that was
%          used (or scalar of 1 if no pre-emphasis was done).
%
% ------------------------------------------------------------------------
% Written by Yonatan Vaizman, 2014.
function [A,g,e,preemph] = lpc_analysis_by_frames(wav,sr,M,do_preemph,hoplen,winlen)

preemph     = [1 -0.9];

if ischar(wav)
    wav_file    = wav;
    [w,sr_orig] = wavread(wav_file);
    w           = mean(w,2);
    w           = resample(w,sr,sr_orig);
else
    w           = wav;
end

if do_preemph
    % Pre-emphasize:
    w           = filter(preemph,1,w);
else
    preemph     = 1;
end

L           = length(w);
n_frames    = floor(L/hoplen);

%% Prepare the silence threshold for this wave:
squares     = w.^2;
seg_dur     = 0.005; % In seconds;
seg_len     = round(seg_dur * sr);
smoother    = (1/seg_len)*ones(1,seg_len);
powers      = filter(smoother,1,squares);
max_power   = max(powers);
thresh      = 10^(-6)*max_power; % Silence threshold will be 60dB below max

%%
A           = [ones(1,n_frames);zeros(M,n_frames)];
g           = zeros(1,n_frames);
e           = zeros(L,1);

window      = hamming(winlen);
for fi = 1:n_frames
    from    = (fi-1)*hoplen + 1;
    to      = from+winlen-1;
    if to > L
        to  = L;
        window  = hamming(to-from+1);
    end
    frame   = w(from:to);
    wframe  = frame .* window;
    % Simple check if there is any power:
    epsilon = thresh;
    if var(wframe) <= epsilon
        % Leave the default filter, gain and error:
        continue;
    end
    
    % Estimate the linear prediction coefficients:
    a       = lpc(wframe,M);
    A(:,fi) = a;
    resid   = filter(a,1,wframe);
    %%
    gain    = sqrt(mean(resid.^2));
    g(fi)   = gain;
    ei      = resid / gain;
    e(from:to)  = e(from:to) + ei;
end

end