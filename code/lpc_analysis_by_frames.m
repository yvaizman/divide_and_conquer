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
% params: Struct, containing:
% params.sr: Sampling rate (in Hz). If a wave vector was given, this will be
%     regarded as the sampling rate of the vector. If a filename was given,
%     the loaded wave will be resampled to this desired sampling rate.
% params.M: Positive integer. Order of the filters (number of poles).
% params.do_preemph: Boolean. Should we do a pre-emphasis? 
%     If true: a high-pass FIR 1-zero filter is applied to the given signal
%     before the frame-by-frame analysis.
% params.hoplen: Positive integer. Number of samples to hop from frame to frame.
% params.winlen: Positive integer. Size of each frame.
% params.preemph: vector. The coefficients of the pre-emphasis FIR filter 
%    to be used (set to 1 if no pre-emphasis is desired).
% params.use_window: boolean. If true, use hamming window for each frame.
% params.criterion: string. The objective to minimize for the
%    prediction-error signal. Either 'L2' or 'acL1'.
% params.lags: vector of positive integers. If criterion 'acL1' is chosen,
%    this parameter specifies for which lags the error-signal's
%    autocorrelation function should be minimized (in absolute value).
% params.version: string. If creterion 'acL1' is chosen, this parameter
%   should be either 'nonsymmetric' (default), 'symmetric' for the
%   symmetric autocorrelation matrices version, or 'convex' for the convex
%   relaxation of the probelm.
%
% Output:
% ------
% results. struct, containing:
% results.A: ((M+1) x n_frames) matrix. The fitted all-pole filters. Each column
%    has the coefficients of a frame's fitted filter (direct form).
% results.g: (1 x n_frames) vector. The gain of every frame.
% results.e: (L x 1) vector. The prediction-error/excitation/residual signal.
% results.params: Struct with the parameters used.
% results.period: (1 x n_frames) vector. Estimated period (in samples) for
%   each frame, or -1 if no estimate is available.
% results.f0: (1 x n_frames) vector. Estimated fundamental frequency (in
%   Hz) for each frame, or -1 in no estimate is available.
%
% ------------------------------------------------------------------------
% Written by Yonatan Vaizman, 2014.
function [results] = lpc_analysis_by_frames(wav,params)

% Default values for parameters:
if ~isfield(params,'sr')
    params.sr = 22050;
end
if ~isfield(params,'M')
	params.M = 20;
end
if ~isfield(params,'do_preemph')
	params.do_preemph = true;
end
if ~isfield(params,'use_window')
    params.use_window = true;
end
if ~isfield(params,'hoplen')
	params.hoplen = 1024;
end
if ~isfield(params,'winlen')
    params.winlen = 2*params.hoplen;
end
if ~isfield(params,'preemph')
    params.preemph = [1 -0.9];
end
if ~isfield(params,'criterion')
    params.criterion = 'L2';
end
if ~isfield(params,'lags')
    params.lags = 0:50;
end

%% Parameters for pitch estimation:
yin_params.sr       = params.sr;
yin_params.wsize    = round(0.75*params.winlen);
yin_params.hop      = 64;

if ischar(wav)
    wav_file    = wav;
    [w,sr_orig] = wavread(wav_file);
    w           = mean(w,2);
    w           = resample(w,params.sr,sr_orig);
else
    w           = wav;
end

if params.do_preemph
    % Pre-emphasize:
    w               = filter(params.preemph,1,w);
else
    params.preemph  = 1;
end

L           = length(w);
n_frames    = floor(L/params.hoplen);

%% Prepare the silence threshold for this wave:
squares     = w.^2;
seg_dur     = 0.005; % In seconds;
seg_len     = round(seg_dur * params.sr);
smoother    = (1/seg_len)*ones(1,seg_len);
powers      = filter(smoother,1,squares);
max_power   = max(powers);
thresh      = 10^(-6)*max_power; % Silence threshold will be 60dB below max

%%
A           = [ones(1,n_frames);zeros(params.M,n_frames)];
g           = zeros(1,n_frames);
e           = zeros(L,1);

aperiod     =  ones(1,n_frames);
periods     = -ones(1,n_frames);

window      = hamming(params.winlen);
window      = window / mean(window.^2);
for fi = 1:n_frames
    if ~mod(fi,100)
        disp(['frame ' num2str(fi)]);
    end
    from    = (fi-1)*params.hoplen + 1;
    to      = from+params.winlen-1;
    if to > L
        to  = L;
        window  = hamming(to-from+1);
        window  = window / mean(window.^2);
    end
    frame   = w(from:to);
    if params.use_window
        wframe  = frame .* window;
    else
        wframe  = frame;
    end
    % Simple check if there is any power:
    epsilon = thresh;
    if var(wframe) <= epsilon
        % Leave the default filter, gain and error:
        continue;
    end
    
    % Estimate the linear prediction coefficients:
    switch (params.criterion)
        case 'L2'
            a       = lpc(wframe,params.M);
        case 'L1'
            a       = burg_l1_minimization(wframe,params);
%            a       = minimize_absolute_error(wframe,params);
        case 'acL1'
            a       = minimize_error_autocorrelation_L1(wframe,params);
        otherwise
            error(['Unsupported optimization criterion: ' params.criterion]);
    end
    A(:,fi) = a;
    
    % Apply the inverse filter to get the residual signal (the prediction
    % error signal):
    resid   = filter(a,1,wframe);
    %%
    gain    = sqrt(mean(resid.^2));
    gain    = max([gain,thresh]);
    g(fi)   = gain;
    ei      = resid;
%    ei      = ei / gain;
    e(from:to)  = e(from:to) + ei;
    
    %% Estimate the pitch from the prediction error signal of the frame:
    yinr            = yin(ei,yin_params);
    apinds          = ~isnan(yinr.ap);
    if sum(apinds) > 0
        aperiod(fi) = mean(yinr.ap(apinds));
    end
    perinds         = ~isnan(yinr.period);
    if sum(perinds) > 0
        periods(fi) = mean(yinr.period(perinds));
    end
    
end

f0              = -ones(1,n_frames);
pos_per         = periods > 0;
f0(pos_per)     = params.sr ./ periods(pos_per);

results         = params;

results.A       = A;
results.e       = e;
results.g       = g;
results.params  = params;
results.aperiod = aperiod;
results.periods = periods;
results.f0      = f0;

end