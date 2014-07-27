function [features] = calc_temporal_features(w,sr,gain_sr,gain_filterbank)

% Calculate gains (local root mean square):
powers      = w.^2;

% The smoothing filter should be longer than the gain-subsampling period:
s_len       = round(2*sr/gain_sr);
smoother    = hamming(s_len); smoother  = smoother / sum(smoother);
rms         = filter(smoother,1,powers);
rms         = rms / max(rms); % Normalize
gains       = resample(rms,gain_sr,sr);

[features,responses] = temporal_filterbank_features(gains,gain_filterbank);

end