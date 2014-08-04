function [filter_mag,excite_mag,filter_mfs,excite_mfs] = separate_excitation_filter_with_logdft_smoothing(S,freq_period_cutoff_hz,sr,mel_mat)

dc              = S(1,:);
S               = S(2:end,:);
time_winlen     = 2*(size(S,1)); % (time samples). Window len of time window that generated this DFT
freq_resolution = sr/time_winlen; % (Hz). Frequency hop from component to the next in this DFT
freq_period_cut = freq_period_cutoff_hz / freq_resolution; % (freq components)

smoother_len    = 50;
w_cutoff        = 1/freq_period_cut;
smoother        = fir1(smoother_len,w_cutoff);
delay           = round(smoother_len/2)-1;

logabsS         = log10(abs(S));
stretched       = [ones(delay,1)*logabsS(1,:);logabsS;];
smoothed        = filtfilt(smoother,1,stretched);
filter_log_spec = smoothed((delay+1):end,:);
excite_log_spec = logabsS - filter_log_spec;

% Add back the DC component:
filter_log_dc   = filter_log_spec(1,:);
filter_log_spec = [filter_log_dc;filter_log_spec];
excite_log_dc   = log10(abs(dc));
excite_log_spec = [excite_log_dc;excite_log_spec];

filter_mag      = 10.^(filter_log_spec);
excite_mag      = 10.^(excite_log_spec);

if length(mel_mat) > 1
    filter_mfs  = 20*log10(mel_mat*filter_mag);
    excite_mfs  = 20*log10(mel_mat*excite_mag);
else
    filter_mfs  = -1;
    excite_mfs  = -1;
end

end