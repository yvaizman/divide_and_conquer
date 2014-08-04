function [period_salience] = multipitch_estimation_with_esacp(spec,params)

spec_size   = length(spec);
NFFT        = (spec_size-1)*2;
freqs       = params.sr*(0:spec_size)/NFFT;

% Assume the spectrum in spec is already flattened.

% The cutoffs for the two filters:
low_f       = 70; % Hz
mid_f       = 1000; % Hz
low_band    = (freqs>low_f & freqs<mid_f);
hi_band     = (freqs>mid_f);

% Filter in the frequency domain:
low_spec                = spec;
low_spec(~low_band)     = 0;
hi_spec                 = spec;
hi_spec(~hi_band)       = 0;

hi_time                 = transform_non_redundant_spectrum_to_time(hi_spec);
% Half-wave rectify:
hi_time_halfwave        = hi_time;
hi_time_halfwave(hi_time_halfwave<0)    = 0;
% Got back to frequency domain (assuming the time signal is already windowed):
hi_dft_halfwave         = fft(hi_time_halfwave);
% Low-pass filter the rectified high band signal:
hi_spec_halfwave_lp     = hi_dft_halfwave(1:spec_size);
hi_spec_halfwave_lp(~low_band)          = 0;

%% Compress in the frequency domain and transform back to the time domain
%% to get the generalized autocorrelation functions:
low_comp                = abs(low_spec) .^ params.compress;
hi_comp                 = abs(hi_spec_halfwave_lp) .^ params.compress;

low_ac                  = transform_non_redundant_spectrum_to_time(low_comp);
hi_ac                   = transform_non_redundant_spectrum_to_time(hi_comp);

%% Summarized autocorrelation function (SACF):
sum_ac                  = low_ac + hi_ac;

%% Enhanced summarized autocorrelation function (ESACF):
period_salience         = enhance_autocorrelation_function(sum_ac,params.max_multiple);

end

function [time_sig] = transform_non_redundant_spectrum_to_time(spec)

sym         = [spec;conj(spec((end-1):-1:2))];
time_sig    = ifft(sym);
time_sig    = real(time_sig);

end