function [period_salience,selected_periods,selected_f0s] = multipitch_estimation_with_esacp(time_sig,params)

% Assume the spectrum in spec is already flattened.


%% Filter the signal:
low_sig     = filter(params.blo,params.alo,time_sig);
hi_sig      = filter(params.bhi,params.ahi,time_sig);

% The high band signal will now be half-wave rectified and low pass
% filtered:
hirec       = hi_sig;
hirec(hirec<0)  = 0;
hireclo     = filter(params.blo,params.alo,hirec);

%% Calculate the periodicity measure from both channels and combine:
lo_ac       = generalized_autocorrelation(low_sig,params.compress);
hi_ac       = generalized_autocorrelation(hireclo,params.compress);

sacf        = lo_ac + hi_ac;

period_salience     = enhance_autocorrelation_function(sacf,params.max_multiple);

% Select top from candidate periods:
period_vals         = params.min_period:params.max_period;
cand_period_sal     = period_salience(1+period_vals);
selected_periods    = select_top_salience_candidates(period_vals,cand_period_sal,params.mix_size);
selected_f0s        = params.sr ./ selected_periods;

end

function [genac] = generalized_autocorrelation(time_sig,compress)

freq_sig    = fft(time_sig);
freq_comp   = abs(freq_sig) .^ compress;
genac       = ifft(freq_comp);
genac       = real(genac);

end