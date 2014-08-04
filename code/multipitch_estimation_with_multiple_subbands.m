% function [sal_mat] = multipitch_estimation_with_multiple_subbands(spec,params)
%
% Analyze a short time magnitude spectrum by separate subbands,
% and produce a saliency function from each subband, for multipitch
% estimation.
%
% Input:
% -----
% spec: (spec_size x 1) vector. The complex DFT,
% without the redundant components.
% params: struct. Containing:
% params.salience_of_what: string. Either 'f0' for salience functions of
% fundamental frequency or 'period' for salience functions of periods (in
% time samples).
% params.dft_freqs: vector of length spec_size of the frequencies (in Hz)
% fitting the DFT bins provided in the given spectrum.
%
% If analysing f0 must provide:
% params.min_fi: The index (in the given DFT spectrum) of the minimal f0
% candidate.
% params.max_fi: The index (in the given DFT spectrum) of the maximal f0
% candidate.
%
% If analysing period must provide:
% params.min_period: The minimal candidate period.
% params.max_period: The maximal candidate period.
%
% params.subbands: (n_subs x 2). For every analysis subband, the starting
% frequency index and the length (number of frequency components from the
% DFT) of the subband.
% params.compress: scalar. The power to compress magnitudes by.
% params.compress_in_freq: boolean. If true do the power compression in the
% frequency domain, otherwise do the compression in the time (lag) domain.
% ------------------------------------------------------------------------
% Written by Yonatan Vaizman, July 2014.
function [sal_mat,selected_f0s_max,selected_f0s_avr] = multipitch_estimation_with_multiple_subbands(spec,params)

spec_size   = size(spec,1);
switch (params.salience_of_what)
    case 'f0'
        n_cand      = params.max_fi - params.min_fi + 1; % Candidates for f0
    case 'period'
        n_cand      = params.max_period - params.min_period + 1; % Candidates for period
end
n_subs      = size(params.subbands,1);

% Should we compress the signal already in the given (frequency) domain:
if params.compress_in_freq
    spec    = abs(spec) .^ params.compress;
end

sal_mat     = zeros(n_cand,n_subs);
for subi = 1:n_subs
    % Subband boundaries:
    from    = params.subbands(subi,1);
    len     = params.subbands(subi,2);
    to      = min([spec_size,from+len-1]);
    subband = spec(from:to,1);
    if (length(subband)<len)
        subband     = [subband;zeros(len-length(subband),1)];
    end
    
    switch (params.salience_of_what)
        case 'f0'
            sub_sal = get_freq_salience_from_subband(subband,params);
        case 'period'
            sub_sal = get_period_salience_from_subband(subband,params);
    end
    
    sal_mat(:,subi)     = sub_sal;
end

%% Select top f0s candidates:
fi_vals             = params.min_fi:params.max_fi;
f0_vals             = params.dft_freqs(fi_vals);

max_sal             = max(sal_mat,[],2);
selected_f0s_max    = select_top_salience_candidates(f0_vals,max_sal,params.mix_size);

avr_sal             = mean(sal_mat,2);
selected_f0s_avr    = select_top_salience_candidates(f0_vals,avr_sal,params.mix_size);

end


function [time_sig] = transform_from_freq_to_time(freq_sig)
% Assume we have the no-DC, non-redundant components of DFT.
% Add a zero DC and the symmetric components:
sym                     = [0;freq_sig;conj(freq_sig((end-1):-1:1))];
time_sig                = real(ifft(sym));
end

function [freq_salience] = get_freq_salience_from_subband(subband,params)

% Get basic periodicity measure:
if params.compress_in_freq
    % Then subband is already compressed magnitude.
    freq_ac             = autocorr(subband,params.max_fi-1); % The lags are of basic frequency resolution of the DFT
else
    % Then we need to transform to time domain, compress and transform
    % back:
    time_sig            = transform_from_freq_to_time(subband);
    % Compress in the time domain:
    time_sig            = abs(time_sig) .^ params.compress;
    % Return to frequency domain to get the generalized
    % frequency-autocorrelation:
    freq_ac             = abs(fft(time_sig));
    % Remove the higher components:
    freq_ac             = freq_ac(1:(params.max_fi+1));
end
freq_ac                 = reshape(freq_ac,length(freq_ac),1);

% Get rid of period-multiples:
enhanced_freq_ac        = enhance_autocorrelation_function(freq_ac,params.max_multiple);

% Get the relevant indices:
freq_salience           = enhanced_freq_ac((params.min_fi):end);

end

function [period_salience] = get_period_salience_from_subband(subband,params)

if params.compress_in_freq
    % Then the subband was already compressed.
    % We now need to symmetrize it and transform to the time-lag domain:
    time_ac             = transform_from_freq_to_time(subband);
    % Remove the higher components:
    time_ac             = time_ac(1:(params.max_period+1));
else
    % Then first transform to time domain:
    time_sig            = transform_from_freq_to_time(subband);
    % Then compress:
    time_sig            = time_sig .^ params.compress;
    % And then calculate the time autocorrelation:
    time_ac             = autocorr(time_sig,params.max_period);
end

time_ac                 = reshape(time_ac,length(time_ac),1);

% Get rid of period-multiples:
enhanced_time_ac        = enhance_autocorrelation_function(time_ac,params.max_multiple);

% Get the relevant indices:
period_salience         = enhanced_time_ac((params.min_period+1):end);

end