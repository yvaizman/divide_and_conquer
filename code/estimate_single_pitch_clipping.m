function [f0,period,f0_ac,period_ac] = estimate_single_pitch_clipping(en,params)

if ~isfield(params,'sr')
    error('Must provide parameter params.sr, in units of Hz');
end
if ~isfield(params,'min_f')
    params.min_f    = 20; % Hz
end
if ~isfield(params,'max_f')
    params.max_f    = 8000; % Hz
end

params.max_lag      = round(params.sr / params.min_f);
params.min_lag      = round(params.sr / params.max_f);

% First try estimating without clipping:
[f0_ac,period_ac] = estimate_f0_from_autocorrelation(en,params);

%% Now try with the clipped signal:
max_abs     = max(abs(en));
thresh      = 0.3*max_abs;

clipped     = en;
clipped(en > -thresh & en < thresh)     = 0;
clipped(en > thresh)                    = 1;
clipped(en < -thresh)                   = -1;

[f0,period] = estimate_f0_from_autocorrelation(clipped,params);

end

function [f0,period] = estimate_f0_from_autocorrelation(en,params)

ac              = autocorr(en,params.max_lag+1);
[peaks,locs]    = findpeaks(ac);
if length(peaks) <= 0
    f0          = -1;
    period      = -1;
    return;
end

period      = locs(1)-1;
f0          = params.sr / period;

end