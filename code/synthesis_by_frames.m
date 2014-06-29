function [synth_w,synth_ge] = synthesis_by_frames(params)

if ~isfield(params,'A')
    error('You must provide a time series of filters params.A');
end
if ~isfield(params,'g')
    error('You must provide a time series of gains params.g');
end
if ~isfield(params,'excitation_source')
    params.excitation_source    = 'periods';
else
    if ~any(strcmp(params.excitation_source,{'periods','e'}))
        error(['Unsupported excitation source: ' params.excitation_source]);
    end
end
if ~isfield(params,params.excitation_source)
    error(['Missing excitation source: params.' params.excitation_source]);
end
if ~isfield(params,'preemph')
    params.preemph  = 1;
end
if ~isfield(params,'hoplen')
    params.hoplen   = 1024;
end
if ~isfield(params,'winlen')
    params.winlen   = 2*params.hoplen;
end

A           = params.A;
g           = params.g;

switch (params.excitation_source)
    case 'e'
        e   = params.e;
    case 'periods'
        e   = synthesize_excitation_from_periods(params.periods,params.hoplen);
end

L           = length(e);
n_frames    = length(g);

% % Prepare a smoother gain wave:
% gains       = reshape(repmat(g,params.hoplen,1),L,1);
% smoother    = hamming(100);
% smoother    = smoother / sqrt(mean(smoother.^2));
% gains       = filter(smoother,1,gains);

synth_ge    = zeros(L,1);
synth_w     = zeros(L,1);
window      = hamming(params.winlen);
window      = window / mean(window.^2);
for fi = 1:n_frames
    from    = (fi-1)*params.hoplen + 1;
    to      = from+params.winlen-1;
    if to > L
        to      = L;
        window  = hamming(to-from+1);
        window  = window / mean(window.^2);
    end
    
    eframe  = e(from:to);
    a       = A(:,fi);
    gain    = g(fi);
    geframe = gain * eframe;
%    gframe  = gains(from:to);
%    geframe = gframe .* eframe;
    sframe  = filter(1,a,geframe);
    wframe  = sframe .* window;
    
    synth_ge(from:to)   = synth_ge(from:to) + geframe;
    synth_w(from:to)    = synth_w(from:to) + wframe;
end

% De-emphasize the pre-emphasis that was done:
%synth_w     = filter(1,params.preemph,synth_w);

end