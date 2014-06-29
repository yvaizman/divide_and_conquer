% function [e] = synthesize_excitation_from_periods(periods,framelen)
%
% Synthesize an excitation signal.
%
% Input:
% -----
% periods: (1 x n_frames) vector. For each frame, the pitch period (in samples)
%   of the periodic excitation in this frame.
% framelen: positive integer. Length (in samples) of each frame.
%
% Output:
% ------
% e: ((n_frames*framelen) x 1) vector. The synthesized signal.
% ------------------------------------------------------------------------
% Written by Yonatan Vaizman, 2014.
function [e] = synthesize_excitation_from_periods(periods,framelen)

n_frames    = length(periods);
L           = framelen * n_frames;
e           = zeros(L,1);

leftover    = 0;
for fi = 1:n_frames
    period          = periods(fi);
    period          = round(period); % to create consistent lags between pulses
    if period <= 0
        frame       = randn(framelen,1);
    else
        first       = max([1,period - leftover]);
        inds        = round(first:period:framelen);
        last        = inds(end);
        leftover    = framelen - last;

        % Create the frame:
        frame       = zeros(framelen,1);
        frame(inds) = 1;
    end
    
    % Normalize frame to have unit power:
    gain            = sqrt(mean(frame.^2));
    frame           = frame / gain;
    
    % Update the frame:
    fr_start    = (fi-1)*framelen + 1;
    fr_end      = fi*framelen;
    e(fr_start:fr_end)  = frame;
end

end