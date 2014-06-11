% function [S] = filters2spectra(A,B,g,N,w_freqs)
%
% Convert the time series of LTI filters to a representation of a
% spectrogram.
%
% Input:
% -----
% A: ((p+1) x T). The p-poles denominator coefficients of the LTI filter
% (first coefficient is expected to be 1) for every time frame 1...T
% B: (q x T). The q-zeros numerator coefficients of the LTI filter
% for every time frame 1...T
% g: (1 x T). The amplitude gain of every time frame 1...T.
% N: integer. If positive: How many frequency points should we sample of the
% spectrum of each frame. If negative: disregard and use parameter w_freqs.
% w_freqs: vector of normalized [0,2pi] frequencies in which to sample the
% frequency response.
%
% Output:
% ------
% S: (N x T) complex. A complex spectrum representation (frequency response)
% for every time frame 1...T.
%
% ---------------------------------------------------------------------
% Written by Yonatan Vaizman. Feb 2014.
function [S] = filters2spectra(A,B,g,N,w_freqs)

if N < 0
    n_freqs = length(w_freqs);
else
    n_freqs = N;
end
[q,T]       = size(B);
S           = zeros(n_freqs,T);
for t = 1:T
    a       = A(:,t);
    b       = B(:,t);
    gain    = g(t);
    % Get the frequency response of the filter:
    if N < 0
        freq_res    = freqz(b,a,w_freqs);
    else
        freq_res    = freqz(b,a,N);
    end
    % Apply the gain of the frame:
    freq_res    = gain*freq_res;
    S(:,t)      = freq_res;
end

end