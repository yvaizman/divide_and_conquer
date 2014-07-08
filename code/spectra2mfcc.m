% function [mfcc,mfs] = spectra2mfcc(S,mel_mat,dct_mat)
%
% Calculate Mel frequency spectra (MFS) and 
% Mel frequency cepstral coefficients (MFCC) from spectra.
%
% Input:
% -----
% S: (n_freqs x T) complex. Collection (possibly time series) of T spectra,
%   each describing the frequency content in uniformly sampled frequencies.
% mel_mat: (n_mel_bins x n_freqs). Mel scale triangular filters to
%   sum up energies in mel scaled frequency bands.
% dct_mat: (num_ceps_coeffs x n_mel_bins). Discrete cosine transform
%   filters to apply to log-mel-bins in order to get cepstrum.
%
% Output:
% ------
% mfcc: (num_ceps_coeff x T). The MFCCs of the spectra.
% mfs: (n_mel_bins x T). The logarithm of energies in mel bins.
%
% ------------------------------------------------------------------------
% Written by Yonatan Vaizman, 2014.
function [mfcc,mfs] = spectra2mfcc(S,mel_mat,dct_mat)

powers      = abs(S).^2;
mel         = mel_mat * powers;
mfs         = 10*log10(mel);
mfcc        = dct_mat * mfs;

end