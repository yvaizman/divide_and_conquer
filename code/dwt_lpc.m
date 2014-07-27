% function [a,en] = dwt_lpc(xn,params)
%
% Find linear prediction parameters, by optimizing partial coefficients
% for different scales (matching different frequency bands).
% This function utilizes the discrete wavelet transform (DWT) to analyze
% the given signal in different scales, and for each scale it dedicates
% some of the poles of the final filter.
%
% Input:
% -----
% xn: (L x 1) vector. The observed signal x(n).
% params: strutc. Containing:
% params.poles_per_level: (1 x N) vector of non-negative integers. This
%   vector indicates the desired depth of the analysis (N-1), and for each
%   level, how many poles we want to fit according to the detail signal of that scale.
%   The last number indicates the number of poles to fit to the
%   approximation signal of the last level.
% params.wname: string. Wavelet name for the standard wavelet to use.
% params.calc_error: boolean. Should we also calculate the residual signal
%   after filtering with the inverse of the estimated all-pole filter.
%
% Output:
% ------
% a: ((M+1) x 1) vector. a(1) = 1. The fitted optimal FIR filter to get the
%   prediction error from the observed signal.
%   (M = sum(params.poles_per_level)).
% en: (L x 1) vector. The residual prediction error signal e(n).
%   In case params.calc_error is false, en will be simply set to -1.
%
% ------------------------------------------------------------------------
% Written by Yonatan Vaizman, 2014.
function [a,poles_per_level,filt_per_level,origfilt_per_level,gain_per_level,en] = dwt_lpc(xn,params)

if ~isfield(params,'poles_per_level')
    params.poles_per_level  = 2*ones(1,8);
end
if ~isfield(params,'wname')
    params.wname            = 'db45';
end
if ~isfield(params,'calc_error')
    params.calc_error       = false;
end

N               = length(params.poles_per_level)-1;
if any(params.poles_per_level < 0)
    error('number of poles per each level must all be non-negative');
end
M               = sum(params.poles_per_level);
poles           = [];%zeros(1,M);
[Lo_D,Hi_D,Lo_R,Hi_R] = wfilters(params.wname);

[ydwt,dec]      = wavedec(xn,N,params.wname);
poles_per_level = cell(1,N+1);
filt_per_level  = cell(1,N+1);
gain_per_level  = zeros(1,N+1);
origfilt_per_level  = cell(1,N+1);

% Add poles for each level's detail coefficients:
pole_i          = 0;
current_scale_poles = poles;
for level = 1:N
    detail      = detcoef(ydwt,dec,level);
    m           = params.poles_per_level(level);
    
    if params.independent_scales
        level_signal        = detail;
    else
        % Adjust the poles collected so far to fit the current scale:
        current_scale_poles = current_scale_poles.^2;
        % Apply the "so-far" inverse filter, to get the residual of the signal
        % from the details of this level:
        current_scale_filter    = poly(current_scale_poles);
        level_signal        = filter(current_scale_filter,1,detail);
    end
    
    % Calculate additional partial filter for this level's signal:
    [poles_lev,raw_poles_lev,gain,a_lev]  = get_poles_from_level(level_signal,level,m,true,params);

    poles_per_level{level}  = poles_lev;
    filt_per_level{level}   = real(poly(poles_lev));
    gain_per_level(level)   = gain;
    origfilt_per_level{level}   = a_lev;
    
%    poles(pole_i+(1:m)) = poles_lev;
    poles                   = [poles;poles_lev];
    if ~params.independent_scales
        current_scale_poles     = [current_scale_poles;raw_poles_lev];
    end
    pole_i                  = pole_i + m;
end

% Calculate the partial filter for the last level's approximation:
approx                      = appcoef(ydwt,dec,params.wname);
if params.independent_scales
    last_signal             = approx;
else
    % Apply the current occumulated filter (adjusted for the last scale):
    % (no need to adjust poles again, since approx is the same scale as the
    % latest detail signal)
    last_scale_filter       = poly(current_scale_poles);
    last_signal             = filter(last_scale_filter,1,approx);
end
m                           = params.poles_per_level(end);
[poles_lev,raw_poles_lev,gain,a_lev]      = get_poles_from_level(last_signal,N,m,false,params);

poles_per_level{end}        = poles_lev;
filt_per_level{end}         = real(poly(poles_lev));
gain_per_level(end)         = gain;
origfilt_per_level{end}     = a_lev;

%poles((pole_i+1):end)   = poles_lev;
poles                       = [poles;poles_lev];

% Get the final filter coefficients in direct form:
a                           = poly(poles);

% The residual (prediction error) signal:
if params.calc_error
    en                      = filter(a,1,xn);
else
    en                      = -1;
end

end

function [poles_lev,raw_poles_lev,gain,a_lev] = get_poles_from_level(sig_lev,level,m,is_detail,params)

[a_lev,ep]      = lpc(sig_lev,m);
gain            = sqrt(ep);
raw_poles_lev   = roots(a_lev);
% First, work only with one copy of complex-conjugate pairs:
poles_lev       = raw_poles_lev(imag(raw_poles_lev)>=0);

% Adjust the poles to their appropriate position in the original scale:
% 1. Reverse the last 2-downsample:
poles_lev       = sqrt(poles_lev);
% 2. If this is detail coefficients, they represent the high-frequency half
% of the spectrum:
if is_detail
    radii       = abs(poles_lev);
    angles      = pi - angle(poles_lev);
%     signs   = sign(angles);
%     signs(signs==0) = 1;
%     add     = (pi/2)*signs;
%     angles  = angles + add;
    poles_lev   = radii.*exp(1i*angles);
end
% 3. Reverse the [level-1] 2-downsamplings (with low passes between them):
mult            = 2^(level-1);
exponent        = inv(mult);
poles_lev       = poles_lev.^exponent;

% Handle the poles that are not on the real axis, and make sure they have a
% complex conjugate:
epsilon         = 10e-4;
effectively_real    = abs(imag(poles_lev)) < epsilon;
poles_lev(effectively_real) = real(poles_lev(effectively_real));
inds            = find(imag(poles_lev));
additional      = conj(poles_lev(inds));
poles_lev       = [poles_lev;additional];

end
