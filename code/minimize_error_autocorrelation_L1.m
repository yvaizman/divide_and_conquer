% function [a,en,report] = minimize_error_autocorrelation_L1(xn,M,lags)
%
% Find optimal linear prediction parameters, such that the prediction error
% signal's autocorrelation has minimal L1-norm (over the desired lags).
%
% This function applies a partial-gradient-descent approach. In each
% iteration the update of the parameters is according to the gradient with
% respect to only some of the lags of the error-autocorrelation.
%
% Input:
% -----
% xn: (L x 1) vector. The observed signal x(n).
% M: positive integer. The order of the predictor.
% lags: vector of positive integers. The lags over which we want to
%   minimize the absolute autocorrelation of the prediction error signal.
%
% Output:
% ------
% a: ((M+1) x 1) vector. a(1) = 1. The fitted optimal FIR filter to get the
%   prediction error from the observed signal.
% en: (L x 1) vector. The residual prediction error signal e(n).
% report: struct. Containing:
% report.A: ((M+1) x T) matrix. History of the updated filter parameters
%   during the T iterations of the minimization procedure.
% report.deltanorms: (1 x T) vector. History of the L2 norms of the
%   loss-derivatives (used as delta updates to the parameters).
% report.ennorms: (2 x T) vector. History of the L1-norms (first row) and
%   L2-norms (second row) of the residual signal e(n).
% report.losses: (1 x T) vector. History of the partial loss function
%   values during the minimization procedure.
%
% ------------------------------------------------------------------------
% Written by Yonatan Vaizman, 2014.
function [a,en,report] = minimize_error_autocorrelation_L1(xn,M,lags)

maxlag          = max(lags);
ac              = autocorr(xn,M+maxlag);
R0              = toeplitz(ac(1:(M+maxlag+1)));
params          = struct('lags',lags,'ac',ac,'R0',R0);

aa              = zeros(M,1);
a               = [1;aa];
max_iter        = 7000;
etta            = 0.005;
A               = zeros(M+1,max_iter);
deltanorms      = zeros(1,max_iter);
ennorms         = zeros(2,max_iter);
losses          = zeros(1,max_iter);
% Start iterating:
for t = 1:max_iter
    en          = filter(a,1,xn);
    ennorms(t)  = norm(en,1);
    ennorms(2,t)= norm(en,2);

    lagi        = mod(t,length(lags)) + 1;
    params.lags = union(0,lags(lagi));
    [loss,delta]= error_autocorrelation_L1(aa,params);
    losses(t)   = loss;
    deltanorms(t)   = norm(delta);
    
    aa          = aa - etta * delta;
    a           = [1;aa];
    A(:,t)      = a;
end

report.A            = A;
report.deltanorms   = deltanorms;
report.ennorms      = ennorms;
report.losses       = losses;

end

function [loss,deriv]   = error_autocorrelation_L1(aa,params)

lags            = params.lags;
ac              = params.ac;
R0              = params.R0;

M               = length(aa);

loss            = 0;
deriv           = zeros(size(aa));

for k = lags
    % Calculate the autocorrelation vectors and matrices for lag k:
    plus_lags   = (k+1):(k+M);
    dplus       = ac(1 + plus_lags); % r(k+1);r(k+2)...;r(k+M)
    minus_lags  = (k-1):-1:(k-M);
    dminus      = ac(1 + abs(minus_lags)); % r(k-1);r(k-2)...;r(k-M)

    RkM         = R0(1:M,(k+1):(k+M));

    % Helper vector and matrix:
    d           = dplus + dminus;
    Rsym        = RkM + RkM';

    % Calculate the current excitation autocorrelation at lag k:
    % r_e(k)    = r_x(k) + d'*a + a'*RkM*a
    re          = ac(1+k) + d'*aa + aa'*RkM*aa;

    % Updaate the current loss:
    % J(a)      = |r_e(k)|
    loss_k      = abs(re);
    loss        = loss + loss_k;

    % Update the current derivative of the loss wrt a:
    % dJ/da     = sign(r_e(k))*(d + (RkM+RkM')*a)
    deriv_k     = sign(re) * (d + Rsym*aa);
    deriv       = deriv + deriv_k;
end

end