% function [a,en,report] = minimize_error_autocorrelation_L1(xn,params)
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
% params: strutc. Containing:
% params.M: positive integer. The order of the predictor.
% params.lags: vector of positive integers. The lags over which we want to
%   minimize the absolute autocorrelation of the prediction error signal.
% params.version: string. Either 'nonsymmetric', 'symmetric' or 'convex'.
%   If 'symmetric' is selected, use the symmetric
%   (not necessarily PSD) versions of the autocorrelation matrices. 
%   This version of the problem is equivalent to the original, 
%   once we assume that we are dealing with autocorrelations and that they 
%   are symmetric. In this version,
%   instead of R^l we use 0.5(R^l + R^{-l})=0.5(R^l + R^l^T).
%   If 'convex' is selected, use the convex relaxation of the
%   optimization problem, in which the partial autocorrelation matrices are
%   projected to the PSD code (by symmetrizing them, eigenvalue
%   decomposing, zeroing the negative eigenvalues, and re-composing).
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
function [a,en,report] = minimize_error_autocorrelation_L1(xn,params)

if ~isfield(params,'version')
    params.version  = 'nonsymmetric';
end
if ~isfield(params,'max_iter')
    params.max_iter         = 200;
end

M               = params.M;
maxlag          = max(params.lags);
power           = mean(xn.^2);
params.ac       = power*autocorr(xn,M+maxlag);
params.R0       = toeplitz(params.ac(1:(M+maxlag+1)));
params.lags     = union(0,params.lags);
lag0_ind        = find(params.lags==0);
[params]        = prepare_ac_matrices_and_vectors(params);

% Initialize the filter:
switch params.version
    case 'nonsymmetric'
        init_params             = params;
        init_params.version     = 'convex';
        init_params.max_iter    = 70;
        [a,init_e,init_report]  = minimize_error_autocorrelation_L1(xn,init_params);
        aa                      = a(2:end);
        etta                    = 0.01;
    case 'symmetric'
        init_params             = params;
        init_params.version     = 'convex';
        init_params.max_iter    = 70;
        [a,init_e,init_report]  = minimize_error_autocorrelation_L1(xn,init_params);
        aa                      = a(2:end);
        [params]                = use_symmetric_ac_matrices(params);
        etta                    = 0.01;
    case 'convex'
        [params]                = project_ac_matrices_to_psd_cone(params);
        aa                      = zeros(M,1);
        a                       = [1;aa];
        etta                    = 0.005;
end

A               = zeros(M+1,params.max_iter);
deltanorms      = zeros(1,params.max_iter);
ennorms         = zeros(2,params.max_iter);
losses          = zeros(1,params.max_iter);
% Start iterating:
for t = 1:params.max_iter
%     en          = filter(a,1,xn);
%     ennorms(t)  = norm(en,1);
%     ennorms(2,t)= norm(en,2);

    lagi        = mod(t,length(params.lags)) + 1;
    % Set the lags to focus on in this iteration:
    params.partial_lagis    = 1:length(params.lags); %union(lag0_ind,lagi); 
    [loss,delta]= error_autocorrelation_L1(aa,params);
    losses(t)   = loss;
    deltanorms(t)   = norm(delta);
    
    aa          = aa - etta * delta;
    a           = [1;aa];
    A(:,t)      = a;
end

en                  = filter(a,1,xn);
report.A            = A;
report.deltanorms   = deltanorms;
report.ennorms      = ennorms;
report.losses       = losses;
report.used_params  = params;

end

function [params] = use_symmetric_ac_matrices(params)

params.RkM_per_lagi__nonSymmetric   = params.RkM_per_lagi;
for lagi = 1:length(params.RkM_per_lagi)
    RkM             = params.RkM_per_lagi{lagi};
    Rsym            = 0.5 * (RkM + RkM');
    params.RkM_per_lagi{lagi}   = Rsym;
end

end

function [params] = project_ac_matrices_to_psd_cone(params)

params.RkM_per_lagi__original   = params.RkM_per_lagi;
epsilon             = 0; %10e-6;

for lagi = 1:length(params.RkM_per_lagi)
    RkM             = params.RkM_per_lagi{lagi};
    Rsym            = 0.5 * (RkM + RkM');
    [V,D]           = eig(Rsym);
    D(D<0)          = epsilon;
    Rpsd            = V*D*V';
    params.RkM_per_lagi{lagi}   = Rpsd;
end

end

function [params] = prepare_ac_matrices_and_vectors(params)

M                       = params.M;

params.RkM_per_lagi     = cell(1,length(params.lags));
params.Rsym_per_lagi    = cell(1,length(params.lags));
params.d_per_lagi       = cell(1,length(params.lags));

for lagi = 1:length(params.lags)
    k           = params.lags(lagi);
    % Calculate the autocorrelation vectors and matrices for lag k:
    plus_lags   = (k+1):(k+M);
    dplus       = params.ac(1 + plus_lags); % r(k+1);r(k+2)...;r(k+M)
    minus_lags  = (k-1):-1:(k-M);
    dminus      = params.ac(1 + abs(minus_lags)); % r(k-1);r(k-2)...;r(k-M)

    RkM         = params.R0(1:M,(k+1):(k+M));

    % Helper vector and matrix:
    d           = dplus + dminus;
    Rsym        = RkM + RkM';
    
    params.RkM_per_lagi{lagi}   = RkM;
    params.Rsym_per_lagi{lagi}  = Rsym;
    params.d_per_lagi{lagi}     = d;
end

end

function [loss,deriv] = error_autocorrelation_L1(aa,params)

loss            = 0;
deriv           = zeros(size(aa));

for lagi = params.partial_lagis
    k           = params.lags(lagi);
    
    % Get the autocorrelation vectors and matrices for lag k:
    RkM         = params.RkM_per_lagi{lagi};
    Rsym        = params.Rsym_per_lagi{lagi};
    d           = params.d_per_lagi{lagi};

    % Calculate the current excitation autocorrelation at lag k:
    % r_e(k)    = r_x(k) + d'*a + a'*RkM*a
    re          = params.ac(1+k) + d'*aa + aa'*RkM*aa;

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