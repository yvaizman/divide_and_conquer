% function [a,en,report] = minimize_absolute_error(xn,params)
%
% Find optimal linear prediction parameters, such that the prediction error
% signal has minimal L1-norm.
%
% Input:
% -----
% xn: (L x 1) vector. The observed signal x(n).
% params: strutc. Containing:
% params.M: positive integer. The order of the predictor.
% params.method: string. Either 'gradient'.
%   If 'gradient' is selected, a partial-gradient descent method is used,
%   which is practically a form of sign-LMS.
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
function [a,en,report] = minimize_absolute_error(xn,params)

M               = params.M;
L               = length(xn);
block_size      = 10;

if ~isfield(params,'method')
    params.method   = 'gradient';
end
if ~isfield(params,'max_iter')
    params.max_iter         = round(10*L/block_size);
end

%% Prepare the examples:
X   = zeros(M+1,L-M);
for row = 1:(M+1)
    X(row,:)    = xn((M+2-row):(L+1-row));
end

% Initialize the filter:
%aa              = zeros(M,1);
%a               = [1;aa];
al2             = lpc(xn,M)';
a               = al2;
aa              = a(2:end);
etta            = 0.001;

en              = filter(a,1,xn);

A               = zeros(M+1,params.max_iter);
deltanorms      = zeros(1,params.max_iter);
ennorms         = zeros(2,params.max_iter);
et              = zeros(block_size,params.max_iter);
% Start iterating:
for t = 1:params.max_iter
    % Take a block (mini-batch) of examples:
    inds        = (block_size*(t-1)):(block_size*t-1);
    inds        = mod(inds,L-M) + 1;
    n_inds      = M + inds;
    
    Xt          = X(:,inds);
    
    % Get the prediction error for this sample:
    E           = a' * Xt;
    en(n_inds)  = E';
    et(:,t)     = E';
    
%    en              = filter(a,1,xn);
    ennorms(1,t)    = norm(en,1);
    ennorms(2,t)    = norm(en,2);
    
    % Update the filter:
    weights     = sign(E);
    gradient    = Xt(2:end,:)*weights';
    deltanorms(t) = norm(gradient);
    aa          = aa - etta*gradient;
    a           = [1;aa];
    A(:,t)      = a;
end

losses              = sum(abs(et));

en                  = filter(a,1,xn);
report.A            = A;
report.deltanorms   = deltanorms;
report.ennorms      = ennorms;
report.losses       = losses;
report.used_params  = params;

end
