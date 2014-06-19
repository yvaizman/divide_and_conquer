% function [a,en,report] = burg_l1_minimization(xn,params)
%
% Find optimal linear prediction parameters, such that the prediction error
% signal has minimal L1-norm. Use a Burg-type algorithm, as described by
% Denoel and Solvay, 1985.
%
% Input:
% -----
% xn: (L x 1) vector. The observed signal x(n).
% params: strutc. Containing:
% params.M: positive integer. The order of the predictor.
%
% Output:
% ------
% a: ((M+1) x 1) vector. a(1) = 1. The fitted optimal FIR filter to get the
%   prediction error from the observed signal.
% en: (L x 1) vector. The residual prediction error signal e(n).
% ks: (1 x M) vector. The reflection coefficients of the lattice filter.
%
% ------------------------------------------------------------------------
% Written by Yonatan Vaizman, 2014.
function [a,en,ks] = burg_l1_minimization(xn,params)

M               = params.M;
%% Check for degenerate cases:
if var(xn) <= 0
    a           = [1;zeros(M,1)];
    en          = zeros(size(xn));
    ks          = poly2rc(a);
end

%% Prepare the order-recursive forward error and backward error containers,
%  for orders 1...M:
F               = repmat(xn,1,M);
B               = repmat(xn,1,M);
% For order 0 the forward error signal and the backward error signal are
% equal to the observed signal xn:
f               = xn;
b               = xn;

% The reflection coefficients:
ks              = zeros(1,M);

% Begin the order-recursive procedure:
for m = 1:M
    % Calculate all the ratios and their weights:
    q_f         = -f(2:end) ./ b(1:end-1); % -f_{m-1}[n]/b_{m-1}[n-1]
    w_f         = abs(b(1:end-1));
    q_b         = 1./q_f;
    w_b         = abs(f(2:end));
    
    q           = [q_f;q_b];
    w           = [w_f;w_b];
    
    total_w     = sum(w);
    half_w      = total_w/2;

    % Sort the ratios:
    [qs,order]  = sort(q,'ascend');
    ws          = w(order);
    cum_w       = cumsum(ws);
    
    % Look for the first value to be larger than half the total weight:
    ind         = find(cum_w>=half_w,1,'first');
    % This is the optimal value of the m'th reflection coefficient:
    k           = qs(ind);
    ks(m)       = k;
    
    % Calculate the new forward and backward prediction errors:
    F(2:end,m)  = f(2:end) + k*b(1:end-1);
    B(2:end,m)  = k*f(2:end) + b(1:end-1);
    f           = F(:,m);
    b           = B(:,m);
end

% The final resulted error signal:
en              = f;

% The final optimal filter, in prediction coefficients:
a               = rc2poly(ks)';

end
