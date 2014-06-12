function [a,en,A,deltanorms,ennorms,losses,a_opt] = iterative_minimization(xn,M,criterion,lags)

switch (criterion)
    case 'ac_lag1'
        loss_func   = @loss__excitation_autocorr_lag1;
    case 'acL1'
        maxlag      = max(lags);
        ac          = autocorr(xn,M+maxlag);
        R0          = toeplitz(ac(1:(M+maxlag+1)));
        params      = struct('lags',lags,'ac',ac,'R0',R0);
        loss_func   = @loss__excitation_abs_autocorr_lags;
    otherwise
        error(['!!! unsupported criterion for minimization: ' criterion]);
end


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
    params.lags = union([0],lags(lagi));
    [loss,delta]= loss_func(xn,aa,params);
    losses(t)   = loss;
    deltanorms(t)   = norm(delta);
    
    aa          = aa - etta * delta;
    a           = [1;aa];
    A(:,t)      = a;
end

a_opt           = A(:,end);

end

function [en,enac] = apply_filter_and_get_ac(xn,a)

en = filter(a,1,xn);
enac = autocorr(en,1000);

end

function [loss,deriv]   = loss__excitation_abs_autocorr_lags(xn,aa,params)

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

function [aa_opt,a_opt] = immediate_optimal_filter__excitation_autocorr_lags(xn,M,lags)
maxlag          = max(lags);
ac              = autocorr(xn,M+maxlag);
d               = zeros(M,1);
R               = zeros(M,M);

R0              = toeplitz(ac(1:(M+maxlag+1)));

for k = lags
    % Calculate the autocorrelation vectors and matrices for lag k:
    plus_lags   = (k+1):(k+M);
    dplus       = ac(1 + plus_lags); % r(k+1);r(k+2)...;r(k+M)
    minus_lags  = (1-k):(M-k);
    dminus      = ac(1 + abs(minus_lags)); % r(k-1);r(k-2)...;r(k-M)
    
    RkM         = R0(1:M,(k+1):(k+M));
    
    % Update the total ac vector and matrix:
    d           = d + dplus + dminus;
    R           = R + RkM + RkM';
end

invmat          = inv(R);
aa_opt          = - invmat * d;

a_opt           = [1;aa_opt];
end
function [aa_opt] = immediate_optimal_filter__excitation_autocorr_lag1(xn,M)
ac              = autocorr(xn,M+1);
d2              = ac(3:M+2); % r2;r3...;r_{M+1}
d0              = ac(1:M); % r0;r1...;r_{M-1}
R0              = toeplitz(ac(1:M+1));
R1M             = R0(1:end-1,2:end);

invmat          = inv(R1M + R1M');
aa_opt          = - invmat * (d2 + d0);
end

function [loss,deriv]   = loss__excitation_autocorr_lag1(xn,aa,param)
M               = length(aa);
ac              = autocorr(xn,M+1);
r1              = ac(2);
d2              = ac(3:M+2); % r2;r3...;r_{M+1}
d0              = ac(1:M); % r0;r1...;r_{M-1}
R0              = toeplitz(ac(1:M+1));
R1M             = R0(1:end-1,2:end);

loss            = r1 + d2'*aa + d0'*aa + aa'*R1M*aa;
deriv           = d2 + d0 + (R1M + R1M')*aa;

end