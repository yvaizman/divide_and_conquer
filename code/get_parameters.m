function [params] = get_parameters(method,params)

switch (method)
    case 'esacf'
        params  = get_esacf_params(params);
    case 'multiband'
        params  = get_multiband_params(params);
end

end

function [params] = add_basic_params(params)
%params.sr = 22050;
params.max_multiple = 10;
params.compress = 0.5;
end

function [params] = get_multiband_params(params)
params.salience_of_what='f0';
params.min_fi=4;
params.max_fi=150;
params.subbands = [[2:100:1000]',200*ones(10,1)];
params.compress_in_freq = true;

params  = add_basic_params(params);
end

function [params] = get_esacf_params(params)
% The cutoffs for the two filters:
low_f       = 70; % Hz
mid_f       = 1000; % Hz
hi_f        = 10000; % Hz

low_w       = 2*low_f / params.sr;
mid_w       = 2*mid_f / params.sr;
hi_w        = 2*hi_f / params.sr;

% Design the filters:
[blo,alo]   = butter(2,[low_w,mid_w],'bandpass');
[bhi,ahi]   = butter(2,[mid_w,hi_w],'bandpass');

params.blo  = blo;
params.alo  = alo;
params.bhi  = bhi;
params.ahi  = ahi;

params.min_period   = 30;
params.max_period   = 300;

params  = add_basic_params(params);
end