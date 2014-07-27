function [filterbank] = get_gamma_filterbank(filt_size,sr)

times_int       = 0:(filt_size-1); % In samples
times           = times_int / sr; % In seconds
[n_vals,b_vals] = get_gamma_filters_parameters();
bank_size       = length(n_vals);
filterbank      = zeros(bank_size,filt_size);

for fi = 1:bank_size
    n           = n_vals(fi);
    b           = b_vals(fi);
    gamma_filt  = times.^(n-1) .* exp(-b*times);
    gamma_filt  = gamma_filt / sqrt(sum(gamma_filt.^2));
    filterbank(fi,:)    = gamma_filt;
end

end

function [n_vals,b_vals] = get_gamma_filters_parameters()

n_vals      = [1.2*ones(8,1);1.5*ones(8,1);2*ones(8,1);3*ones(8,1)];
b_vals      = [...
    0.2;0.25;0.333;0.5;1;2;4;10;...
    0.5;0.625;0.833;1.25;2.5;5;10;25;...
    1;1.25;1.67;2.5;5;10;20;50;
    2;2.5;3.33;5;10;20;40;100;...
    ];

end