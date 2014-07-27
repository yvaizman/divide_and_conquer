function [spec2f0s] = get_harmonic_matrix(spec_size,min_comp,max_comp)

N_noDC  = spec_size - 1;
n_f0s   = max_comp - min_comp + 1;

H       = zeros(n_f0s,N_noDC);
for fi = 1:n_f0s
    f0  = min_comp - 1 + fi;
    H(fi,f0:f0:end) = 1;
end

sums        = sum(H,2);
H           = H ./ (sums*ones(1,N_noDC));
spec2f0s    = [zeros(n_f0s,1),H];


end