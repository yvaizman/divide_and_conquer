function [f0_mat,f0_comps] = mark_f0_in_spec(spec)

[spec_size,num_frames] = size(spec);

max_powers  = max(spec);
threshes    = 10e-4*max_powers;
over_thresh = (spec - (ones(spec_size,1)*threshes) > 0);
is_peak     = (...
    spec > [zeros(1,num_frames);spec(1:end-1,:)] &...
    spec > [spec(2:end,:);zeros(1,num_frames)]);

candidate   = over_thresh & is_peak;

f0_mat      = zeros(size(spec));
f0_comps    = zeros(1,num_frames);
for fi = 1:num_frames
    f0ind   = find(candidate(:,fi),1,'first');
    f0_comps(fi)        = f0ind;
    f0_mat(f0ind,fi)    = 1;
end

end