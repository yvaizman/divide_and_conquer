function [] = evaluate_multipitch_estimation(specs,gt_octaves,params)

n_trials                    = 100;

params.smooth_cutoff_hz     = 250;
[params.spec_size,params.N] = size(specs);
params.NFFT                 = 2*(params.spec_size-1);
params.dft_freqs            = params.sr * (0:(params.spec_size-1)) / params.NFFT;

params                      = get_parameters('esacf',params);
params                      = get_parameters('multiband',params);

mix_inds_per_trial          = cell(1,n_trials);
scores_per_trial            = zeros(6,n_trials);
octupdn_trial               = zeros(2,n_trials);
scores1b_per_trial          = zeros(4,n_trials);

for ti = 1:n_trials
    [scores,mix_inds,octupdn_matches,scores1b]    = ...
        do_random_mix(specs,gt_octaves,params);
    mix_inds_per_trial{ti}  = mix_inds;
    scores_per_trial(:,ti)  = scores;
    octupdn_trial(:,ti)     = octupdn_matches;
    scores1b_per_trial(:,ti)= scores1b;
end

end

function [scores,mix_inds,octupdn_matches,scores1b] = do_random_mix(specs,gt_octaves,params)

mix_inds        = select_m_inds(params.N,params.mix_size)
gt              = gt_octaves(mix_inds);
mix_spec        = mean(specs(:,mix_inds),2);

[filter_mag,excite_mag] = separate_excitation_filter_with_logdft_smoothing(...
    mix_spec,params.smooth_cutoff_hz,params.sr,-1);
excite_spec     = excite_mag .* exp(1i*angle(mix_spec));
sym             = [excite_spec;conj(excite_spec((end-1):-1:2))];
time_sig        = ifft(sym);
time_sig        = real(time_sig);

% Estimate the pitches:
[period_salience,selected_periods,selected_f0s] = multipitch_estimation_with_esacp(time_sig,params);
[esacf_pp,esacf_rr] = evaluate_multipitch_estimation_score(gt,selected_f0s);

[sal_mat,selected_f0s_max,selected_f0s_avr] = multipitch_estimation_with_multiple_subbands(excite_spec,params);
[mbmax_pp,mbmax_rr,octup_matches,octdn_matches] = evaluate_multipitch_estimation_score(gt,selected_f0s_max);
[mbavr_pp,mbavr_rr] = evaluate_multipitch_estimation_score(gt,selected_f0s_avr);

params_1band    = params;
params_1band.subbands = [2,params.spec_size];
[sal_mat1b,selected_f0s_max1b] = multipitch_estimation_with_multiple_subbands(excite_spec,params_1band);
[mbmax_pp1b,mbmax_rr1b,octup_matches1b,octdn_matches1b] = evaluate_multipitch_estimation_score(gt,selected_f0s_max1b);


scores              = [esacf_pp;esacf_rr;mbmax_pp;mbmax_rr;mbavr_pp;mbavr_rr];
octupdn_matches     = [octup_matches;octdn_matches];
scores1b            = [mbmax_pp1b;mbmax_rr1b;octup_matches1b;octdn_matches1b];
end

function [precision,recall,octup_matches,octdn_matches] = evaluate_multipitch_estimation_score(gt_octaves,est_freqs_hz)

est_octaves     = log2(est_freqs_hz / 440);
oct_thresh      = 1/24; % quarter of a tone

[precision,recall] = estimation_score(gt_octaves,est_octaves,oct_thresh);
[pp1,rr1,octup_matches] = estimation_score(gt_octaves+1,est_octaves,oct_thresh);
[pp2,rr2,octdn_matches] = estimation_score(gt_octaves-1,est_octaves,oct_thresh);

end

function [precision,recall,num_matches] = estimation_score(gt_octaves,est_octaves,oct_thresh)
n_gt            = length(gt_octaves);
n_est           = length(est_octaves);
match_mat       = zeros(n_gt,n_est);
for gti = 1:n_gt
    gt          = gt_octaves(gti);
    diffs       = reshape(est_octaves - gt,1,n_est);
    matches     = abs(diffs) < oct_thresh;
    match_mat(gti,:)    = matches;
end

num_matches     = sum(match_mat(:));
% Precision is average (over the estimated values) of correctly estimated:
precision       = mean(max(match_mat,[],1));
% Recall is average (over ground truth values) of detected:
recall          = mean(max(match_mat,[],2));
end

function [inds] = select_m_inds(N,m)

inds            = ceil(N*rand(1,m));

end