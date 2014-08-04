function [selected] = select_top_salience_candidates(candidates,salience,m)

[peaks,locs]            = findpeaks(salience);
if (length(peaks) <= 0)
    locs                = 1:length(salience);
    peaks               = salience;
end
[vals,inds]             = sort(peaks,'descend');
num_found               = length(inds);
if num_found < m
    num_missing         = m - num_found;
    inds                = [inds(1)*ones(num_missing);inds];
end
top_inds                = locs(inds(1:m));
selected                = candidates(top_inds);

end