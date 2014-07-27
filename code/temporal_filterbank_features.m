function [features,responses] = temporal_filterbank_features(gains,filterbank)

L           = length(gains);
bank_size   = size(filterbank,1);
responses   = zeros(bank_size,L);

for fi = 1:bank_size
    responses(fi,:)     = filter(filterbank(fi,:),1,gains);
end

features    = max(responses,[],2);

end