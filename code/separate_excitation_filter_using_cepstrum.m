function [lowceps_S,hiceps_S,lowceps_mfs,hiceps_mfs,lowceps,hiceps] = separate_excitation_filter_using_cepstrum(S,num_lowceps,mel_mat,do_real)

spec_size   = size(S,1);
symS        = [S;conj(S(end:-1:2,:))];
if do_real
    symS    = abs(symS);
end
logS        = log(symS);
cepstra     = ifft(logS);

lowceps     = cepstra; lowceps((num_lowceps+1):end,:)=0;
hiceps      = cepstra; hiceps(1:num_lowceps,:)=0;

lowceps_S   = exp(fft(lowceps));
hiceps_S    = exp(fft(hiceps));

lowceps_S   = lowceps_S(1:spec_size,:);
hiceps_S    = hiceps_S(1:spec_size,:);


lowceps_mfs = 20*log10(mel_mat*(abs(lowceps_S).^2));
hiceps_mfs  = 20*log10(mel_mat*(abs(hiceps_S).^2));

end
