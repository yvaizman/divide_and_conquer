function [salience] = spec2f0_salience(absS_noDC,fs)

spec_size   = size(absS_noDC,1);
max_har     = 1000;
max_comp    = 300;

sal_mat     = get_salience_filter_mat(spec_size,max_comp,max_har,fs);

salience    = sal_mat * absS_noDC;

end

function [sal_mat] = get_salience_filter_mat(spec_size,max_comp,max_har,fs)

NFFT        = 2*spec_size;
dfreq       = fs / NFFT;

alpha       = 52; %Hz.
betta       = 320; %Hz.
sal_mat     = zeros(max_comp,spec_size);

for f0i = 1:max_comp
    f0      = dfreq*f0i;
    for m = 1:max_har
        hari    = f0i*m;
        if hari > spec_size
            break;
        end
        har     = dfreq*hari;
        gim     = (f0 + alpha) / (har + betta);
        sal_mat(f0i,hari)  = gim;
    end
end

end