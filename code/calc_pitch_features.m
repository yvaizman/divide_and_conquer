function [pitch_features] = calc_pitch_features(seg,yin_params)

r           = yin(seg,yin_params);
periodicity = 1 - r.ap;

[periodicity_peak]   = max(periodicity);

[max_power,ind]         = max(r.pwr);
period_peak             = r.period(ind);
reloct_peak             = r.f0(ind);
f0_peak                 = yin_params.sr / period_peak;

thresh              = 0.4*max_power;
is_periodic         = ~isnan(r.pwr) & r.pwr > thresh;

periods             = r.period(is_periodic);periods=periods(~isnan(periods));
relocts             = r.f0(is_periodic);relocts=relocts(~isnan(relocts));
f0s                 = yin_params.sr ./ periods;f0s=f0s(~isnan(relocts));

med_period          = median(periods);
med_reloct          = median(relocts);
med_f0              = median(f0s);
std_period          = std(periods);
std_reloct          = std(relocts);
std_f0              = std(f0s);

pitch_features      = [...
    periodicity_peak;period_peak;reloct_peak;f0_peak;...
    med_period;med_reloct;med_f0;...
    std_period;std_reloct;std_f0];

end