function [enhanced_ac] = enhance_autocorrelation_function(ac,max_multiple)

ac(ac < 0)              = 0;
for m = 2:max_multiple
    % Add zeros:
    expand              = upsample(ac,m);
    % Linear interpolation:
    expand              = filter([0.5,1,0.5],1,expand);
    % Get rid of the delay of the smoother:
    expand              = expand(2:end);
    % Remove expanded version from the autocorrelation function:
    ac                  = ac - expand(1:length(ac));
    ac(ac<0)            = 0;
end

enhanced_ac             = ac;

end
