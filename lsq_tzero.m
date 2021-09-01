function [Ith] = lsq_tzero(exfrac, Ith, Iexp)

Ith = Ith *(exfrac*100);

Ith = Ith - Iexp;

end

