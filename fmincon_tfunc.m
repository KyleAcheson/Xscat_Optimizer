function [I] = fmincon_tfunc(weights, Ith, Iexp, Nts, Ntraj, Nq, FLAGxfrac);

weights(1:Ntraj) = weights(1:Ntraj) / sum(weights(1:Ntraj));

if FLAGxfrac == 1
    exfrac = weights(end);
else
    exfrac = 1;
end

I = zeros(Nq, Nts);

for ts=1:Nts
    for tr=1:Ntraj
        A = weights(tr) * exfrac * Ith(tr, 1:Nq, ts);
        I(1:Nq, ts) = I(1:Nq, ts) + A'; % weights(tr) * Ith(tr, 1:Nq, ts);
    end
end

I = I - Iexp;
I = sum(sum(I.^2));

end

