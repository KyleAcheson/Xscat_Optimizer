function [I] = lsq_tfunc(weights, Ith, Iexp, Nts, Ntraj, Nq, CM, FLAGxfrac, FLAGexclude, ex_trajs)

if FLAGexclude == 1
    weights(ex_trajs) = 0;
end

weights(1:Ntraj) = weights(1:Ntraj) / sum(weights(1:Ntraj));

if FLAGxfrac == 1
    exfrac = weights(end);
else
    exfrac = 1;
end

I = zeros(Nq, Nts);

for ts=1:Nts
    for tr=1:Ntraj
        temp = weights(tr) * Ith(tr, 1:Nq, ts);
        I(1:Nq, ts) = I(1:Nq, ts) + temp.';
    end
end

I = I*(exfrac*100);
I = I - Iexp;
I = I .* CM;

end

