function [Ith_bin, Tth_bin, Iexp, TE, CM, T0_exp] = match_signal(Tth, dt, Ith, Iexp, Texp, T0, T0_exp, CM, Confidence_Tol, qmin, qmax, Ntraj, Nq, FLAGtdelay)

% MATCHES EXPERIMENTAL AND THEORY TEMPORAL AXES - FOR SIGNAL AND CONF. MAT

% INPUTS:
% Tth - theory time vector
% dt - theory dt
% Ith - theory signal (every trajectory)
% Iexp - experimental signal
% Texp - experimental time delays
% T0 - theory T0 to centre on
% T0_exp - experimental time delay selected as time zero to fit to T0
% CM - confidence matrix
% Confidence_Tol - threshold to remove data points from fit
% qmin/ qmax - min/ max q index in experiment
% Ntraj - number trajectories
% Nq - number points in q space
% FLAGtdelay - 0 = no binning, just sample time delays, 1 = binning of
% theory data according to diff(Texp)

% OUTPUTS:
% Ith_bin - theory signal (every trajectory) at matched exp. time delays
% Tth_bin - new theory time vector matched to exp. delays
% Iexp - experimental signal centered on time zero
% TE - selected experimental time delays
% CM - new confidence matrix, matching dimensions of Ith_bin and Iexp
% T0_exp - time delay on which T0 is centered in exp


T0th = find(Tth == T0); % T0 index in theory
Tlen_Th = floor(length(Tth)*dt); % time duration convoluted theory
pad_back = Tth(find(Tth == T0)) - Tth(1);

T0v(1:length(Texp)) = T0_exp;
Tfv(1:length(Texp)) = T0_exp + Tlen_Th; % duration of convoluted signal added

[~, T0_closest] = min(abs(T0v - Texp)); % closest T0 value in experiment vs requested T
[~, Tf_closest] = min(abs(Tfv - Texp)); % closest Tf value in experiment vs requested T

while ceil(Texp(Tf_closest) - Texp(T0_closest)) < Tlen_Th
    Tf_closest = Tf_closest + 1;
end

T0_exp = Texp(T0_closest); % experimental T0 centre

Tstart_exp = T0_exp - pad_back;
T0v(1:length(Texp)) = Tstart_exp;
[~, Ts_closest] = min(abs(T0v - Texp)); % closest T value in experiment vs requested T
Tstart_exp = Texp(Ts_closest);
Tf_exp = Texp(Tf_closest); % experimental Tf - end time
Tlen_exp = ceil(Tf_exp - Tstart_exp);

while Tlen_Th < Tlen_exp  % cut experiment time grid so it is always equal or shorter than theory
    Tf_closest = Tf_closest - 1;
    Tf_exp = Texp(Tf_closest);
    Tlen_exp = ceil(Tf_exp - Tstart_exp);
end
   
tt = length(Tth);
while Tlen_Th > Tlen_exp
    Tth = Tth(1:tt);
    tt = tt - 1;
    Tlen_Th = Tlen_Th - dt;
end

TE = Texp(Ts_closest:Tf_closest);
Iexp = Iexp(qmin:qmax, Ts_closest:Tf_closest);
CM = CM(qmin:qmax, Ts_closest:Tf_closest);
CM = CM./max(CM(:)); % normalise to max value
CM(CM < Confidence_Tol) = 0; % values below this will be set to zero weight


disp(['Observed values with a confidence less than: ', num2str(Confidence_Tol), ' are set to 0 weighting.'])
disp(['Length of theory time vector: ', num2str(Tlen_Th) , '. Length of experimental time vector, ' num2str(Tlen_exp)])      
disp(['EXPERIMENTAL T0 CENTRED ON: ', num2str(T0_exp), '. Signal included back from, ' num2str(TE(1))]);

if FLAGtdelay == 0
    bins = ceil([diff(TE)]);
    bin_len = floor(bins/dt);
    bin_len = [1 bin_len];
    bin_inds = cumsum(bin_len);

    Nts = length(bin_inds);
    Ith_bin = zeros(Ntraj, Nq, Nts);
    for ts=1:Nts
        Ith_bin(:, :, ts) = Ith(:, :, bin_inds(ts));
        Tth_bin(ts) = Tth(bin_inds(ts));
    end

elseif FLAGtdelay == 1
    
bins = ceil([diff(TE)]);
bin_len = floor(bins/dt); % indexes for each bin in time (theory)
bin_len = [1 bin_len];
bin_len = cumsum(bin_len);

lb_final = bin_len(end-1)+1; % HACK 
ub_final = bin_len(end)+1;

bin_lb = [bin_len(1:end-1)];
bin_ub = [bin_len(2:end)];
timestep = (bin_ub - bin_lb);

Ith_bin = zeros(Ntraj, Nq,length(bin_lb));

for i=1:length(bin_len)-1
    lb = bin_lb(i);
    ub = bin_ub(i);
    Ith_bin(:,:,i) = sum(Ith(:,:,lb:ub), 3)./timestep(i);
    times(i,1) = Tth(lb);
    times(i,2) = Tth(ub);
end

Tth_bin = times(:,1);
Iexp = Iexp(:, 1:end-1);
TE = TE(1:end-1);
CM = CM(:, 1:end-1);
end

end
