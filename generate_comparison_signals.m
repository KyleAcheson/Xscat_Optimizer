clear all
close all

FLAGelec = 1; % 0 = x-ray, 1 = electron
FLAGinel = 1; % include inelastic corrections
FLAGsignal = 0; % dI/I = 0
FLAGconfmat = 0; % 0 = no conf matrix, 1 = use
FLAGtdelay = 0; % 0 = no binning (just sample convolved signal), 1 = bin theory data

xfrac = 0.0338; % excitation fraction
pulse = 100; % pulse FWHM for convolution
atmnum = [6 16 16]; % atomic numbers
kin = 1.8751e+03; % kin ued
dt = 0.5; % time step in theory
time = 0:dt:1000; % time vec
T0_exp = 80; % T0 in experiment
T0 = 0; % T0 in theory
Confidence_Tol = 0; % confidence tol - 0 = include all

fpath_exp = 'FixedExperimental.mat'; % path to experimental data
fpath_traj = 'Filtered_Trajs_Final.mat'; % path to theory data

[Texp, Iexp, q_exp, CM] = load_experiment(fpath_exp, FLAGconfmat); % need experimental time vec, signal and q range
[Q, multiplicity] = load_trajectories(fpath_traj); % need geometries and spin of trajectories

[natom, ~, Ntraj, Nts] = size(Q);

qmin = min(q_exp); qmax = max(q_exp);
[~, qmin_closest] = min(abs(qmin - q_exp)); % closest q values in experiment vs requested q range
[~, qmax_closest] = min(abs(qmax - q_exp));
q_exp = q_exp(qmin_closest:qmax_closest); % selected q range for experiment (inv Angstrom)
Nq = length(q_exp);

[pdW] = theory_signal(Q, kin, atmnum, q_exp, time, FLAGelec, FLAGinel, FLAGsignal);

if pulse ~= 0
    [conv_signal, tconv, sigma] = convolute(pdW, time, pulse, 0);
else
    conv_signal = pdW; tconv = time;
end

conv_signal = conv_signal * (100 * xfrac);

Ntc = length(tconv);
conv_signal = reshape(conv_signal, [1, Nq, Ntc]);

[binned_signal, Tth_bin, Iexp, TE, CM, T0_exp] = match_signal(tconv, dt, conv_signal, Iexp, Texp, T0, T0_exp, CM, Confidence_Tol, qmin_closest, qmax_closest, 1, Nq, FLAGtdelay);

binned_signal = squeeze(binned_signal);
conv_signal = squeeze(conv_signal);

figure
[QQ, TT] = meshgrid(q_exp, tconv);
mesh(QQ, TT, conv_signal.');
axis tight
xlim([0 8]);
view(0, 90);
%title('Pulse = 230 fs')
ylabel(['Time (fs)'], 'interpreter', 'latex')
xlabel(['s (\AA', '$ ^{-1} $', ')'], 'interpreter', 'latex')





