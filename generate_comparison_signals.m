clear all
% close all

FLAGopt = 1; % 0 = use equally weighted trajectories (no opt data), 1 = use outputs from optimisation
FLAGelec = 1; % 0 = x-ray, 1 = electron
FLAGinel = 1; % include inelastic corrections
FLAGsignal = 0; % dI/I = 0, dsM = 1
FLAGconfmat = 0; % 0 = no conf matrix, 1 = use
FLAGtdelay = 0; % 0 = no binning (just sample convolved signal), 1 = bin theory data
FLAGqextend = 0; % 0 = stick to original experimental q vector, 1 = extend q vector to max_q


if FLAGopt == 1
    load('Optimised_OUTPUT.mat', 'weights') % file to load weights from if using
end

xfrac = 0.0338; % excitation fraction
pulses = [0 50 80 100 150 230]; % pulse FWHM for convolution
atmnum = [6 16 16]; % atomic numbers
kin = 1.8751e+03; % kin ued
max_q = 12;
dt = 0.5; % time step in theory
time = 0:dt:1000; % time vec
T0_exp = 80; % T0 in experiment
T0 = 0; % T0 in theory
Confidence_Tol = 0; % confidence tol - 0 = include all

rsignal = {};
bsignal = {};
times = {};
TEs= {};
qs = {};

for p=1:length(pulses);

    pulse = pulses(p);

fpath_exp = '/home/kyle/2TB_HDD/OPTDATA/INPUTS/FixedExperimental.mat'; % path to experimental data
fpath_traj = '/home/kyle/2TB_HDD/OPTDATA/INPUTS/Filtered_Trajs_Final.mat'; % path to theory data
fout = 'dsMsignals'

[Texp, Iexp, q_exp, CM] = load_experiment(fpath_exp, FLAGconfmat); % need experimental time vec, signal and q range
[Q, multiplicity] = load_trajectories(fpath_traj); % need geometries and spin of trajectories

[natom, ~, Ntraj, Nts] = size(Q);

qmin = min(q_exp); qmax = max(q_exp);
[~, qmin_closest] = min(abs(qmin - q_exp)); % closest q values in experiment vs requested q range
[~, qmax_closest] = min(abs(qmax - q_exp));
q_exp = q_exp(qmin_closest:qmax_closest); % selected q range for experiment (inv Angstrom)
dq = q_exp(2)-q_exp(1);

if FLAGqextend == 1
    q_exp = qmin:dq:max_q;
end

Nq = length(q_exp);

if FLAGopt == 0 % do not use optimised weights - assume equal weighting
    [pdW] = theory_signal(Q, kin, atmnum, q_exp, time, FLAGelec, FLAGinel, FLAGsignal);

    if pulse ~= 0
        [conv_signal, tconv, sigma] = convolute(pdW, time, pulse, 0);
    else
        conv_signal = pdW; tconv = time;
    end
    
else  % use weights from optimisation (must load file)
    pdW = zeros(Ntraj, Nq, Nts);
    for i=1:Ntraj
        [pdW(i,:,:)] = theory_signal(Q(:,:,i,:), kin, atmnum, q_exp, time, FLAGelec, FLAGinel, FLAGsignal);
    end
    
    if pulse ~= 0
        for i=1:Ntraj
            [temp_signal(i,:,:), tconv, sigma] = convolute(squeeze(pdW(i,:,:)), time, pulse, 0);
        end
        
    else
        temp_signal = pdW; tconv = time;
    end
    
    Nts = length(tconv);
    conv_signal = zeros(Nq, Nts);
    
    for ts=1:Nts
        for tr=1:Ntraj
            temp = weights(tr) * temp_signal(tr, 1:Nq, ts);
            conv_signal(1:Nq, ts) = conv_signal(1:Nq, ts) + temp.';
        end
    end
    
end
  
              
conv_signal = conv_signal * (100 * xfrac);

Ntc = length(tconv);
conv_signal = reshape(conv_signal, [1, Nq, Ntc]);

[binned_signal, Tth_bin, Iexp, TE, CM, T0_exp] = match_signal(tconv, dt, conv_signal, Iexp, Texp, T0, T0_exp, CM, Confidence_Tol, qmin_closest, qmax_closest, 1, Nq, FLAGtdelay);

binned_signal = squeeze(binned_signal);
conv_signal = squeeze(conv_signal);

%figure
%[QQ, TT] = meshgrid(q_exp, tconv);
%mesh(QQ, TT, conv_signal.');
%axis tight
%xlim([0 8]);
%view(0, 90);
%title('All Time')
%ylabel(['Time (fs)'], 'interpreter', 'latex')
%xlabel(['s (\AA', '$ ^{-1} $', ')'], 'interpreter', 'latex')

% figure
% [QQ, TT] = meshgrid(q_exp, TE);
% mesh(QQ, TT, binned_signal.');
% axis tight
% xlim([0 8]);
% view(0, 90);
% title('Exp Time Delays')
% ylabel(['Time (fs)'], 'interpreter', 'latex')
% xlabel(['s (\AA', '$ ^{-1} $', ')'], 'interpreter', 'latex')

%save(fout, '-v7.3')

rsignal{p} = conv_signal;
bsignal{p} = binned_signal;
times{p} = tconv;
TEs{p} = TE;
qs{p} = q_exp;

clear conv_signal; clear binned_signal; clear tconv; clear TE; clear temp_signal;

end
