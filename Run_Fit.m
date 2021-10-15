clear all
close all


%%%%%%%%%%%%%% READ BEFORE RUNNING THIS CODE! %%%%%%%%%%%%%%

% Main wrapper for fitting algorithm. Works for both UED and XRS.
% Currently two target functions - minimisation of constrained nonlinear multivariable function (fmincon),
% and non-linear least squares (lsqnonlin) - lsqnonlin recomended. 
% You can always add your own target function to optimise and/or modify
% algorithms in section 6 of fit_traj_main.m 
% Current target functions are based on sum of squared errors.

% Can include excitation fraction in the optimisation, or scan over a series of ranges.
% In addition, scan a range of experimental and theoretical time shifts.
% One can run a global fit (optimises weights, excitation fraction and T0
% in one shot). Alternativley, set FLAG_T0 = 1 AND choose qlims for
% integration to perform an independent fit of T0 based on the assumption
% of equal weights (this fits the initial rise that results from
% convolution with a gaussian of x fs). 

% EDIT the functions that load the experimental data and trajectories as required for your system.
% Be sure to return the correct variables from these functions:
% Experiment needs - Scattering signal, time bins, q vector, and a confidence matrix. 
% If you have no confidence matrix, just return a matrix of ones (same dimensions as the signal).
% Confidence matrix could be no. of photon hits per pixel, standard deviations from bootstrapping, standard error etc.
% Trajectory loading function must return geometries (IN ANGSTROM!!) and
% spin multiplicities if you would like to group the trajectories into classes of spins and fit each class.
% If you would like to fit individual trajectories - just return a matrix
% of zeros/ ones with length number of trajectories, it does not matter as long as FLAGtfunc = 0.

% BEWARE!!
% It is possible to return some key outputs from the fit_traj_main main function - see documentation for that function.
% You can find documentation for every important function used at the top of the corrosponding .m file (just like this one).
% Compton intensities for inealstic corrections for your atoms may not be paramaterised by default.
% You may need to tabulate the parameters for your atom in compton_intensities.m
% Same goes for the form factors - f_functions.m and f_functions_electron.m


diary OPT_weights_std_20.diary %EDIT WITH EACH SUBSEQUENT RUN


%%%%% SETUP FLAGS %%%%% EDIT BEFORE RUNNING

FLAGpolar = 0; % 1 - include polarisation correction, 0 - do not
FLAGinel = 1; % 1 - include inelastic terms
FLAGelec = 1; % 0 for X-ray, 1 for UED.
FLAGsignal = 0; % 0 = fit with dI/I, 1 = fit using dsM
FLAGtdelay = 0; % 0 = no binning on theory (delays only), 1 = bin theory
FLAGopt = 2; % 0 - fmincon (IP), 1 - fmincon (AS), 2 - lsqnonlin
FLAGtfunc = 0; % 0 - individual trajs, 1 - singlet, triplet, non-diss classes, 2 - non-diss & diss classes
FLAGxfrac = 1; % Include xfrac in optimisation, must set to a scalar guess in input
FLAGconfmat = 1; % Include a confidence matrix that assigns a measure of confidence in the data at each point. 1 = Include, 0 = Exclude.
FLAGexclude = 0; % 1 = exclude certain number trajectories from opt - specified in ex_traj, 0 = exclude none.
Npar = 2; % number of processors to run using. 1 = serial execution 
DEBUG = 0; % 1 for debugging info
FLAG_T0 = 0; % 0 = perform global fit, 1 = perform independent T0 fit
FLAG_wtype = 1; % how initial weights are generated. 0 = N random weights on interval [0, 1]
                % 1 = N random weights on interval [mean_weight-weight_variation, mean_weight+weight_variation]
                % 2 = N weights generated on a harmonic function w some force constant k

%%%%% INPUT PARAMS %%%%% EDIT BEFORE RUNNING

OPT_Tol = [1.d-5 1.d-5 1.d-6]; % Thresh for opt: 1 - Func. Term. Tol, 2 - lower bound on stepsize, 3 - Max dx for finite differances
OPT_Bounds = [0.0, 10.0]; % Bounds for opt
Confidence_Tol = [0]; %, 0.45, 0.50, 0.55, 0.60, 0.65]; % Values of confidence less than this thresh are given 0 weight in the optimisation. Set = 0 to include all.
T0_exp = [80]; % Experimental time zero to be selected/ scanned over
T0_theory = 0; % Theory value of T0 - WARNING: RECOMENDED TO KEEP CONSTANT AT 0 AND SCAN T0_exp
xfrac = [0.03]; % excitation fracs to scan - in percentage points.
pulses = [230]; % pulse FWHM for convolution
atmnum = [6 16 16]; % atomic numbers
kin = 1.8751e+03; % kin ued
dt = 0.5; % time step in theory simulations
Tlen = 1000; % max time in fs to fit over
q_range = [1, 12]; % q range to use for fitting
qlims = [2.8 4.2]; % limits for integration if fitting T0
ex_traj = [14 20 36 93 96 98 100 137 176 197]; % trajectory numbers to exclude from opt
ninit_conds = 5; % number of initial guess conditions for coefficients
weight_std = 20; % std. dev. on average (1/ntraj) weight for sampling in accordance with FLAG_wtype = 1

%%%%% PATHS TO EXPERIMENTAL AND THEORETICAL DATA %%%%% EDIT BEFORE RUNNING

fpath_exp = '/home/kyle/2TB_HDD/OPTDATA/INPUTS/FixedExperimental.mat'; % path to experimental data
fpath_traj = '/home/kyle/2TB_HDD/OPTDATA/INPUTS/Filtered_Trajs_Final.mat'; % path to theory data
fname = 'OPT_Weights_std_20'; % Prefix of .mat output file data will be saved to.

% Edit these two functions to load relevent data
[Texp, Iexp, q_exp, CM] = load_experiment(fpath_exp, FLAGconfmat); % need experimental time vec, signal and q range
[Q, multiplicity] = load_trajectories(fpath_traj); % need geometries and spin of trajectories


%%%%%%%%%%%%%% END OF INPUT SPECIFICATION! %%%%%%%%%%%%%%

time_start = tic;

% Some basic checks prior to running
if FLAG_T0==1; disp('PERFORMING OPTIMISATION OF T0'); end
if FLAG_T0==0; disp('PERFORMING GLOBAL OPTIMISATION'); end
if FLAGelec==0; disp('WILL CALCULATE XRAY SCATTERING'); end    
if FLAGelec==1; disp('WILL CALCULATE ELECTRON SCATTERING'); end 
if FLAGsignal==0; disp('ALL CALCULATED SIGNALS ARE dI/I'); end
if FLAGsignal==1; disp('ALL CALCULATED SIGNALS ARE dsM'); end
if FLAGxfrac==1; disp('EXCITATION FRACTION INCLUDED IN OPT'); end
if Npar>1; disp(['Parallel calculation, Npar=' num2str(Npar)]); end
if DEBUG==1; disp('RUNNING IN DEBUG MODE'); end
if FLAGxfrac==1 & isscalar(xfrac) ~= 1; error('XFRAC MUST BE A SCALAR GUESS IF INCLUDING IN OPT'); end
if FLAGconfmat == 0; Confidence_Tol = 0; end
if FLAGexclude == 1; disp(['Excluding Trajectories...', num2str(ex_traj)]); end


% Loop over each time shift/ xfrac and execute optimisation algorithm
for i=1:length(pulses)
for j=1:length(xfrac)
   for k=1:length(T0_exp)
           if Npar>1 delete(gcp('nocreate')); end
           fout = [fname, '_PULSE_', num2str(pulses(i)), '_xfrac_', num2str(xfrac(j)), '_T0E_', num2str(T0_exp(k)), '.mat']; 
           disp(['-------------PULSE SCAN ITER: ', num2str(i), '-----------']);
           disp(['-------- EXFRAC SCAN ITERATION: ', num2str(j), ' --------']);
           disp(['-------- TIME SCAN ITERATION: ', num2str(k), ' --------']);
           pulse = pulses(i);
           exfrac = xfrac(j);
           Tshift_exp = T0_exp(k);
         
           fit_traj_main3(exfrac, Tshift_exp, T0_theory, Tlen, q_range, Texp, dt, Iexp, q_exp,...
                         Q, multiplicity, pulse, atmnum, kin, fout, FLAGpolar, FLAGinel, FLAGelec,...
                         FLAGopt, FLAGtfunc, Npar, OPT_Tol, OPT_Bounds, DEBUG, FLAGxfrac, CM, Confidence_Tol,...
                         FLAGexclude, ex_traj, FLAGsignal, ninit_conds, FLAGtdelay, qlims, FLAG_T0, FLAG_wtype, weight_std);
       end
   end
end


if Npar>1; delete(gcp('nocreate')); end
telapsed_total = toc(time_start);
disp(['Time elapsed for TOTAL calculation (s):' num2str(telapsed_total)]);
disp(['Time elapsed for TOTAL calculation (min):' num2str(telapsed_total/60)]);
disp(['Time elapsed for TOTAL calculation (hrs):' num2str(telapsed_total/3600)]);
disp(['-------------------- FINISHED ALL CALCULATIONS! --------------------']);

diary off;
