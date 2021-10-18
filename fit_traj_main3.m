function [weight_final, weight_init, exfrac_final, Ff, Fi] = fit_traj_main3(exfrac, T0_exp, T0, Tlen, q_range, Texp, dt, Iexp, q_exp, Q, multiplicity, pulse, atmnum, kin, fout, FLAGpolar, FLAGinel, FLAGelec, FLAGopt, FLAGtfunc, Npar, OPT_Tol, OPT_Bounds, DEBUG, FLAGxfrac, CM, Confidence_Tol, FLAGexclude, ex_trajs, FLAGsignal, ninit_conds, FLAGtdelay, qlims, FLAG_T0, FLAG_wtype, weight_ub, prev_weights)

% INPUTS:
% exfrac - excitation fraction in percentage units - either a guess to be optimised or an explicit weight
% T0_exp - experimental time zero - will select bin value closest to the one requested - in fs
% T0 - theory time zero shift relative to experimental - i.e. T0_exp + T0 in fs
% Tlen - time range in fs to fit over
% q_range - minimum q value and maximum q value to fit over in inverse Angstroms - i.e. q_range = [0 6.2];
% Texp - experimental time vector (binned)
% dt - time step in theory calculations - must match time step in trajectories
% Iexp - experimental signal - binned. Must be dimensions (q, time)
% q_exp - experimental momentum transfer vector in inverse angstrom
% Q - Nuclear coordinates from trajectories - in angstrom
% multiplicity - spin multiplicity of trajectories - length(Ntraj) - 0 = non-diss, 1 = singlet, 2 = triplet
% pulse - fwhm of pump/probe pulses in fs for convolution - i.e. pulse = [150 0]. Will be summed
% atmnum - atomic numbers of each nuclei
% kin - incident wave vector
% fout - name of file to save output data to
% FLAGpolar - FLAG for including polarisation correction, 0 = none, 1 = corrected
% FLAGinel - include constant inelastic term calculated from comptom intensities, 1 = include, 0 = do not
% FLAGelec - 0 for calculating Xray Scattering, 1 for calculating Electron Diffraction signals
% FLAGopt - type of optimisation : 0 - fmincon (IP), 1 - fmincon (AS), 2 - lsqnonlin 
% FLAGtfunc - tfunc to fit: 0 - individual trajs, 1 - singlet, triplet, non-diss classes, 2 - non-diss & diss classes 
% Npar - Number of cores to run in parallel
% OPT_tol - thresh for opt: 1 - Func. Term. Tol, 2 - lower bound on stepsize, 3 - Max dx for finite differances
% OPT_Bounds - bounds for opt e.g. [1, 10] in form [lower_bound, upper_bound]
% FLAGxfrac - include excitation fraction as an additional factor in the opt - using exfrac as a init guess
% BinSize - in event of non-linear bins, specify the bin spacing you would like to spline experiment onto
% CM  - confidence matrix giving each q and time value a weight in opt
% Confidence_Tol - Confidence value below which data points not included
% FLAGexclude - Flag to exclude certain trajectories in optimisation - trajs given in ex_trajs
% ex_trajs - Numbers of trajectories to be excluded from optimisation


% OUTPUTS (IF WANTED TO RETURN TO ANOTHER WRAPPER):
% binned_signal - percentage difference signal - binned to match experimental bins
% conv_signal - pre-binned but convoluted (if done) signal
% weight_final - final weights after opt, weight(end) is the optimised exfrac if included
% Ff - values of the target function after optimisation
% TE - experimental time bins over range used in opt
% Ttheory - theory time vector over range used in opt
% qAng - momentum transfer vector over range used in opt

tstart = tic;

% 1 -  INIT 

global au2ang ang2au
au2ang = 0.52917721092d0;
ang2au = 1/au2ang;

[natm, ~, Ntraj, Nts] = size(Q);

if natm ~= length(atmnum)
    error('NUMBER OF ATOMS IN TRAJS AND NUMBER OF ATOMIC NUMBERS INCONSITENT.')
end


%% 2 - Calculate IAM for each trajectory

disp(['>>>> CALCULATING THEORY PATTERN OVER SELECTED TIME FRAME. <<<<']);

time = 0:dt:floor(Nts*dt); % Theoretical time vector for calculating IAM signal and convolving


qmin(1:length(q_exp)) = q_range(1); 
qmax(1:length(q_exp)) = q_range(2);

[~, qmin_closest] = min(abs(qmin - q_exp)); % closest q values in experiment vs requested q range
[~, qmax_closest] = min(abs(qmax - q_exp));
qAng = q_exp(qmin_closest:qmax_closest); % selected q range for experiment (inv Angstrom)
Nq = length(qAng);
    

if length(time) ~= Nts; error("Theory time vector does not match trajectory dimensions"); end

if Npar>1
    parpool(Npar);
end

exfrac_init = 1;

switch FLAGtfunc
    case 0  % individual trajs
        pdW = zeros(Ntraj, Nq, Nts);  % percentage difference
        if Ntraj < 2 && Npar > 1
            warning('DO NOT RUN IN PAR FOR ONE TRAJ.')
            Npar = 1;
        end

        if Npar > 1 % parallel case
            parfor ntr=1:Ntraj
                Qtr = Q(:,:,ntr,:);
                [pdW(ntr,:,:)] = theory_signal(Qtr, kin, atmnum, qAng, time, FLAGelec, FLAGinel, FLAGsignal);
            end
        else
            for ntr=1:Ntraj % run in serial
                Qtr = Q(:,:,ntr,:);
                [pdW(ntr,:,:)] = theory_signal(Qtr, kin, atmnum, qAng, time, FLAGelec, FLAGinel, FLAGsignal);
            end
        end
        
        nclass = Ntraj; % number of classes in opt same as number of trajs
        
    case {1, 2} % bound, singlet, triplet || bound, diss
        %keyboard
        Q_classes = sort_traj_classes(Q, multiplicity, FLAGtfunc);
        nclass = length(Q_classes);
        pdW = zeros(nclass, Nq, Nts);
       
       
        
        if Npar > 1
            parfor c=1:nclass
                Qtr = Q_classes{c};
                [pdW(c,:,:)] = theory_signal(Qtr, kin, atmnum, qAng, time, FLAGelec, FLAGinel, FLAGsignal);
            end
        else
            for c=1:nclass
                Qtr = Q_classes{c};
                [pdW(c,:,:)] = theory_signal(Qtr, kin, atmnum, qAng, time, FLAGelec, FLAGinel, FLAGsignal);
            end
        end
end


if DEBUG == 1 % inspect raw data set
  
    [QQ, TT] = meshgrid(qAng, time);
    [pdW_avg] = theory_signal(Q, kin, atmnum, qAng, time, FLAGelec, FLAGinel);
    
    figure
    mesh(QQ,TT,pdW_avg.')
    title('Average Theoretical Signal')
    axis tight
    view(0,90)
    saveas(gcf, 'Avg_Theory_Signal_DEBUG.pdf')
    
end

 
%% 3 - Polarisation Correction

switch FLAGpolar
    case 0 
        disp('FLAGpolar=0, do not account for polarisation')
        P = ones(length(qAng),1); % polarisation factor
    case 1
        if FLAGelec == 1 error('NO POLARISATION FOR ELECTRON BEAM'); end
        disp('FLAGpolar=1, account for polarisation (rotational average)')
        theta = asin(qAng/max(qAng));
        P = 0.5d0 * (1.d0 + cos(theta).^2).'; % polarization factor - rot.avg.
end

for c=1:nclass
    for ts=1:Nts
        pdW(c,1:Nq,ts) = P(1:Nq) .* pdW(c,1:Nq,ts)';
    end
end


if DEBUG == 1 && FLAGpolar == 1 % inspect polarised data set
    figure
    plot(qAng, P, 'r')
    title('Inspect Polarisation')
    xlabel('q (A^{-1})')
    ylabel('Polarisation Factor (P)')
    saveas(gcf, 'Polarisation_Factor_DEBUG.pdf')
end
       
%% 4 - Convolute Theory Data Set

disp(['>>>> PERFORMING CONVOLUTION ON THEORY DATASET. <<<<']);

conv_signl = zeros(nclass, Nq, Nts);

if Npar > 1
    
    tconv = cell(1,nclass);
    sigma = cell(1,nclass);
    
    parfor ntr=1:nclass
        [conv_signal(ntr,:,:), tconv{ntr}, sigma{ntr}] = convolute(squeeze(pdW(ntr,:,:)),time,pulse,0);
    end
    
    tconv = tconv{1};
    sigma = sigma{1};
    
else
    for ntr=1:nclass
        [conv_signal(ntr,:,:), tconv, sigma] = convolute(squeeze(pdW(ntr,:,:)),time,pulse,0);
    end
end

Nts = length(tconv);

if DEBUG == 1 % inspect convoluted data set
    [QQ, TT] = meshgrid(qAng, tconv);
    
    avg_conv = zeros(Nq, length(tconv));
    for i=1:Ntraj
        avg_conv = avg_conv + squeeze(conv_signal(i,:,:));
    end
    avg_conv = avg_conv./Ntraj;

    figure
    mesh(QQ, TT, avg_conv.')
    title('Average Convoluted Difference Signal')
    axis tight
    view(0,90)
    saveas(gcf, 'Avg_Convolved_Signal_DEBUG.pdf')
end

%% 4 - Match Theory and Experimental in Time

disp(['>>>> MATCHING THEORY AND EXPERIMENTAL DATA GRIDS. <<<<']);

[binned_signal, Tth_bin, Iexp, TE, CM, T0_exp] = match_signal(tconv, dt, conv_signal, Iexp, Texp, T0, T0_exp, CM, Confidence_Tol, qmin_closest, qmax_closest, nclass, Nq, FLAGtdelay);

if 0>1
Tth_bin = Tth_bin(4:end);  % HACK WARNING
binned_signal = binned_signal(:,:,4:end);
TE = TE(4:end);
Iexp = Iexp(:, 4:end);
CM = CM(:, 4:end);
end

Nts = length(Tth_bin);

if DEBUG==1
    [QQ, TT] = meshgrid(qAng, Tth_bin);
    figure
    mesh(QQ, TT, (squeeze(sum(binned_signal, 1))./Ntraj).')
    title('Equally Weighted (Pre-Optimised) Theoretical Signal')
    view(0,90)
    axis tight
    saveas(gcf, 'Avg_PreOpt_Matched_Signal_DEBUG.pdf')
end


%% 5 - Integrating over q limits - for fitting T0 rise independantly

if FLAG_T0 == 1
    
    disp(['>>>> Integrating Signal Over q Limits (T0 Fitting). <<<<']);

    Ith = squeeze(sum(binned_signal, 1))./Ntraj;

    qlb(1:length(qAng)) = qlims(1);
    qub(1:length(qAng)) = qlims(2);

    disp(['Requested q integration range: ', num2str(qlims(1)), ' : ', num2str(qlims(2))])

    [~, qlb_closest] = min(abs(qlb - qAng)); % closest T0 value in experiment vs requested T
    [~, qub_closest] = min(abs(qub - qAng)); % closest Tf value in experiment vs requested T

    qlb = qAng(qlb_closest);
    qub = qAng(qub_closest);

    disp(['Selected q integration range: ', num2str(qlb), ' : ', num2str(qub)])

    qlb_ind = find(qAng == qlb); qub_ind = find(qAng == qub);
    
    dq = qAng(2)-qAng(1);

    Ith_Integrated = zeros(1,Nts);
    Iexp_Integrated = zeros(1,Nts);

    for ts=1:Nts
        for qq=qlb_ind:qub_ind
            Ith_Integrated(ts) = Ith_Integrated(ts) + Ith(qq, ts);
            Iexp_Integrated(ts) = Iexp_Integrated(ts) + Iexp(qq, ts);
        end
    end

    Ith_Integrated = Ith_Integrated .* dq;
    Iexp_Integrated = Iexp_Integrated .* dq;
    
    if DEBUG == 1
        figure
        plot(Tth_bin, Iexp_Integrated, '-b')
        hold on
        plot(Tth_bin, Ith_Integrated*(100*exfrac), '-r')
        axis tight
        legend('Iexp Integrated', 'Ith Integrated * xfrac (guess)')
        title('Integrated Intensity for Individual T0 Fit. Ith Multiplied by initial guess xfrac.')
        saveas(gcf, 'T0_Integrated_Intensity_DEBUG.pdf')
    end
    
end



%% 6 - Fit each trajectory to data and store in weight vector

disp(['>>>> FITTING TO EXPERIMENTAL DATA. <<<<']);

if FLAG_T0 == 0 % GLOBAL FITTING

    weight_init = [];
    
    switch FLAG_wtype
        case 0

            for i=1:ninit_conds
                w_add(i, 1:nclass) = rand(1, nclass);
                w_add(i, 1:nclass) = w_add(i, 1:nclass)/sum(w_add(i, 1:nclass)); % normalize to unity;
                w_add(i, nclass+1) = exfrac; % include exfrac as additional term
            end
            
        case 1 % generate a distribution of weights around the mean with some std. dev.
            
            mean_weight = 1/nclass; % averaged weights
            %weight_lb = mean_weight - (mean_weight*weight_std);

            %if weight_lb < 0
            %    weight_lb = 0;
            %end
            
            weight_lb = 0;
            
            %weight_ub = mean_weight + (mean_weight*weight_std);
            
            OPT_Bounds(1) = weight_lb; % change bounds according to
            OPT_Bounds(2) = weight_ub; % calculated bounds using std. dev. given
            
            [w_add, ~] = randfixedsum(nclass, ninit_conds, 1, weight_lb, weight_ub);

            
            [~, npw] = size(prev_weights);

            if npw > 0
                w_add = [w_add  prev_weights];
                ninit_conds = ninit_conds + npw;
            end
            
            
        case 2
            
    end

    weight_init = [];

    weight_init = [weight_init; w_add]; % initial guess weights for opt
    

    OPT_Verbose = 1;  % Verbose printing of output for each step. Change to zero to turn off.

    switch FLAGopt
        case 0 
            opt = optimset('fmincon');
            opt = optimset(opt, 'Algorithm','interior-point');
            Aeq(1,1:nclass) = 1; % one of each weight - should sum to 1
            beq = 1; 
            lb(1,1:nclass) = OPT_Bounds(1); % remove bounds and ensure normarlised within tfunc
            ub(1,1:nclass) = OPT_Bounds(2); 

            if FLAGxfrac == 1 % add constraints for optimised xfrac
                Aeq(nclass+1) = 0; % ex.frac not included in norm
                lb(nclass+1) = 0;
                ub(nclass+1) = 100; % max percentage for ex.frac
            end

            tfunc = 'fmincon_tfunc';
        case 1
            opt = optimset('fmincon');
            opt = optimset(opt,'Algorithm','active-set');
            Aeq (1,1:nclass) = 1; % one of each weight - sum to 1
            beq = 1; 
            lb(1,1:nclass) = OPT_Bounds(1); % remove bounds and ensure normarlised within tfunc
            ub(1,1:nclass) = OPT_Bounds(2); 

                Aeq(nclass+1) = 0;
                lb(nclass+1) = 0;
                ub(nclass+1) = 1; % max exfrac

        case 2
             opt = optimset('lsqnonlin');
             tfunc = 'lsq_tfunc';
             Aeq (1,1:nclass) = 1; % one of each weight
             beq = 1;  % sum to one
             lb(1:nclass) = OPT_Bounds(1); % no norm in tfunc 
             ub(1:nclass) = OPT_Bounds(2); % change to mean +/- std. dev.
             Aeq(nclass+1) = 0; 
             lb(nclass+1) = 0;
             ub(nclass+1) = 1; % max exfrac
    
    end

    opt = optimset(opt,'TolFun',OPT_Tol(1),'TolX',OPT_Tol(2),'DiffMaxChange',OPT_Tol(3),'FunValCheck','on');

    if OPT_Verbose == 1   % output info for each step     
        opt = optimset(opt,'Display','iter','Diagnostic','on'); 
    end

    weight_final = zeros(size(weight_init));

    if OPT_Bounds(1) ~= 0 warning('LOWER BOUND OF OPT NON-ZERO!'); end
    if OPT_Bounds(2) ~= 1 warning('UPPER BOUND OF OPT NOT UNITY! IGNORE IF TARGET FUNC INCLUDES NORMALISATION'); end
    opt  % print optimisation configuration

    % lb = [];
    % ub = [];
    if FLAG_wtype == 0
        Aeq = []; % HACK: if normalising within target function - Aeq and beq are not needed.
        beq = []; 
    end

    if isempty(Aeq) || isempty(beq) warning('NO LINEAR EQUALITY CONSTRAINTS SET'); end
    if isempty(lb) || isempty(ub) warning('NO BOUNDS ON OPT SET - OK IF TFUNC INCLUDES NORMALISATION'); end
    
    for i=1:ninit_conds
        disp(['Optimize ' num2str(i) '/' num2str(ninit_conds)]);
        switch FLAGopt
            case {0, 1} % fmincon
                Fi(i) = feval(tfunc, weight_init(i,:), binned_signal, Iexp, length(TE), nclass, Nq, FLAGxfrac)./Nts;
                [weight_final(i,:), z1, z2, flag] = fmincon(tfunc, weight_init(i,:), [], [], Aeq, beq, lb, ub, [], opt, binned_signal, Iexp, Nts, nclass, Nq, FLAGxfrac);
                weight_final(i, 1:nclass) = weight_final(i, 1:nclass) / sum(weight_final(i, 1:nclass)); % normalise
                exfrac_final(i) = weight_final(i, end); % seperate xfrac if optimising
                Fi(i) = feval(tfunc, weight_init(i,:), binned_signal, Iexp, length(TE), nclass, Nq, FLAGxfrac)./Nts;
                Ff(i) = feval(tfunc, weight_final(i,:), binned_signal, Iexp, length(TE), nclass, Nq, FLAGxfrac)./Nts;
            case 2 % lsq
               
                [weight_final(i,:), z1, z2, flag] = lsqnonlin(tfunc, weight_init(i,:), lb, ub, opt, binned_signal, Iexp, Nts, nclass, Nq, CM, FLAGxfrac, FLAGexclude, ex_trajs, FLAG_wtype);
                %weight_final(i, 1:nclass) = weight_final(i, 1:nclass) / sum(weight_final(i, 1:nclass));
                exfrac_final(i) = weight_final(i, end);
                if flag < 0 warning('Optimisation Failed'); end;
                fi = feval(tfunc, weight_init(i,:), binned_signal, Iexp, length(TE), nclass, Nq, CM, FLAGxfrac, FLAGexclude, ex_trajs, FLAG_wtype);
                ff = feval(tfunc, weight_final(i,:), binned_signal, Iexp, length(TE), nclass, Nq, CM, FLAGxfrac, FLAGexclude, ex_trajs, FLAG_wtype);
                Fi(i) = sum(fi(1:numel(fi)).^2);
                Ff(i) = sum(ff(1:numel(ff)).^2);
                clear fi ff
        end

        disp(['% INIT FUNC VAL: ', num2str(Fi(i))])
        disp(['% OPT FUNC VAL: ', num2str(Ff(i))])
    
    end


else   % FITTING OF T0 ONLY
    
    opt = optimset('lsqnonlin');
    tfunc = 'lsq_tzero';
    Aeq = [];
    beq = [];
    lb(1) = 0;
    ub(1) = 1; % max exfrac

    opt = optimset(opt,'TolFun',OPT_Tol(1),'TolX',OPT_Tol(2),'DiffMaxChange',OPT_Tol(3),'FunValCheck','on');

    OPT_Verbose = 1;

    if OPT_Verbose == 1   % output info for each step     
        opt = optimset(opt,'Display','iter','Diagnostic','on'); 
    end

        [exfrac_final, z1, z2, flag] = lsqnonlin(tfunc, exfrac , lb, ub, opt, Ith_Integrated, Iexp_Integrated);
        if flag < 0 warning('Optimisation Failed'); end;
        fi = feval(tfunc, exfrac, Ith_Integrated, Iexp_Integrated);
        ff = feval(tfunc, exfrac_final, Ith_Integrated, Iexp_Integrated);
        Fi = sum(fi(1:numel(fi)).^2)./Nts;
        Ff = sum(ff(1:numel(ff)).^2)./Nts;
        clear fi ff 

        disp(['% INIT FUNC VAL: ', num2str(Fi) ])
        disp(['% OPT FUNC VAL: ', num2str(Ff) ])
    
end


%% 7 - Save fitting parameters 

weight_final = weight_final(:, 1:nclass);
weight_init = weight_init(:, 1:nclass);
telapsed = toc(tstart);
disp(['Time elapsed for ITER   (s):' num2str(telapsed)]);
disp(['Time elapsed for ITER (min):' num2str(telapsed/60)]);
disp(['Time elapsed for ITER (hrs):' num2str(telapsed/3600)]);

save(fout, 'weight_final', 'weight_init', 'exfrac_final', 'Fi', 'Ff', '-v7.3')
disp(['SAVING DATA TO :', fout]);
disp(['-------------------- FINISHED SCAN ITER! --------------------']);


end

