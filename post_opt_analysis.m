clear all
close all

% NEEDS REFACTORING - WRITTEN VERY QUICKLY 

fprefix = 'OPTall_inel_q1_5_lsqn_xfrac';
xfrac = [2.9];
init_T = [20];
FLAGplot = 1;

%diary lsqnonlin_2_Analysis.diary
%diary TEST_ANAL2.diary

iter = 0;
for a=1:length(xfrac)
    for b=1:length(init_T)
    iter = iter + 1;
    
    %fpath = [fprefix, '_', num2str(xfrac(a)), '_T0_', num2str(init_T(b)), '.mat'];
    fpath = 'OPT_Weights_std_0.7_PULSE_230_xfrac_0.03_T0E_80.mat'
    disp(['>>> LOADING DATA FROM: ', fpath])
    load(fpath)
    
    disp(['% REQUESTED EXPERIMENTAL Q RANGE: ', num2str(q_range(1)), ' TO ', num2str(q_range(2))]);
    disp(['% SELECTED EXPERIMENTAL Q RANGE: ', num2str(qAng(1)), ' TO ', num2str(qAng(end))]);
    disp(['% REQUESTED EXPERIMENTAL TIME RANGE: ', num2str(T0), ' TO ', num2str(T0+1000)]);
    disp(['% SELECTED EXPERIMENTAL TIME RANGE: ', num2str(TE(1)), ' TO ', num2str(TE(end))]);

    disp('>>>> Analysis of Optimisation <<<<')

    if FLAGxfrac == 1
        for i=1:ninit_conds
            if weight_final(i,end) > 100 || weight_final(i,end) < 0.0
                warning('Xfrac value not within bounds.')
            end
        end
    end

    disp(['INIT FUNC VALUES: ', num2str(Fi(1)), ', ', num2str(Fi(2)), ', ', num2str(Fi(3))])
    disp(['FINAL FUNC VALUES: ', num2str(Ff(1)), ', ', num2str(Ff(2)), ', ', num2str(Ff(3))])        
    best_opt = find(Ff == min(Ff));
    disp(['BEST OPT - NUM. : ', num2str(best_opt), ' FURTHER ANALYSIS BASED ON THIS OPT']);
    
    if sum(weight_final(best_opt,1:nclass)) ~= 1.0 warning('Final Weights NOT Normalised!.'); end

    c = 0;
    for i=1:nclass
        if weight_final(best_opt,i) > mean_weight
            c = c + 1;
            important_trajs(c) = i;
            traj_spin(c) = multiplicity(i);
            traj_weights(c) = weight_final(best_opt, i);
        end
    end

    disp(['>>> TRAJS WITH > 1% WEIGHTING: '])
    for i=1:c
        disp(['NO.', num2str(i), ' TRAJ: ', num2str(important_trajs(i)), ', SPIN: ', num2str(traj_spin(i)), ', WEIGHT: ', num2str(traj_weights(i))])
    end
    disp(['IN TOTAL: ', num2str(c), ' TRAJS HAVE > 1% WEIGHTING. THESE MAKE UP ', num2str(sum(traj_weights)), ' OF THE TOTAL WEIGHT']) 

 
    if FLAGxfrac == 1 
        opt_xfrac = weight_final(best_opt, end);
        disp(['>>> EXFRAC OPT :', num2str(opt_xfrac)]);
    else
        opt_xfrac = xfrac;
    end
    
    singlet = 0;
    triplet = 0;
    bound = 0;
    for i=1:c
        if traj_spin(i) == 1
            singlet = singlet + traj_weights(i);
        elseif traj_spin(i) == 3 
            triplet = triplet + traj_weights(i);
        elseif traj_spin(i) == 0
            bound = bound + traj_weights(i);
        end
    end
    
    disp(['>>> FROM TRAJS WITH WEIGHTS > 1% : '])
    disp(['Singlet Weight: ', num2str(singlet*100)])
    disp(['Triplet Weight: ', num2str(triplet*100)])
    disp(['Branching Ratio: ', num2str(triplet/singlet)])
    disp(['Bound Weight: ', num2str(bound*100)])
    disp(['Total: ', num2str((singlet+triplet+bound)*100)])
    
    singlet = 0;
    triplet = 0;
    bound = 0;
    for i=1:nclass
        if multiplicity(i) == 1
            singlet = singlet + weight_final(best_opt, i);
        elseif multiplicity(i) == 3 
            triplet = triplet + weight_final(best_opt, i);
        elseif multiplicity(i) == 0
            bound = bound + weight_final(best_opt, i);
        end
    end
    
    disp(['>>> FROM ALL TRAJS : '])
    disp(['Singlet Weight: ', num2str(singlet*100)])
    disp(['Triplet Weight: ', num2str(triplet*100)])
    disp(['Branching Ratio: ', num2str(triplet/singlet)])
    disp(['Bound Weight: ', num2str(bound*100)])
    disp(['Total: ', num2str((singlet+triplet+bound)*100)])
    
    
    theory_fit = zeros(Nq,length(TE));
    for i=1:nclass
        theory_fit = theory_fit + (squeeze(pdW_bin(i,:,:) * weight_final(best_opt, i)));
    end
    theory_fit = theory_fit * opt_xfrac; % weighted (optimised) average of all trajectories
    
    theory_prefit = ((squeeze(sum(pdW_bin, 1)))./Ntraj)*opt_xfrac; % pre-optimised (equal weight) trajs
        
    
    theory_gt1_fit = zeros(Nq, length(TE)); % only important optimised (> 1% weight) trajs
    theory_gt1_prefit = zeros(Nq, length(TE)); % only important pre-optimised (> 1% weight) trajs
    c = 0;
    for i=1:nclass
        if ismember(i, important_trajs)
            c = c + 1;
            Qmain(:,:,c,:) = Q(:,:,i,:);
            multiplicity_main(c) = multiplicity(i);
            theory_gt1_fit = theory_gt1_fit + (squeeze(pdW_bin(i,:,:) * traj_weights(c)));
            theory_gt1_prefit = theory_gt1_prefit + squeeze(pdW_bin(i,:,:));
        end
    end
    theory_gt1_fit = theory_gt1_fit * opt_xfrac;
    theory_gt1_prefit = (theory_gt1_prefit./length(important_trajs)) * opt_xfrac;
   
    
    DT = [0 dt_exp];
    inds = [1];
    t = 0;
    for i=1:length(DT)
        t = t + DT(i);
        if t > 180 & t < 260 % HACK - HARDCODED RANGES
            inds = [inds i]; % time indexes to bin over in lineouts - 200 fs spacing
            t = 0;
        end
    end
        
    avg_dt = mean(diff(TE(inds))); % average bin length - should be ~ 200 fs
    
    lineout_exp = zeros(Nq, length(inds)-1);
    lineout_fit = zeros(Nq, length(inds)-1);
    lineout_prefit = zeros(Nq, length(inds)-1);
    
    for j=1:length(inds)-1
        i = j + 1;
        [~, norm_exp] = size(Iexp(:, inds(j):inds(i)));
        [~, norm_fit] = size(theory_fit(:, inds(j):inds(i)));
        [~, norm_prefit] = size(theory_fit(:, inds(j):inds(i)));
        lineout_exp(:,j) = (sum(Iexp(:, inds(j):inds(i)), 2))./norm_exp;
        lineout_fit(:,j) = (sum(theory_fit(:, inds(j):inds(i)), 2))./norm_fit;
        lineout_prefit(:,j) = (sum(theory_prefit(:, inds(j):inds(i)), 2))./norm_prefit;
    end   
 
     % PLOTTING
            
     if FLAGplot == 1
         
        [QQ, TT] = meshgrid(qAng, TE);
        minq = min(qAng);
        maxq = max(qAng);
        
        f = figure;
        mesh(QQ,TT,Iexp') % EXPERIMENTAL SIGNAL
        pt = ['Experimental Signal - T0 =  ', num2str(round(min(TE))), ' fs'];
        ylabel(['Time (fs)'], 'interpreter', 'latex')
        xlabel(['s (\AA', '$ ^{-1} $', ')'], 'interpreter', 'latex')
        ylim([min(TE), max(TE)])
        xlim([minq maxq])
        view(0,90)
        pt = [pt '.fig']
        savefig(f, pt)
        f = figure;
        mesh(QQ,TT,theory_prefit') %THEORY PREFIT
        pt = ['Averaged Pre-Optimised Signal - T0 =  ',  num2str(round(min(TE))), ' fs'];
        ylabel(['Time (fs)'], 'interpreter', 'latex')
        xlabel(['s (\AA', '$ ^{-1} $', ')'], 'interpreter', 'latex')
        view(0,90)
        ylim([min(TE), max(TE)])
        xlim([minq maxq])
        pt = [pt '.fig'];
        savefig(f, pt)
        f = figure;
        mesh(QQ,TT,theory_fit') % THEORY FITTED
        pt = ['Averaged Optimised Signal - T0 = ',  num2str(round(min(TE))), ' fs'];
        ylabel(['Time (fs)'], 'interpreter', 'latex')
        xlabel(['s (\AA', '$ ^{-1} $', ')'], 'interpreter', 'latex')
        view(0,90)
        ylim([min(TE), max(TE)])
        xlim([minq maxq])
        pt = [pt '.fig'];
        %savefig(f, pt)
        f = figure;
        mesh(QQ,TT,(abs(theory_fit-Iexp))') % RESIDUAL
        pt = ['Residual - T0 = ',  num2str(round(min(TE))), ' fs'];
        ylabel(['Time (fs)'], 'interpreter', 'latex')
        xlabel(['s (\AA', '$ ^{-1} $', ')'], 'interpreter', 'latex')
        view(0,90)
        ylim([min(TE), max(TE)])
        xlim([minq maxq])
        pt = [pt '.fig'];
        %savefig(f, pt)
        
        
        

        f = figure;
        time_labels = 0:ceil(avg_dt):1000; % TIME AVERAGED LINEOUTS
        pt = [num2str(ceil(avg_dt)), ' fs Averaged Line Outs - T0 = ', num2str(round(min(TE))), ' fs'];
        for i=1:length(inds)-1
            j = i + 1;
            subplot(length(inds)-1,1,i)
            plot(qAng, lineout_fit(:,i), '-r')
            hold on
            plot(qAng, lineout_prefit(:,i), '-b')
            plot(qAng, lineout_exp(:,i), '--k')
            title([num2str(time_labels(i)), ' - ', num2str(time_labels(j)), ' fs' ], 'FontSize', 8)
            xlim([minq maxq])
            ylim([-0.5 0.5])
            if i == 3
                ylabel(['$\Delta$sM'], 'interpreter', 'latex')
            end
            if i == length(inds)-1
                xlabel(['s (\AA', '$ ^{-1} $', ')'], 'interpreter', 'latex')
                legend('Theory (Optimised)', 'Theory (Averaged)', 'Experiment', 'Orientation', 'horizontal', 'Location', 'south')
            end
            hold on
        end
        pt = [pt '.fig'];
        %savefig(f, pt)
        
        f = figure;
        time_labels = 0:ceil(avg_dt):1000; % LINEOUTS
        pt = ['Line Outs - T0 = ', num2str(round(min(TE))), ' fs'];
        for i=1:length(inds)
            subplot(length(inds),1,i)
            plot(qAng, theory_fit(:,i), '-r')
            hold on
            plot(qAng, theory_prefit(:,i), '-b')
            plot(qAng, Iexp(:,i), '--k')
            title([num2str(time_labels(i)), ' fs' ], 'FontSize', 8)
            xlim([minq maxq])
            ylim([-0.5 0.5])
            if i == 4
                ylabel(['$\Delta$sM'], 'interpreter', 'latex')
            end
            if i == length(inds)
                xlabel(['s (\AA', '$ ^{-1} $', ')'], 'interpreter', 'latex')
                legend('Theory (Optimised)', 'Theory (Averaged)', 'Experiment', 'Orientation', 'horizontal', 'Location', 'south')
            end
            hold on
        end
        pt = [pt '.fig'];
        %savefig(f, pt)
        
        distances(Q, Ttheory, dt, [150 0], squeeze(weight_final(best_opt,:)) , nclass, multiplicity, 2)
        
        keyboard
        distances(Qmain, Ttheory, dt, [150 0], traj_weights , length(traj_weights), multiplicity_main, 1)
        
     end

    disp(['-------------- END OF ANALYSIS FOR FILE ', num2str(iter), ' --------------'])
    
    end
end

%diary off
