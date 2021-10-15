clear all
close all

pulses = [230];
T0s = [80];
fpath = '/home/kyle/2TB_HDD/OPTDATA/INPUTS/FixedExperimental.mat';
%fpath = '/media/2TB_HDD/OPTDATA/LOCAL_VER/EXTEND_T0/RESULTS/MANY_NEGATIVE_TIME_DELAY/FINAL_FIXED/OldExperimental.mat';

Np = length(pulses);
Nt0 = length(T0s);

xfracs = zeros(Nt0, Np);
errf = zeros(Nt0, Np);

for rr=1:Nt0
    for cc=1:Np
        pulse = pulses(cc);
        pulse = sum(pulse);
        T0 = T0s(rr);
        file = ['T0FIT_FIXED_PULSE_', num2str(pulse), '_xfrac_0.03_T0E_', num2str(T0), '_T0T_0.mat']
        load(file); %, 'exfrac_final', 'Ff', 'T0_exp', 'Ith_Integrated', 'binned_signal', 'qAng', 'Ntraj', 'Nts', 'qlb', 'qub', 'Tth_bin', 'Integrated_Iexp');
        
        xfracs(rr,cc) = exfrac_final;
        errf(rr,cc) = Ff;
        dq = qAng(2)-qAng(1);
        
        [~, Iexp_all, q_exp_all, ~] = load_experiment(fpath, 0);
                
        
        Tind = find(TE == T0_exp);
        Tind_all = find(Texp == T0_exp);
        ind = find(Texp == TE(1));
        ind_end = ind + length(TE)-1;
        Texp_all = Texp(1:ind_end);
        
        Iexp_all = Iexp_all(:, 1:ind_end);
        
        padback = ceil(TE(1)-Texp_all(1));
        delays = (ceil(diff(Texp_all(1:ind))));
        delays = fliplr(delays);
        Tth_bin2 = fliplr(Tth_bin);
        
        
        Nte = length(Texp_all);
        
        qlb_ind = find(q_exp_all == qlb);
        qub_ind = find(q_exp_all == qub);
        
        Integrated_Iexp_all = zeros(1, Nte);
        
        for ts = 1:Nte
            for qq=qlb_ind:qub_ind
                Integrated_Iexp_all(ts) = Integrated_Iexp_all(ts) + Iexp_all(qq, ts) *dq ;
            end
        end
        
        qlb_ind = find(qAng == qlb);
        qub_ind = find(qAng == qub);
        
        Integrated_Ith_all = zeros(Ntraj, Nts);
        
        for tr=1:Ntraj
            for ts=1:Nts
                for qq=qlb_ind:qub_ind
                    Integrated_Ith_all(tr, ts) = Integrated_Ith_all(tr,ts) + binned_signal(tr,qq,ts) *dq;
                end
            end
        end
        
        Integrated_Ith_all2 = fliplr(Integrated_Ith_all);
        Ith_Integrated = fliplr(Ith_Integrated);
        count = 0;
        for i=length(Tth_bin)+1:length(Tth_bin)+length(delays)
            count = count + 1;
            Tth_bin2(i) = Tth_bin2(i-1) - delays(count);
            Integrated_Ith_all2(:, i) = 0;
            Ith_Integrated(i) = 0;
        end
         
        Tth_bin = fliplr(Tth_bin2);
        Integrated_Ith_all = fliplr(Integrated_Ith_all2);
        Ith_Integrated = fliplr(Ith_Integrated);
        
        
        Integrated_Ith_all = Integrated_Ith_all*(100*exfrac_final);
        Ith_Integrated = Ith_Integrated*(100*exfrac_final)*dq;
        Inegrated_Iexp_all = Integrated_Iexp_all;
        
        
        Ith_std = std(Integrated_Ith_all, 1, 1);
        
        x2 = [Tth_bin, fliplr(Tth_bin)];
        pos_err = Ith_Integrated + Ith_std;
        neg_err = Ith_Integrated - Ith_std;
        inbetween = [neg_err, fliplr(pos_err)];

        lw = 1;
        ms = 2;
        pw = 3.65;
        ph = 2.35;
        
        f = figure;
        f.PaperUnits = 'inches';
        f.Units = 'inches';
        f.PaperPosition = [0 0 pw ph];
        f.PaperSize = [pw ph];
        f.OuterPosition = [0 0 pw ph];
        f.InnerPosition = [0 0 pw ph];
        
        
        plot(Tth_bin, Ith_Integrated, '-or', 'MarkerFaceColor', 'r', 'LineWidth', lw, 'MarkerSize', ms)
        hold on
        fill(x2, inbetween, 'r', 'FaceAlpha', 0.2)
        plot(Tth_bin, Integrated_Iexp_all, '-ob', 'MarkerFaceColor', 'b', 'LineWidth', lw, 'MarkerSize', ms);
        ax = gca;
        ax.FontSize = 6;
        ax.XAxis.MinorTick = 'on';
        
       
        xlabel('$t$ (fs)', 'interpreter', 'latex', 'FontSize', 10)
        ylabel('$\%\Delta I^{\mathrm{int}}$', 'interpreter', 'latex', 'FontSize', 10);
        xline(0, '--k')
        axis tight
        %title(['T0 = ', num2str(T0_exp), ' fs , Pulse Width = ', num2str(sum(pulse)), ' fs'])
        legend('Theory ($w_i = 1/N_{\mathrm{trj}}$)', 'Std. dev. (theory)', 'Experiment', 'interpreter', 'latex', 'Location', 'NorthWest', 'FontSize', 6)
        pt = ['T0FIXED_TE_', num2str(round(T0_exp)), '_Conv_', num2str(sum(pulse)), '.fig'];
        exportgraphics(f, 'T0FIT.pdf', 'BackgroundColor', 'none', 'ContentType', 'vector')
        keyboard
        
        ind = find(Texp_all == T0_exp);
        TEp = Texp_all(1:ind-8);
        Integrated_Iexpp = Integrated_Iexp_all(1:ind-8);
        
        p = polyfit(TEp, Integrated_Iexpp, 1);
        p1 = polyval(p, TEp);
   
        figure      

        plot(TEp, Integrated_Iexpp, '-ob', 'MarkerFaceColor', 'b')
        axis tight
        hold on
        plot(TEp, p1, '-r')
        legend('Experiment', 'Fit', 'Location', 'SouthWest')
        
        
    end
end


colors = [254 196 79; 236 112 20; 153 52 4]./255;

[P, T] = meshgrid(pulses, T0s);

xfracs = xfracs*100;

lw = 1.5;
dx = 0.025;
np = 2;
alphabet = ('a':'z').';
chars = num2cell(alphabet(1:np));
chars = chars.';
charlbl = strcat('(',chars,')');

f1 = figure;
h = gobjects(np); 
tl = tiledlayout(np, 1, 'TileSpacing','compact', 'Padding','normal');
tl.Units = 'inches';

ax{1} = nexttile;
%ax{1}.OuterPosition = [0.25 0.25 3 3];
s1 = plot(P(1,:), xfracs(1,:), '-o', 'Color', colors(1,:), 'MarkerFaceColor', colors(1,:), 'LineWidth', lw);
xticks([150 170 190 210 230 250]);
xticklabels({'150', '170', '190', '210', '230', '250'});
hold on 
s2 = plot(P(1,:), xfracs(2,:), '-s', 'Color', colors(2,:), 'MarkerFaceColor', colors(2,:), 'LineWidth', lw);
xticks([150 170 190 210 230 250]);
xticklabels({'150', '170', '190', '210', '230', '250'});
hold on 
s3 = plot(P(1,:), xfracs(3,:), '-^', 'Color', colors(3,:), 'MarkerFaceColor', colors(3,:), 'LineWidth', lw);
xticks([150 170 190 210 230 250]);
xticklabels({'150', '170', '190', '210', '230', '250'});
axis tight
ylabel('Excitation Fraction (%)', 'FontSize', 10)
ylim([2 6])
yticks(0:1:6)
ax{1}.YAxis.MinorTick = 'on';
ax{1}.YAxis.MinorTickValues = 0:0.25:6;
xlabel('Gaussian Width (fs)', 'FontSize', 10)

ax{2} = nexttile;
%ax{2}.OuterPosition = [0.25 0.25 3 3];
s4 = plot(P(1,:), errf(1,:), '-o', 'Color', colors(1,:), 'MarkerFaceColor', colors(1,:), 'LineWidth', lw);
xticks([150 170 190 210 230 250]);
xticklabels({'150', '170', '190', '210', '230', '250'});
hold on 
s5 = plot(P(1,:), errf(2,:), '-s', 'Color', colors(2,:), 'MarkerFaceColor', colors(2,:), 'LineWidth', lw);
xticks([150 170 190 210 230 250]);
xticklabels({'150', '170', '190', '210', '230', '250'});
hold on 
s6 = plot(P(1,:), errf(3,:), '-^', 'Color', colors(3,:), 'MarkerFaceColor', colors(3,:), 'LineWidth', lw);
xticks([150 170 190 210 230 250]);
xticklabels({'150', '170', '190', '210', '230', '250'});
axis tight
ylabel('Target Function Value', 'FontSize', 10)
ylim([0 10])
yticks(0:2.5:10)
ax{1}.YAxis.MinorTick = 'on';
ax{1}.YAxis.MinorTickValues = 0:1.25:10;
xlabel('Gaussian Width (fs)', 'FontSize', 10)

lg = legend(ax{2}, '-20 fs', '+30 fs', '+80 fs', 'Orientation','Horizontal','NumColumns',3);
lg.Layout.Tile = 'South';
lg.FontSize = 12;
set(lg, 'Box', 'off');


