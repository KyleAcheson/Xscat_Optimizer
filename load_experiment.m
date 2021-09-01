function [Texp, Iexp, q_exp, CM] = load_experiment(fpath, FLAGconfmat)

% THIS FILE IS TO BE EDITIED AS PLEASED - LOAD THE RELEVENT DATA
% IN AND RETURN IT IN THE FORM REQUIRED. 
% IF NO CONFIDENCE MATRIX - SET TO ONES size(Iexp)

% INPUTS:
% fpath - path to file containing experimental data
% FLAGconfmat - 1 = returns actual confidence matrix used in fitting. 0 =
% returns matrix of ones.
% Confidence_Tol - Thresh at which points on the grid will be set to zero
% weight. The confidence matrix is applied to the residual error, the
% error at these indexes will be set to zero.

% OUTPUTS:
% Texp - selected time range to scan over in fit
% Iexp - selected experimental signal to fit - choose a suitable subrange to scan over
% q_exp - q (or s) scattering vector in experiment
% CM - confidence matrix - provides a confidence value for every time/ q point in the fit

load(fpath)

%% EDIT THESE VAIRBALES/ RANGES AS REQUIRED
Texp = flip(time3235);
Iexp = flip(deltai3235, 1)';
q_exp = s;


switch FLAGconfmat
    case 0
        CM = ones(size(Iexp));
    case 1
        CM = flip(distd3235, 1)';
        CM = 1./CM; % take inverse if matrix of standard deviations
        disp(["Using Confidence Matrix In OPT."])
end
        
% Texp = Texp(13:43);  % subset of experimental time range to select
% Iexp = Iexp(:, 13:43);
% CM = CM(:, 13:43);
%%

disp(["Experiment Sliced to Time Range: ", num2str(Texp(1)), " TO ", num2str(Texp(end))])

end

