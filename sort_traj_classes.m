function [vout] = sort_traj_classes(Q, multiplicity, FLAGtfunc)

[natom, ~, Ntraj, Nts] = size(Q);

if Ntraj ~= length(multiplicity) error('Ntraj and Multiplicity vector inconsistent'); end

switch FLAGtfunc
    case 1
        n_bound = sum(multiplicity == 0);
        n_singlet = sum(multiplicity == 1);
        n_triplet = sum(multiplicity == 3);
        Q_bound = zeros(natom, 3, n_bound, Nts);
        Q_singlet = zeros(natom, 3, n_singlet, Nts);
        Q_triplet = zeros(natom, 3, n_triplet, Nts);
        b = 0;
        s = 0;
        t = 0;
        for tr=1:Ntraj
            if multiplicity == 0
                b = b + 1;
                Q_bound(:, :, b, :) = Q(:, :, tr, :);
            elseif multiplicity == 1
                s = s + 1;
                Q_singlet(:, :, s, :) = Q(:, :, tr, :);
            elseif multiplicity == 3
                t = t + 1
                Q_triplet(:, :, t, :) = Q(:, :, tr, :);
            end
        end
        vout{1} = Q_bound;
        vout{2} = Q_singlet;
        vout{3} = Q_triplet;
    
    case 2
        n_bound = sum(multiplicity == 0);
        n_diss = sum(multiplicity ~= 0);
        Q_bound = zeros(natom, 3, n_bound, Nts);
        Q_diss = zeros(natom, 3, n_diss, Nts);
        b = 0 
        d = 0
        for tr=1:Ntraj
            if multiplicity == 0
                b = b + 1
                Q_bound(:, :, b, :) = Q(:, :, tr, :)
               
            else
                d = d + 1;
                Q_diss(:, :, d, :) = Q(:, :, tr, :);
            end
        end
        
        vout{1} = Q_bound;
        vout{2} = Q_diss;
        
end
                
        
end

