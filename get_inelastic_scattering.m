function [f_inel] = get_inelastic_scattering(q, Nq, atnum, FLAGelec)

% Calculates inelastic scattering for specified atoms over q range
% from interpolating functions of compton intensities.
% returns summed contribution over all atoms

Nat = length(atnum);

f_inel(1:Nat,1:Nq) = 0;

F = compton_intensities(); % Generate arrayfun interpolating object

for i=1:Nat
    switch atnum(i) % add more atoms according to atomic number
        case 6 
            f_inel(i,1:Nq) = F.C(q);
        case 16
            f_inel(i,1:Nq) = F.S(q);
    end
end

if FLAGelec == 1;
    f_inel = f_inel ./(q.^4);
end


f_inel = sum(f_inel, 1);

end

