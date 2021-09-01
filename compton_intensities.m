function [F] = compton_intensities()

% interpolates compton intensities over specified q range. 
% REF: International tables of crystalography vol. C section 7.4.3
% https://it.iucr.org/Cb/ch7o4v0001/sec7o4o3/

theta = [0, 0.10, 0.20 ,0.30, 0.40, 0.50, 0.60, 0.70, 0.80, 0.90, 1.00, 1.50, 2.00, 20.0]; % added in bounds (0, 20)
I_C = [0, 1.039, 2.604, 3.643, 4.184, 4.478, 4.690, 4.878, 5.051, 5.208, 5.348, 5.781, 5.930, 6]; % end = no. elec
I_S = [0, 2.151, 4.960, 6.795, 8.002, 8.960, 9.829, 10.626, 11.336, 11.952, 12.472, 13.990, 14.641, 16];

F.C = @(q) arrayfun(@(q) interp1(theta, I_C , q), q);
F.S = @(q) arrayfun(@(q) interp1(theta, I_S , q), q);

end

