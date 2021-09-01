function [pdW] = theory_signal(Q, kin, atmnum, qAng, time, FLAGelec, FLAGinel, FLAGsignal)

% Calculates individual IAM signal for a series of or a single trajectory.

% INPUTS:
% Q - geometries
% exfrac - excitation fraction
% kin - incident wave vector
% atnum - atomic numbers of all atoms in molecule, in order of labels
% qAng - momentum transfer vec in inv angstrom
% time - time vector
% FLAGelec - FLAG for Xray or UED
% FLAGsignal - format of signal output - 0 = dI/I, 1 = dsM

% OUTPUTS:
% pdW - signal in format specificed by FLAGsignal 


[x,y,Ntraj,Nts] = size(Q); % assumes dimensions natom, 3, ntraj, nts

if x ~= length(atmnum); error('Trajectory dimensions do not add up. Must be [natom, 3, nts].'); end
%if Nts ~= length(time); warning('Length of time vector and trajectory time axis do not match - IGNORE IF SKIPPING TIME.'); end

%Q = Q(:,:,:,1:Nts);
tmax = time(end);
au2ang = 0.52917721092d0;
ang2au = 1/au2ang;
kin = kin*au2ang; % wave vector - cs2 UED
Nq = length(qAng); % number of points over theta range 

[FF,fq]  = get_scattering_factors(qAng,atmnum,FLAGelec); % NB: call using q in Angstroms^(-1)
Iat      = sum(fq.^2); % atomic scattering term 

if FLAGinel == 1
    f_inel = get_inelastic_scattering(qAng, Nq, atmnum, FLAGelec); % compton scattering
else
    f_inel = zeros(1,Nq);
end

Wiam_tot = zeros(Nq,Nts);
for ts=1:Nts % loop over time steps
    Wiam = zeros(Ntraj,Nq); % init scattering matrix for N trajs
    D = Q(:,:,:,ts); % geometry at each time-step
    for traj=1:Ntraj % loop over trajs
        Imol = zeros(1,Nq); % init for each traj at each time step
        sinQ = zeros(1,Nq);  
        for a=1:x
            for b=a+1:y
                DQ = qAng(1:Nq)*norm(D(a,1:3,traj)-D(b,1:3,traj)); %qRij
                sinQ(1:Nq) = sin(DQ(1:Nq))./DQ(1:Nq); % Sin(qRij)/qRij
                ind = find(abs(DQ)<1.d-9); % check q=0 limit
                sinQ(ind) = 1.d0;
                Imol(1:Nq) = Imol(1:Nq) + 2.d0*(squeeze(FF(a,b,1:Nq)).').*sinQ(1:Nq); 
            end
        end
        
        switch FLAGsignal
            case 0 % I
                Wiam(traj,1:Nq) = Imol(1:Nq) + Iat + f_inel; 
            case 1 % sM
                Wiam(traj,1:Nq) = (qAng .* Imol(1:Nq) + f_inel) ./ (Iat); %  s*Imol/Iat
        end
    end
    
    Wiam_avg = sum(Wiam,1)./Ntraj; % assumes equal weights
    Wiam_tot(1:Nq,ts) = Wiam_avg;  
end

if Nts > 1
    pdW = zeros(Nq,Nts);
    switch FLAGsignal
        case 0
            for ts=1:Nts % Delta I/ I
                pdW(1:Nq,ts) = ( Wiam_tot(1:Nq,ts) - Wiam_tot(1:Nq,1) ) ./ Wiam_tot(1:Nq,1);
            end
            
        case 1
            for ts=1:Nts % Delta sM
                pdW(1:Nq,ts) = Wiam_tot(1:Nq,ts) - Wiam_tot(1:Nq,1);
            end
    end
end


end

