% adam kirrander 2015-05-14
% get scattering form factors appropriate for traj_main.m
% calculation (elastic x-ray scattering from AI-MCE trajectories)

%  INPUT
% q(1:Nq)          momentum transfer scattering vector (au^-1)
% atmnum(Nat)      atom number for each atom
% FLAGelec         0 xray scattering: return IAM tabulated form-factor
% FLAGelec         1 electron scattering: return IAM tabulated form-factor + nuclear scattering term

%  OUTPUT
% FF(Nat,Nat,Nq)   form-factor product  fq(i,1:Nq).*fq(j,1:Nq)
% fq(Nat,Nq)       elastic form-factor 

% NB form factors in inverse Angstroms

%function [FF,fq] = get_scattering_factors(q,atmnum)
function [FF,fq] = get_scattering_factors(q,atmnum,FLAGelec)
    
Nat = length(atmnum);
Nq  = length(q);

% get elastic scattering form factors (using James' code)
if max(q) > 8*pi     error('Maximum value for q too large');   end
if FLAGelec==0
    F  = f_functions;          % elegant in-line function object-based
elseif FLAGelec==1
    F  = f_functions_electron; % IAM + nuclear charge term
else
    disp(['FLAGelec=' num2str(FLAGelec)]);         
    stop('FLAGelec has invalid value');
end
fq = zeros(Nat,Nq);
% check atom-types
for i = 1:Nat
    switch atmnum(i)
      case 1 % H-atom
        fq(i,1:Nq) = F.H(q);
      case 6 % C-atom
        fq(i,1:Nq) = F.C(q);
      case 7 % N-atom
        fq(i,1:Nq) = F.N(q);
      case 9 % F-atom
        fq(i,1:Nq) = F.F(q);      
      case 16 % S-atom
        fq(i,1:Nq) = F.S(q);      
      case 53 % I-atom
        fq(i,1:Nq) = F.I(q);      
      case 54 % Xe-atom
        fq(i,1:Nq) = F.Xe(q);      
      case default
        error('Atom not parametrized in f_functions')
    end
end % end for-loop

    
% in calculations, f-factors always appear as products
FF = zeros(Nat,Nat,Nq);
for i = 1:Nat
    for j = 1:Nat
        FF(i,j,1:Nq) = fq(i,1:Nq).*fq(j,1:Nq);
    end
end


end % end function




