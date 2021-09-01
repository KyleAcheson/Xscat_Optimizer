function [Q, multiplicity] = load_trajectories(fpath)

% OUTPUTS:
% Q - trajectory geometries [natom, 3, ntraj, nts]
% multiplicity - labels of traj multiplicities

load(fpath);
Q = geometries(:,:,:,:);
multiplicity = multiplicity(:);

end

