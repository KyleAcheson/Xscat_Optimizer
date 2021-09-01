# Xscat_Optimizer

Fits Theoretical X-ray or electron scattering patters to experiment.
Reads in trajectories in `load_trajectories.m` and experimental signal in `load_experiment.m` - modify as needed.
For theory - a series of geometries over N trajectories, and their spin multiplicities are required. 
For experiment you need the experimental difference signal, time bins and q range.


## Types of optimisation
- minimisation of constrained multi-variable function via. interior-point or active-set algorithms.
- non-linear least squares minimisation 

## Options
- Polarisation correction (for X-ray scattering only)
- Inelastic (Compton) correction (assuming geometry independent within Independent Atom Model)
- Fits individual trajectories, or classes of trajectories based on spin multiplicity (must have labels for spin of trajs)
- Include excitation fraction in optimisation - or scan a range
- Scan time zero over range of experimental times
- Single or double convolution of theoretical signal with pump/ probe pulse
- Run in parallel
- Tune optimisation tolerances/ threshs 

## Post optimisation analysis
- Use the `post_opt_analysis.m` script to analyse your optimisation.
- Just edit the file names to loop over for each optimisation run.
- Inspect the optimisation visually using the FLAGplot option.
