RF = 'C:\Users\jongb013\Documents\PHD\2-Programming\'; %RootFolder


addpath(genpath(strcat(RF,'Tools\adcptools'))); %path to ADCPTools

readnew = 1;
if readnew
    adcp = rdi.readDeployment('Rijn'); %path to .mat structs of processed data
else
    load proc_adcp.mat
end
V = rdi.VMADCP(adcp);
V.horizontal_position_provider = HorizontalPositionFromBottomTracking;
[ef, xs] = cross_section_selector(V);

%% Bathy
bathy = BathymetryScatteredPoints(ef, V);
bathy.interpolator.span = 0.01;

mesh_maker = SigmaZetaMeshFromVMADCP(ef, xs, bathy, 'NoExpand', V);
%%
mesh = mesh_maker.get_mesh(resn = 50, resz = 15);

%SolverOptions
opts = SolverOptions(extrapolate_vert = 0, lat_weight_factor = 1);
%opts.force_zero = [1 1 1 1 1];


% Empirical model: TaylorTidal
flow_model = VelocityModel;
%flow_model.s_order = [1 1 1];
%flow_model.n_order = [1 1 1];
%flow_model.sigma_order = [1 1 1];

%Solver options and regularization

% %% Add regularization terms
flow_regs = regularization.Velocity.get_all_regs(mesh, bathy, xs, flow_model, opts, 'NoExpand', V);


lc = 1;
flow_regs(1).weight =  lc;
flow_regs(2).weight =  lc;
flow_regs(3).weight =  lc;
flow_regs(4).weight =  lc;
flow_regs(5).weight =  lc;

%Solve for the flow

flow_solv = LocationBasedVelocitySolver(mesh, bathy, xs, ef, flow_model, opts, 'NoExpand', V, flow_regs); 
flow_solv.rotation = xs.angle;



flow = flow_solv.get_solution();



flow.plot_solution()

