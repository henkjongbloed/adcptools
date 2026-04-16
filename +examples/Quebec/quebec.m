% Script to play for Ariane

%% Path management
RF = 'C:\Users\jongb013\Documents\PHD\2-Programming\'; %RootFolder
%addpath(genpath(strcat(RF,'Tools\adcptools'))); %path to ADCPTools
% addpath(genpath(strcat(RF,'Tools\adcptools'))); %possible other folders

%% Quick documentation walkthrough - comment out

%open_adcptools_documentation()

%% Loading in the data

DP = 'C:\Users\jongb013\Documents\PHD\5-Projects\Ariane\data\quebec'; %path to data
dat = rdi.readDeployment('Quebec_0_0', DP);

%% Load water level datA

load('C:\Users\jongb013\Documents\PHD\5-Projects\Ariane\data\Donnees_validation\Donnees_validation\2009\marégraphes_h_2009_HNE_NMM_3min\marégraphes_h_2009_HNE_NMM_3min\3250Lauzon2009_HNE_NMM_3min.mat')

%% Temporal interpolation: Tidal Constituents

constituents = {'M2', 'M4'};

%% waterlevel

filt = ~isnan(h);
water_level = VaryingWaterLevel(datetime(t(filt), 'ConvertFrom', 'datenum'), h(filt));
water_level.model = TidalScalarModel(constituents = constituents);
water_level.model.scalar_name = 'eta'; % Scalar
water_level.get_parameters();


%% Modify the following code to analyze the data

V = rdi.VMADCP(dat);
V.horizontal_position_provider = HorizontalPositionFromBottomTracking; % possibly modify

V.water_level_object = water_level;

B = BathymetryScatteredPoints(V);

%Bfilt = find(B.known(2,:)>0);

B.interpolator.span = .001;
figure;
B.plot

V.filters = Filter;
%V.shipvel_provider = ShipVelocityFromBT; % possibly modify


%figure;
%hold on
%V.plot_all

[ef, xs] = cross_section_selector(V);

%% Mesh for plotting

mesh_maker = SigmaZetaMeshFromVMADCP(ef, xs, B, 'NoExpand', V);

mesh = mesh_maker.get_mesh(resn = 50, resz = 15);

%% Model

opts = SolverOptions(extrapolate_vert = 0, lat_weight_factor = 10); % possibly modify
%opts.force_zero = [1 1 1 1 1];

% Empirical model: VelocityModel;
flow_model = TaylorTidalVelocityModel; % possibly modify to enter desired empirical model formulation
flow_model.constituents = constituents;

%or TaylorVelocityModel
flow_model.n_order = [1 1 1];
flow_model.s_order = [1 1 1];
flow_model.sigma_order = [1 1 1];


%Solver options and regularization
flow_regs = regularization.Velocity.get_all_regs(mesh, B, xs, flow_model, opts, 'NoExpand', V);


% Bulk regularization parameter % possibly modify
lc = 1000.0;
flow_regs(1).weight =  lc;
flow_regs(2).weight =  lc;
flow_regs(3).weight =  lc;
flow_regs(4).weight =  lc;
flow_regs(5).weight =  lc;


% Solve for the flow
flow_solv = LocationBasedVelocitySolver(mesh, B, xs, ef, flow_model, opts, 'NoExpand', V, flow_regs); 
flow_solv.rotation = xs.angle;
flow = flow_solv.get_solution(); % possibly modify
%figure
%spy(flow.M)
% Plot the state vector
flow.plot_solution()