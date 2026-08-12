%% ADCPtools quick start: inspect a vessel-mounted RDI deployment
% This example uses bundled sample data. It can be run from any MATLAB
% working directory and does not require changing a path for a specific user.

exampleDirectory = fileparts(mfilename('fullpath'));
toolboxRoot = fileparts(exampleDirectory);
sampleDataDirectory = fullfile(toolboxRoot, 'doc', 'sample_data', ...
    'rdi_muara_muntai_bend');

assert(isfolder(sampleDataDirectory), ...
    "ADCPtools:ExampleDataMissing", ...
    "The bundled sample data folder was not found: %s", sampleDataDirectory)

% Make the toolbox available for this MATLAB session. Do not use genpath:
% external dependencies are loaded by ADCPtools only when they are required.
addpath(toolboxRoot)

%% 1. Read raw RDI/WinRiver files
% The deployment prefix 'trans' selects all matching raw and navigation files.
rawDeployment = rdi.readDeployment('trans', sampleDataDirectory);

%% 2. Create the vessel-mounted ADCP object
adcp = rdi.VMADCP(rawDeployment);

fprintf("Loaded %d ensembles \n", adcp.nensembles)

%% 3. Perform an initial quality inspection
% Review these plots before continuing with positioning, filtering, or solving.
figure('Name', 'ADCP orientation')
adcp.plot_orientations

figure('Name', 'Ship track')
adcp.plot_track

figure('Name', 'Depth-averaged velocity')
adcp.plot_track_velocity

figure('Name', 'Velocity profiles')
adcp.plot_velocity

% Next: open the longer MATLAB tutorials with open_adcptools_documentation.
