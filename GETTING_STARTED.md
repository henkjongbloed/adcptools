# Getting started with ADCPtools

This guide gets a first vessel-mounted ADCP (VMADCP) dataset loaded and
inspected. It uses data included in this repository, so no external download
or local path configuration is needed.

## Requirements

- MATLAB R2021a or newer. The toolbox metadata currently targets R2021a.
- A local checkout of this repository.

Some workflows may require additional MATLAB toolboxes. Start with the
included example below; it is the quickest way to identify whether your
MATLAB installation satisfies the requirements for the basic workflow.

## Install the toolbox

Clone or download the repository, then add its root directory to the MATLAB
path. Run the following from the MATLAB Command Window, replacing the example
path with the location of your checkout:

```matlab
toolboxRoot = 'C:\\path\\to\\adcptools';
addpath(toolboxRoot)
savepath  % optional: keep ADCPtools on the path after restarting MATLAB
```

Do not use `addpath(genpath(...))`. The toolbox loads external packages only
when they are needed, and recursively adding every folder can mask MATLAB
functions or introduce version conflicts.

## Run the first tutorial

In MATLAB, open and run:

```matlab
run(fullfile(toolboxRoot, 'examples', 'quick_start_vmadcp.m'))
```

The script loads the bundled RDI/WinRiver sample deployment, creates an
`rdi.VMADCP` object, and opens a few diagnostic plots. It finds both the
toolbox and the sample data relative to its own location; no username-specific
paths need to be edited.

## Use your own data

The first step depends on instrument type and deployment configuration:

```matlab
% RDI vessel-mounted data (a deployment shares a filename prefix)
raw = rdi.readDeployment('deploymentPrefix', 'C:\\path\\to\\deployment');
adcp = rdi.VMADCP(raw);

% RDI moored data
raw = rdi.readADCP('C:\\path\\to\\measurement.PD0');
adcp = rdi.ADCP(raw);

% Nortek data
adcp = nortek.VMADCP('C:\\path\\to\\measurement.SigVM');
```

After constructing an ADCP object, start with `adcp.plot_orientations`,
`adcp.plot_velocity`, and `adcp.plot_backscatter` to inspect the data.

## Next steps

The bundled MATLAB help contains longer walkthroughs for reading data,
positioning, repeat-transect processing, mesh construction, and solving. Open
it from MATLAB with:

```matlab
open_adcptools_documentation
```

The tutorials in `doc/` are the source of that MATLAB help. The `examples/`
folder contains short, runnable entry points intended for new users.
