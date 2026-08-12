# ADCPtools

ADCPtools is a MATLAB toolbox for reading, inspecting, and analysing Acoustic
Doppler Current Profiler (ADCP) data. It supports common ADCP workflows such
as coordinate transformations, backscatter calculations, quality filtering,
visualisation, repeat-transect processing, bathymetry estimation, and
regularised velocity solving.

[![View adcptools on File Exchange](https://www.mathworks.com/matlabcentral/images/matlab-file-exchange.svg)](https://nl.mathworks.com/matlabcentral/fileexchange/115160-adcptools)

## Start here

New to the toolbox? Follow the step-by-step [getting-started guide](GETTING_STARTED.md).
It includes a runnable first analysis using the sample data shipped with this
repository:

```matlab
run('examples/quick_start_vmadcp.m')
```

## Supported data and workflows

- RDI, Nortek, and Sontek ADCP data readers
- Moored and vessel-mounted ADCP datasets
- Instrument and Earth-coordinate transformations
- Navigation, heading, vertical and horizontal positioning
- Data-quality filters and diagnostic plots
- Acoustic backscatter calculations and calibration support
- Repeat-transect, bathymetry, mesh, and regularised velocity workflows

## Documentation

The MATLAB help tutorials can be opened after adding the toolbox root to the
MATLAB path:

```matlab
open_adcptools_documentation
```

The sources for these tutorials live in `doc/`; short portable examples live
in `examples/`.

## Requirements

MATLAB R2021a or newer is currently the supported baseline. See
[GETTING_STARTED.md](GETTING_STARTED.md) for installation and first-use
instructions.

## Project status

ADCPtools is actively being improved toward a more accessible, documented,
and reproducible research toolbox. Contribution, citation, release, and JOSS
submission guidance will be added in dedicated project files.

## License

ADCPtools is distributed under the GNU General Public License, version 3. See
[COPYING.txt](COPYING.txt).
