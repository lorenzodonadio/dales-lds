# DALES - Dutch Atmospheric Large Eddy Simulation

[![DOI](https://zenodo.org/badge/32735454.svg)](https://zenodo.org/badge/latestdoi/32735454)

## Documentation
The following documents are included in the DALES repository (but are not fully up to date):

* [User manual](https://github.com/dalesteam/dales/blob/master/utils/doc/input/dales-manual.pdf)
* Documentation of the configuration options [Namoptions.pdf](https://github.com/dalesteam/dales/blob/master/utils/doc/input/Namoptions.pdf)
* [readme.pdf](https://github.com/dalesteam/dales/blob/master/utils/doc/input/readme.pdf)
* [INSTALL.md](https://github.com/dalesteam/dales/blob/master/INSTALL.md) contains installation instructions
* The [DALES Wiki](https://github.com/dalesteam/dales/wiki/) contains [Installation instructions](https://github.com/dalesteam/dales/wiki/Installation-notes) for various systems and notes for [visualizing and processing](https://github.com/dalesteam/dales/wiki/Visualizing-and-processing-output) the model output

## License
DALES is made available under the terms of the GNU GPL version 3, see the file COPYING for details.

## References
Formulation of the Dutch Atmospheric Large-Eddy Simulation (DALES) and overview of its applications, T. Heus et al, [Geosci. Model Dev., 3, 415-444, 2010](https://doi.org/10.5194/gmd-3-415-2010)

Ouwersloot, H.G., Moene, A.F., Attema, J.J. et al. Large-Eddy Simulation Comparison of Neutral Flow Over a Canopy: Sensitivities to Physical and Numerical Conditions, and Similarity to Other Representations. [Boundary-Layer Meteorol 162, 71–89 (2017)](https://doi.org/10.1007/s10546-016-0182-5)

### CHANGES by Lorenzo Donadio

check this [PR](https://github.com/lorenzodonadio/dales-lds/pull/1)

- Add `modgenstat_lite.f90`
    this new module write only the following four profiles to a netcdf file, a reduced version of modgenstat.f90 module
    
    `rhof` : Full level slab averaged density - kg/m^3 - tt
    `rhobf` : Full level base-state density - kg/m^3 - tt
    `rhobh` : Half level base-state density - kg/m^3 - mt
    `presh` : Pressure at cell center - Pa - tt
    
    It write to a file with this format:
    
    `profiles_lite.xxx.nc`
    
    controlled via the following namelist:

    `namelist/NAMGENSTATLITE/ timeav,lstat`


- Add domain rotation to the Coriolis implementation, in `modforces.f90` controlled via `xyrot` in `DOMAIN` namespace. if xyrot is zero (default, or less than 0.1 to be precise) then there is no rotation applied to the coriolis force. 

- Add ekh and ekm fielddump option `modfielddump.f90`, controlled via `lekh` and `lekm` flags in `NAMFIELDDUMP`

- Refactored `modfielddump.f90` to separate netcdf and binary outputs into subroutines

- Add timeaverage fielddump save option to `modfielddump.f90` controlled by `lsavetimeavg` flag in `NAMEFIELDDUMP`