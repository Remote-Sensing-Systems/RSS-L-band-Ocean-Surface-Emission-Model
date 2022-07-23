# RSS-L-band-Ocean-Surface-Emission-Model
FORTRAN90 code package for RSS L-band Ocean Surface Emission Model of the RSS/NASA Version 5 Salinity Release

This code package contains:
(1).    meissner_wentz_dielectric.f90: 
FORTRAN90 code for emission from specular ocean surface using the 
Meissner Wentz model for the dielectric constant of sea-water.

(2).    L_Band_wind_emissivity_AQV5_module.f90
FORTRAN90 code for emisison from wind roughened ocean-surface

(3).    deW_harm_coeffs_V9A_MI.dat
Binary file containing table of harmonic coefficients for emisison from wind roughened ocean-surface used in (2).

(4).    delta_EW_V5_B.dat
Binary file containing table of SST dependency adjustment for emisison from wind roughened ocean-surface used in (2).

References: 
 1.	Meissner, T.; F. Wentz and D. Le Vine, Aquarius Salinity Retrieval Algorithm Theoretical Basis Document (ATBD), 
       End of Mission Version; RSS Technical Report 120117; December 1, 2017; 
       Available online at ftp://podaac-ftp.jpl.nasa.gov/allData/aquarius/docs/v5/AQ-014-PS-0017_Aquarius_ATBD-EndOfMission.pdf.

 2.    Meissner, T, F. Wentz, and D, Le Vine, 2018, 
       The Salinity Retrieval Algorithms for the NASA Aquarius Version 5 and SMAP Version 3 Releases, Remote Sensing 10, 1121, doi:10.3390/rs10071121. 

 3.    Meissner, T. and F. Wentz, The complex dielectric constant of pure and sea water from microwave satellite observations, 
       IEEE TGRS, 2004, 42(9), 1836 – 1849, doi:10.1109/TGRS.2004.831888.

 4. 	 Meissner, T. and F. Wentz, The emissivity of the ocean surface between 6 and 90 GHz over a large range of wind speeds and Earth incidence angles, 
       IEEE TGRS, 2012, 50(8), 3004 – 3026, doi: 10.1109/TGRS.2011.2179662.   

 5. 	Meissner, T., F. Wentz, F. and L. Ricciardulli, The emission and scattering of L-band microwave radiation  from rough ocean surfaces 
       and wind speed measurements from Aquarius, J. Geophys. Res. Oceans, 2014, 119, doi:10.1002/2014JC009837.
