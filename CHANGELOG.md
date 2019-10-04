# Changelog for the MEMACS public repository

# 2019-10-04, v0.4
## FarField
- Major change to the base grid. Only stored when needed (after changes in grids/coor/pol)
- Fixed setRangeSph. Works for all the spherical grids in their native forms.
- 2D plot and interpolation default grids improved to speed up plotting
- Improved performance for very large field grids
- Fixed AzEl and ElAz poles 
- Added beamwidth and SLL calculator
- Added readMeasurements and readCSTascii
- Added multiple frequency concat function
- Added optional struct input argument that keep all the optional arguments
- Updated several help files
- Updated the testscript

# 2019-06-29, v0.3
## Arrays
- Cleanup and add some calibration functionality

# 2019-05-14, v0.2
## FarField
- Added: Astrogrids, coorType = 'power'
- Updated the astronomical grid functionality 
- included a lightweight 'power' coorType for fields where only the power pattern, and not the actual E-fields, are of interest.

# 2019-05-09, v0.1
## First commit
Create repository




