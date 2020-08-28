# ToDo list for new code

# General
- [ ] write README
- [ ] include property for code version in all classes

# Farfield
- [ ] Bug in minus. Strange results when different number of frequencies passed in objects. Issue in mathSetup
- [ ] Sort out BOR1 efficiency calculations.  Should get them from the BOR expansion, and not from normalising with the total power. 
- [ ] make light weight save/load option - just store the grid, fields, power as .mat (see below also)
- [ ] Sort out save and load to use a struct and not rely on the object itself
- [ ] finish SimpleTaper constructor, write help and testscript
- [ ] Fix and speed up mirrorSymmetricPattern. Might require custom implementations of all grid types. Start with spherical though!
- [ ] Fix rotate - used to work at some stage?  Only works for power coorType. And plenty of issues here with the base after the setRangeSph update.
- [ ] Fix the setOrientation workflow.  Does not seem to go back to [0,0,0] before doing a new orientation. Always rotates from where it is due to base settings in rotate.
- [ ] Check readGRASPgrd, not sure of E1 and E2 order for all cases. Pre-allocate the matrices for speed.
- [ ] writeFITS (DdV)
- [ ] Typical pattern parameters calculator: eff, XP, etc.
- [ ] Make a angular sub-sampler - input an approximate step size (does not interpolate)
- [ ] Make a clearer version of currentForm2Base - a method to sample the current field on a grid (resample?)
- [ ] remove static method redundant validators - the main constructor should catch the errors
- [ ] harden up E-field setters
- [ ] cannot plot DirCos version of an astrogrids base object due to empty obj.th
- [ ] fix roundGrid to always convert grid to PhTh, round in degrees, convert back
- [ ] Delete redundant field points in the setSymmetry functions, don't just check for them (make it optional)
- [ ] Fix plotPrincipleCuts for symmetric fields. Should mirror before plotting
- [ ] readSPH
- [ ] writeFEKOffe
- [ ] writeCSTffs
- [ ] Pattern getters help files
- [ ] Field and frequency setters help files
- [ ] Grid range shifter help files
- [ ] transformTypes help file
- [ ] Base grid function help files
- [ ] Plotting method help files
- [ ] interpolation method help files
- [ ] Maths overloaded methods help files
- [ ] Frequency and field modification help files
- [ ] symmetry handler help files
- [ ] Format and tester help files
- [ ] astronomical methods help files
- [ ] File output help files
- [ ] Farfield file reading methods help files
- [ ] Extend isGrid4Pi to operate on all spherical grids directly for speed
- [ ] expandBORpattern must also be able to handle phi=45deg cut with Ludwig3 Co-Xp fields.
- [ ] Sort out all the CP BOR1 pattern stuff (farFieldFromPowerPattern, expandBORpattern, getBORpattern)
- [ ] Make several example patterns using simple dipoles and powerPattern functions
- [ ] Fix FarField.rotate field components pole at th = 180 
- [ ] Field symmetries: XY plane
- [ ] Gaussian/cosn pattern fitter
- [ ] Gaussian/cosn/etc pattern constructor
- [ ] 2/3D plots should automatically expand BOR0/1 symmetry fields before plotting
- [ ] Overlap integral calculator 
- [ ] Array pattern adder
- [ ] plot on a spherical surface
- [ ] Add semilog scaleMag option for plot (might be handy for error comparison)
- [ ] Fix 3D plot for negative y-(th)axis cases
- [ ] provide 3D plot for users without Antennas toolbox
- [ ] Jones getter
- [ ] Jones plotter (fix up)
- [ ] Stokes getter
- [ ] Stokes plotter
- [ ] Resample of SWE results on different grid
- [ ] Use above resample for interpolation when plotting
- [ ] CBFP/SWE/Zernike interpolation in freq
- [ ] General CBFP/SWE/Zernike interpolation over parameters
- [ ] Zernike expansion
- [ ] Major Rework: Change angle base to degrees and not radians (for angular grids - not DirCos type grids)
- [x] make a legacy struct constructor
- [x] Small bug fix in readNFSscan
- [x] readFITS fixed up a bit
- [x] writeASCII added
- [x] added constant property: version 
- [x] readASCII added
- [x] Fixed issue in catFreq
- [x] Added getGridIndex method
- [x] Added integration rule optional argument for pradInt

# Pnt3D
- [ ] concat - must work for fused and split versions and all mixtures between them
- [x] Added the fuse and split functions
- [x] Updated constructor to take vector input
- [x] fixed plot to handle the vector input (not very hardened, but works for basic cases)

# CoordinateSystem
- [x] Added the fromYZ and fromXZ constructors


# Arrays
## General
- [ ] Make all classes frequency dependent

## PlaneWaveSignal
- [ ] freq, sigPower, sigPhase frequency dependent
- [ ] include polarization information
- [ ] include plot for the direction and polarization vectors

## ArrayImaging
- [ ] make class

## ArrayElements
- [ ] channelErrors frequency dependent
- [ ] different field patterns for each element
- [ ] include plot for the element positions

## ArrayReceiver
- [ ] Full noise coupling matrix
- [ ] LNA noise parameters specification method
- [ ] Gains frequency dependent
- [ ] More receiver (amp) characteristics - IP1/3, 3dB comp, etc

## ArrayDAC
- [ ] Signed binary
- [ ] 2s complement

## ArraySystem
- [ ] Expand plotting option for line styles etc.
- [ ] Automate the time scale to Engineering units

## ArrayBeamformer
- [ ] Implement for analog arrays
- [ ] Include transmit array pattern calculator

## ArrayDBE
- [ ] Sort out power levels in the PSD plots (Units)
- [ ] plotScanBeam must be extended to handle 1D and 2D scan plots
