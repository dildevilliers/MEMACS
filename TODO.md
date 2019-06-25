# ToDo list for new code

# General
- [ ] write README
- [ ] include property for code version in all classes

# Farfield
- [ ] readMeasurements
- [ ] fix roundGrid to always convert grid to PhTh, round in degrees, convert back
- [ ] make a legacy struct constructor
- [ ] Speed up mirrorSymmetricPattern. Might require custom implementations of all grid types. Start with spherical though!
- [ ] Delete redundant field points in the setSymmetry functions, don't just check for them (make it optional)
- [ ] Typical pattern parameters calculator: eff, SLL, XP, Beamwidth, etc.
- [ ] Fix plotPrincipleCuts for symmetric fields. Should mirror before plotting
- [ ] writeFEKOffe
- [ ] Plan help files
- [ ] generalise the PhTh readNFSscan for repeating data - only remove what is repeated
- [ ] Extend isGrid4Pi to operate on all spherical grids directly for speed
- [ ] ReadCSTascii
- [ ] Multiple frequency concat
- [ ] expandBORpattern must also be able to handle phi=45deg cut with Ludwig3 Co-Xp fields.
- [ ] Sort out all the CP BOR1 pattern stuff (farFieldFromPowerPattern, expandBORpattern, getBORpattern)
- [ ] Make several example patterns using simple dipoles and powerPattern functions
- [ ] Check readGRASPgrd, not sure of E1 and E2 order for all cases. Pre-allocate the matrices for speed.
- [ ] Fix FarField.rotate field components pole at th = 180 (DdV)
- [ ] Field symmetries: XY plane (DdV)
- [ ] Gaussian/cosn pattern fitter
- [ ] 2/3D plots should automatically expand BOR0/1 symmetry fields before plotting
- [ ] readFITS (DdV)
- [ ] Overlap integral calculator (DdV)
- [ ] writeCSTffs
- [ ] Fix pole field sign change for setYrange gridType PhTh for the 3 cases that are not done yet.
- [ ] Test the sym/pos and 180/360 plotting order rules.  Should be X and then Y shifts always - force this in the code somehow.
- [ ] Fix AzEl and ElAz poles in getELudwig2EA and getELudwig2AE: should not be 0
- [ ] Array pattern adder
- [ ] plot on a spherical surface
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
- [x] Set coorType = power for astrogrids by default
- [x] Fix default behaviour of the 2-D plotType to have a finite step. It should run out the box with no step provided.
- [x] More robust automatic interpolation grid calculator (rotate has a basic version hardcoded)
- [x] Refactor setXgrid and setYgrid. 
- [x] include the elevation type angles for setXgrid and setYgrid

# Arrays
## General
- [ ] Make all classes frequency dependent

## PlaneWaveSignal
- [ ] freq, sigPower, sigPhase frequency dependent
- [ ] include polarization information
- [ ] include plot for the direction and polarization vectors

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
