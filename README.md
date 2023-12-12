# pointing_calibration
Given expected and actual angles to a source, calculate the pointing errors for antenna offsets

Point the antenna where you expect Cas A to drift through the beam and log the receive power levels. When the power level is
maximized you are "on" the source. You know where Cas A is supposed to be in the sky at that time and you know where the
antenna is pointed. Do this over the whole sky to try and figure the least squares fit.

There is a problem with this process since the "on" isn't actually peaked on the source, it's just somewhere in the beam.
I've not yet set up a system for following Cas A and peaking on it, with simple X's or conical scans, since ambushing the 
star is very much easier.
