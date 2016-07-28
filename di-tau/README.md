Fit parameters for HLT_DoubleMediumIsoPFTau35_Trk1_eta2p1_Reg trigger, expressed per single tau leg.  
The measurement is performed using the entire ICHEP2016 dataset, corresponding to 12.9/fb 

Efficiencies are computed wrt byIsolationMVArun2v1DBoldDMwLT.  
- 'cumulative' stands for efficiencies for event with tau iso >= WP, i.e. taus that pass medium
- 'binned' stands for efficiencies for events with tau iso == WP, i.e. taus that pass medium but not tight 

Parameters are to be passed to this fit function  
https://github.com/rmanzoni/triggertools/blob/master/objects/FitFunctions.py#L68  
that corresponds to an (approximate) convolution of a CrystalBall resolution and a step function  
