# Tau trigger scale factors

see more details here  
https://indico.cern.ch/event/605406/contributions/2482759/attachments/1414963/2165912/tau_leg_triggers2016.pdf 

* for lepton + tau channels scale factors can be retrieved using SFreader class in getSF.py
   * it computes SF from TGraphs stored in `mu-tau/trigger_sf_mt.root` or `ele-tau/trigger_sf_et.root`
   * can optionally interpolate bewteen bins
   
* for di-tau, we use an analytical function to model trigger efficiencies, see `di-tau/README.md`


TGraps and the function parameters can also be found here:
```
/afs/cern.ch/work/m/manzoni/public/rereco2016triggerSF/trigger_sf_mt.root
/afs/cern.ch/work/m/manzoni/public/rereco2016triggerSF/trigger_sf_et.root
/afs/cern.ch/work/m/manzoni/public/rereco2016triggerSF/fitresults_tt_moriond2017.json
```

