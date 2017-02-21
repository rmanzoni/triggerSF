#!/usr/bin/env python
import ROOT
import numpy as np

ROOT.TH1.SetDefaultSumw2()
ROOT.gStyle.SetOptStat(False)

class SFReader():
    '''
    Class to get trigger SF for lep-tau channels.
    Reads efficiency TGraphs.
    '''
    
    def __init__(self, file, interpolate=False):
        self.interpolate = interpolate
        self.file = ROOT.TFile.Open(file, 'read')
    
    def prepareGraph(self, graph):
        '''
        Remove empty bins at the end of the spectrum
        '''
        npoints = graph.GetN()
        X = np.array([graph.GetX()[i] for i in range(npoints)])
        Y = np.array([graph.GetY()[i] for i in range(npoints)])
        ylast = Y[-1]
        jj = 1
        while ylast == 0.:
            jj += 1
            ylast = Y[-jj]
        if jj>1:
            X = X[:-jj+1]
            Y = Y[:-jj+1]
            graph = ROOT.TGraph(len(X), X, Y)
        return graph
    
    
    def getY(self, pt, graph):
        '''
        Get Y given X.
        Can be either binned or linearly etrapolated between bins, depending
        on how the class was istantiated.
        '''
        npoints = graph.GetN()
        X = np.array([graph.GetX()[i] for i in range(npoints)])
        Y = np.array([graph.GetY()[i] for i in range(npoints)])
        Xedges = [0.]
        for i in range(npoints):
            Xedges.append(X[i]+(X[i] - Xedges[-1]))
        Xedges = np.array(Xedges)
        ibin = max(np.where(pt>=Xedges)[0])
        # allow interpolation only if it's not for the last bin
        if self.interpolate and ibin<len(X):
            y = max(min(graph.Eval(pt), 1.), 0.)
        elif ibin<len(X):
            y = Y[ibin]
        else:
            y = Y[-1]
        return y    
    
    def _getWeight(self, isdata, tau_pt, tau_eta, tau_isocut='MediumIso', genuine=True, tau_dm=None):
        self.file.cd()
        istau = 'genuine'*genuine + 'fake'*(not genuine)
        eta = 'barrel'*(abs(tau_eta)<1.5) + 'endcap'*(abs(tau_eta)>=1.5)     
        data = 'data'*isdata + 'mc'*(not isdata)   
        tojoin = [data, istau, eta, tau_isocut]
        if tau_dm is not None:
            tojoin.append('dm%d'%tau_dm)
        name = '_'.join(tojoin)
        graph = self.file.Get(name)
        if not graph:
            tojoin = [data, istau, eta, tau_isocut]
            name = '_'.join(tojoin)
            graph = self.file.Get(name)
        #print name
        graph = self.prepareGraph(graph)
        return self.getY(tau_pt, graph)      
    
    def getDataWeight(self, tau_pt, tau_eta, tau_isocut='MediumIso', genuine=True, tau_dm=None):
        return self._getWeight(True, tau_pt, tau_eta, tau_isocut, genuine, tau_dm)

    def getMCWeight(self, tau_pt, tau_eta, tau_isocut='MediumIso', genuine=True, tau_dm=None):
        return self._getWeight(False, tau_pt, tau_eta, tau_isocut, genuine, tau_dm)
        
    def getSF(self, tau_pt, tau_eta, tau_isocut='MediumIso', genuine=True, tau_dm=None):
        data = self.getDataWeight(tau_pt, tau_eta, tau_isocut, genuine, tau_dm)
        mc = self.getMCWeight(tau_pt, tau_eta, tau_isocut, genuine, tau_dm)
        return data / max(1.e-30, mc)


if __name__ == '__main__':
    file = 'ele-tau/trigger_sf_et.root'
#     file = 'mu-tau/trigger_sf_mt.root'
    
    reader = SFReader(file, interpolate=False)

    taus = [
        (50.0303226283    , -0.690436537219, 10),
        (59.5032264387    , -1.95157183364 ,  0),
        (43.1410575009    ,  0.801236923794,  1),
        (35.2404971103    ,  1.89180734956 ,  0),
        (45.5487299326    , -1.33507357888 ,  1),
        (49.3098431225    , -0.55660978639 ,  1),
        (24.3615004894    ,  0.786721853764,  1),
        (59.8846279899    , -0.132008199188,  0),
        (30.2513219915    , -1.56642180495 , 10),
        (31.4383629152    ,  1.61171893841 ,  1),
        (400.4383629152   ,  1.61171893841 ,  1),
        (500.4383629152   ,  1.61171893841 ,  1),
        (10000.4383629152 ,  1.61171893841 ,  1),
    ]

    genuine = False
    
    for tau in taus:
        data = reader.getDataWeight(tau[0], tau[1], tau_isocut='VLooseIso', genuine=genuine, tau_dm=tau[2])
        mc = reader.getMCWeight    (tau[0], tau[1], tau_isocut='VLooseIso', genuine=genuine, tau_dm=tau[2])
        sf = reader.getSF          (tau[0], tau[1], tau_isocut='VLooseIso', genuine=genuine, tau_dm=tau[2])
        
        print 'data eff %.4f\t MC eff %.4f\t SF %.4f' %(data, mc, sf)

    import sys
    sys.exit(0)

    print '========================================'
    print '==== testing, only in my local area ===='
    print '========================================\n\n\n'


    filenames = [
        '/afs/cern.ch/work/m/manzoni/diTau2015/CMSSW_8_0_25/src/CMGTools/H2TauTau/cfgPython/single_tau/doneMC/DYJetsToLL_M50_LO_ext/TauTreeProducer/skim.root',
        '/afs/cern.ch/work/m/manzoni/diTau2015/CMSSW_8_0_25/src/CMGTools/H2TauTau/cfgPython/single_tau/doneMC/QCD_HT1000to1500/TauTreeProducer/skim.root',
        '/afs/cern.ch/work/m/manzoni/diTau2015/CMSSW_8_0_25/src/CMGTools/H2TauTau/cfgPython/single_tau/doneMC/QCD_HT1000to1500_ext/TauTreeProducer/skim.root',
        '/afs/cern.ch/work/m/manzoni/diTau2015/CMSSW_8_0_25/src/CMGTools/H2TauTau/cfgPython/single_tau/doneMC/QCD_HT100to200/TauTreeProducer/skim.root',
        '/afs/cern.ch/work/m/manzoni/diTau2015/CMSSW_8_0_25/src/CMGTools/H2TauTau/cfgPython/single_tau/doneMC/QCD_HT1500to2000/TauTreeProducer/skim.root',
        '/afs/cern.ch/work/m/manzoni/diTau2015/CMSSW_8_0_25/src/CMGTools/H2TauTau/cfgPython/single_tau/doneMC/QCD_HT2000toInf/TauTreeProducer/skim.root',
        '/afs/cern.ch/work/m/manzoni/diTau2015/CMSSW_8_0_25/src/CMGTools/H2TauTau/cfgPython/single_tau/doneMC/QCD_HT2000toInf_ext/TauTreeProducer/skim.root',
        '/afs/cern.ch/work/m/manzoni/diTau2015/CMSSW_8_0_25/src/CMGTools/H2TauTau/cfgPython/single_tau/doneMC/QCD_HT200to300/TauTreeProducer/skim.root',
        '/afs/cern.ch/work/m/manzoni/diTau2015/CMSSW_8_0_25/src/CMGTools/H2TauTau/cfgPython/single_tau/doneMC/QCD_HT200to300_ext/TauTreeProducer/skim.root',
        '/afs/cern.ch/work/m/manzoni/diTau2015/CMSSW_8_0_25/src/CMGTools/H2TauTau/cfgPython/single_tau/doneMC/QCD_HT300to500/TauTreeProducer/skim.root',
        '/afs/cern.ch/work/m/manzoni/diTau2015/CMSSW_8_0_25/src/CMGTools/H2TauTau/cfgPython/single_tau/doneMC/QCD_HT300to500_ext/TauTreeProducer/skim.root',
        '/afs/cern.ch/work/m/manzoni/diTau2015/CMSSW_8_0_25/src/CMGTools/H2TauTau/cfgPython/single_tau/doneMC/QCD_HT500to700/TauTreeProducer/skim.root',
        '/afs/cern.ch/work/m/manzoni/diTau2015/CMSSW_8_0_25/src/CMGTools/H2TauTau/cfgPython/single_tau/doneMC/QCD_HT500to700_ext/TauTreeProducer/skim.root',
        '/afs/cern.ch/work/m/manzoni/diTau2015/CMSSW_8_0_25/src/CMGTools/H2TauTau/cfgPython/single_tau/doneMC/QCD_HT700to1000/TauTreeProducer/skim.root',
        '/afs/cern.ch/work/m/manzoni/diTau2015/CMSSW_8_0_25/src/CMGTools/H2TauTau/cfgPython/single_tau/doneMC/QCD_HT700to1000_ext/TauTreeProducer/skim.root',
        '/afs/cern.ch/work/m/manzoni/diTau2015/CMSSW_8_0_25/src/CMGTools/H2TauTau/cfgPython/single_tau/doneMC/QCD_Pt120to170/TauTreeProducer/skim.root',
        '/afs/cern.ch/work/m/manzoni/diTau2015/CMSSW_8_0_25/src/CMGTools/H2TauTau/cfgPython/single_tau/doneMC/QCD_Pt15to30/TauTreeProducer/skim.root',
        '/afs/cern.ch/work/m/manzoni/diTau2015/CMSSW_8_0_25/src/CMGTools/H2TauTau/cfgPython/single_tau/doneMC/QCD_Pt170to300/TauTreeProducer/skim.root',
        '/afs/cern.ch/work/m/manzoni/diTau2015/CMSSW_8_0_25/src/CMGTools/H2TauTau/cfgPython/single_tau/doneMC/QCD_Pt170to300_ext/TauTreeProducer/skim.root',
        '/afs/cern.ch/work/m/manzoni/diTau2015/CMSSW_8_0_25/src/CMGTools/H2TauTau/cfgPython/single_tau/doneMC/QCD_Pt50to80/TauTreeProducer/skim.root',
        '/afs/cern.ch/work/m/manzoni/diTau2015/CMSSW_8_0_25/src/CMGTools/H2TauTau/cfgPython/single_tau/doneMC/QCD_Pt80to120/TauTreeProducer/skim.root',
        '/afs/cern.ch/work/m/manzoni/diTau2015/CMSSW_8_0_25/src/CMGTools/H2TauTau/cfgPython/single_tau/doneMC/WJetsToLNu_LO/TauTreeProducer/skim.root',
        '/afs/cern.ch/work/m/manzoni/diTau2015/CMSSW_8_0_25/src/CMGTools/H2TauTau/cfgPython/single_tau/doneMC/WpTauNuM1000/TauTreeProducer/skim.root',
        '/afs/cern.ch/work/m/manzoni/diTau2015/CMSSW_8_0_25/src/CMGTools/H2TauTau/cfgPython/single_tau/doneMC/WpTauNuM1200/TauTreeProducer/skim.root',
        '/afs/cern.ch/work/m/manzoni/diTau2015/CMSSW_8_0_25/src/CMGTools/H2TauTau/cfgPython/single_tau/doneMC/WpTauNuM1400/TauTreeProducer/skim.root',
        '/afs/cern.ch/work/m/manzoni/diTau2015/CMSSW_8_0_25/src/CMGTools/H2TauTau/cfgPython/single_tau/doneMC/WpTauNuM1600/TauTreeProducer/skim.root',
        '/afs/cern.ch/work/m/manzoni/diTau2015/CMSSW_8_0_25/src/CMGTools/H2TauTau/cfgPython/single_tau/doneMC/WpTauNuM1800/TauTreeProducer/skim.root',
        '/afs/cern.ch/work/m/manzoni/diTau2015/CMSSW_8_0_25/src/CMGTools/H2TauTau/cfgPython/single_tau/doneMC/WpTauNuM2000/TauTreeProducer/skim.root',
        '/afs/cern.ch/work/m/manzoni/diTau2015/CMSSW_8_0_25/src/CMGTools/H2TauTau/cfgPython/single_tau/doneMC/WpTauNuM2200/TauTreeProducer/skim.root',
        '/afs/cern.ch/work/m/manzoni/diTau2015/CMSSW_8_0_25/src/CMGTools/H2TauTau/cfgPython/single_tau/doneMC/WpTauNuM2400/TauTreeProducer/skim.root',
        '/afs/cern.ch/work/m/manzoni/diTau2015/CMSSW_8_0_25/src/CMGTools/H2TauTau/cfgPython/single_tau/doneMC/WpTauNuM2600/TauTreeProducer/skim.root',
        '/afs/cern.ch/work/m/manzoni/diTau2015/CMSSW_8_0_25/src/CMGTools/H2TauTau/cfgPython/single_tau/doneMC/WpTauNuM2800/TauTreeProducer/skim.root',
        '/afs/cern.ch/work/m/manzoni/diTau2015/CMSSW_8_0_25/src/CMGTools/H2TauTau/cfgPython/single_tau/doneMC/WpTauNuM3000/TauTreeProducer/skim.root',
        '/afs/cern.ch/work/m/manzoni/diTau2015/CMSSW_8_0_25/src/CMGTools/H2TauTau/cfgPython/single_tau/doneMC/WpTauNuM3200/TauTreeProducer/skim.root',
        '/afs/cern.ch/work/m/manzoni/diTau2015/CMSSW_8_0_25/src/CMGTools/H2TauTau/cfgPython/single_tau/doneMC/WpTauNuM3400/TauTreeProducer/skim.root',
        '/afs/cern.ch/work/m/manzoni/diTau2015/CMSSW_8_0_25/src/CMGTools/H2TauTau/cfgPython/single_tau/doneMC/WpTauNuM3600/TauTreeProducer/skim.root',
        '/afs/cern.ch/work/m/manzoni/diTau2015/CMSSW_8_0_25/src/CMGTools/H2TauTau/cfgPython/single_tau/doneMC/WpTauNuM3800/TauTreeProducer/skim.root',
        '/afs/cern.ch/work/m/manzoni/diTau2015/CMSSW_8_0_25/src/CMGTools/H2TauTau/cfgPython/single_tau/doneMC/WpTauNuM400/TauTreeProducer/skim.root',
        '/afs/cern.ch/work/m/manzoni/diTau2015/CMSSW_8_0_25/src/CMGTools/H2TauTau/cfgPython/single_tau/doneMC/WpTauNuM4000/TauTreeProducer/skim.root',
        '/afs/cern.ch/work/m/manzoni/diTau2015/CMSSW_8_0_25/src/CMGTools/H2TauTau/cfgPython/single_tau/doneMC/WpTauNuM4200/TauTreeProducer/skim.root',
        '/afs/cern.ch/work/m/manzoni/diTau2015/CMSSW_8_0_25/src/CMGTools/H2TauTau/cfgPython/single_tau/doneMC/WpTauNuM4400/TauTreeProducer/skim.root',
        '/afs/cern.ch/work/m/manzoni/diTau2015/CMSSW_8_0_25/src/CMGTools/H2TauTau/cfgPython/single_tau/doneMC/WpTauNuM4600/TauTreeProducer/skim.root',
        '/afs/cern.ch/work/m/manzoni/diTau2015/CMSSW_8_0_25/src/CMGTools/H2TauTau/cfgPython/single_tau/doneMC/WpTauNuM4800/TauTreeProducer/skim.root',
        '/afs/cern.ch/work/m/manzoni/diTau2015/CMSSW_8_0_25/src/CMGTools/H2TauTau/cfgPython/single_tau/doneMC/WpTauNuM5000/TauTreeProducer/skim.root',
        '/afs/cern.ch/work/m/manzoni/diTau2015/CMSSW_8_0_25/src/CMGTools/H2TauTau/cfgPython/single_tau/doneMC/WpTauNuM5200/TauTreeProducer/skim.root',
        '/afs/cern.ch/work/m/manzoni/diTau2015/CMSSW_8_0_25/src/CMGTools/H2TauTau/cfgPython/single_tau/doneMC/WpTauNuM5400/TauTreeProducer/skim.root',
        '/afs/cern.ch/work/m/manzoni/diTau2015/CMSSW_8_0_25/src/CMGTools/H2TauTau/cfgPython/single_tau/doneMC/WpTauNuM5600/TauTreeProducer/skim.root',
        '/afs/cern.ch/work/m/manzoni/diTau2015/CMSSW_8_0_25/src/CMGTools/H2TauTau/cfgPython/single_tau/doneMC/WpTauNuM5800/TauTreeProducer/skim.root',
        '/afs/cern.ch/work/m/manzoni/diTau2015/CMSSW_8_0_25/src/CMGTools/H2TauTau/cfgPython/single_tau/doneMC/WpTauNuM600/TauTreeProducer/skim.root',
        '/afs/cern.ch/work/m/manzoni/diTau2015/CMSSW_8_0_25/src/CMGTools/H2TauTau/cfgPython/single_tau/doneMC/ZpTTM1250/TauTreeProducer/skim.root',
        '/afs/cern.ch/work/m/manzoni/diTau2015/CMSSW_8_0_25/src/CMGTools/H2TauTau/cfgPython/single_tau/doneMC/ZpTTM1500/TauTreeProducer/skim.root',
        '/afs/cern.ch/work/m/manzoni/diTau2015/CMSSW_8_0_25/src/CMGTools/H2TauTau/cfgPython/single_tau/doneMC/ZpTTM1750/TauTreeProducer/skim.root',
        '/afs/cern.ch/work/m/manzoni/diTau2015/CMSSW_8_0_25/src/CMGTools/H2TauTau/cfgPython/single_tau/doneMC/ZpTTM2000/TauTreeProducer/skim.root',
        '/afs/cern.ch/work/m/manzoni/diTau2015/CMSSW_8_0_25/src/CMGTools/H2TauTau/cfgPython/single_tau/doneMC/ZpTTM2500/TauTreeProducer/skim.root',
        '/afs/cern.ch/work/m/manzoni/diTau2015/CMSSW_8_0_25/src/CMGTools/H2TauTau/cfgPython/single_tau/doneMC/ZpTTM3000/TauTreeProducer/skim.root',
        '/afs/cern.ch/work/m/manzoni/diTau2015/CMSSW_8_0_25/src/CMGTools/H2TauTau/cfgPython/single_tau/doneMC/ZpTTM3500/TauTreeProducer/skim.root',
        '/afs/cern.ch/work/m/manzoni/diTau2015/CMSSW_8_0_25/src/CMGTools/H2TauTau/cfgPython/single_tau/doneMC/ZpTTM4000/TauTreeProducer/skim.root',
        '/afs/cern.ch/work/m/manzoni/diTau2015/CMSSW_8_0_25/src/CMGTools/H2TauTau/cfgPython/single_tau/doneMC/ZpTTM500/TauTreeProducer/skim.root',
        '/afs/cern.ch/work/m/manzoni/diTau2015/CMSSW_8_0_25/src/CMGTools/H2TauTau/cfgPython/single_tau/doneMC/ZpTTM750/TauTreeProducer/skim.root',
    ]
    
    t1 = ROOT.TChain('tree')
    
    print 'loading trees ...'
    
    for i, fname in enumerate(filenames):
        print '\t%d/%d' %(i+1, len(filenames))
        if fname.endswith('.url'):
            with open(fname) as ff:
                fname = ff.readlines()[0].rstrip()
        t1.Add(fname)
    
    tot = t1.GetEntries()
    print '... loaded trees, total events:', tot

    bins = np.array([0., 10., 20., 22.5, 25., 27.5, 30., 32.5, 35., 37.5, 40., 42.5, 45., 50., 55., 60., 70., 90., 120., 200., 500., 1000.])

    h_barrel = ROOT.TH2F('tau_leg_et_barrel', 'tau_leg_et_barrel', len(bins)-1, bins, 3, 0, 3)
    h_endcap = ROOT.TH2F('tau_leg_et_endcap', 'tau_leg_et_endcap', len(bins)-1, bins, 3, 0, 3)

    h_barrel_norm = ROOT.TH2F('tau_leg_et_barrel_norm', 'tau_leg_et_barrel_norm', len(bins)-1, bins, 3, 0, 3)
    h_endcap_norm = ROOT.TH2F('tau_leg_et_endcap_norm', 'tau_leg_et_endcap_norm', len(bins)-1, bins, 3, 0, 3)
    
    for hh in [h_barrel, h_endcap, h_barrel_norm, h_endcap_norm]:
        hh.GetXaxis().SetTitle('#tau p_{T} [GeV]')
        hh.GetYaxis().SetTitle('#tau decay mode')
    
    
    for ii, tau in enumerate(t1):
        if ii%11==0:
            continue
#         if ii >= 5000:
#             break
        if ii%10000 == 0: print '===> %d/%d' %(ii, tot)
        # data = reader.getDataWeight(tau.tau_pt, tau.tau_eta, tau_isocut='VLooseIso', genuine=(abs(tau.tau_gen_pdgId) in [11, 13, 15]), tau_dm=tau.tau_decayMode)
        # mc = reader.getMCWeight    (tau.tau_pt, tau.tau_eta, tau_isocut='VLooseIso', genuine=(abs(tau.tau_gen_pdgId) in [11, 13, 15]), tau_dm=tau.tau_decayMode)
        sf = reader.getSF          (tau.tau_pt, tau.tau_eta, tau_isocut='VLooseIso', genuine=(abs(tau.tau_gen_pdgId) in [11, 13, 15]), tau_dm=tau.tau_decayMode)
        
        if abs(tau.tau_eta) < 1.5:
            h_barrel     .Fill(tau.tau_pt, tau.tau_decayMode*(tau.tau_decayMode<2) + 2*(tau.tau_decayMode==10), sf)
            h_barrel_norm.Fill(tau.tau_pt, tau.tau_decayMode*(tau.tau_decayMode<2) + 2*(tau.tau_decayMode==10), 1.)
        else:        
            h_endcap     .Fill(tau.tau_pt, tau.tau_decayMode*(tau.tau_decayMode<2) + 2*(tau.tau_decayMode==10), sf)
            h_endcap_norm.Fill(tau.tau_pt, tau.tau_decayMode*(tau.tau_decayMode<2) + 2*(tau.tau_decayMode==10), 1.)
        
        #print 'data eff %.4f\t MC eff %.4f\t SF %.4f' %(data, mc, sf)
        
    h_barrel.Divide(h_barrel_norm)
    h_endcap.Divide(h_endcap_norm)


    h_barrel.Draw('colz')
    ROOT.gPad.SaveAs('tau_leg_et_barrel.pdf')

    h_endcap.Draw('colz')
    ROOT.gPad.SaveAs('tau_leg_et_endcap.pdf')
    











