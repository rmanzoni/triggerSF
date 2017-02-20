#!/usr/bin/env python
import ROOT
import numpy as np

class SFReader():
    '''
    Class to get trigger SF for lep-tau channels.
    Reads efficiency TGraphs.
    '''
    
    def __init__(self, file, interpolate=False):
        self.interpolate = interpolate
        self.file = ROOT.TFile.Open(file, 'read')
        self.file.cd()
    
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
        if self.interpolate:
            return max(min(graph.Eval(pt), 1.), 0.)
        npoints = graph.GetN()
        X = np.array([graph.GetX()[i] for i in range(npoints)])
        Y = np.array([graph.GetY()[i] for i in range(npoints)])
        Xedges = [0.]
        for i in range(npoints):
            Xedges.append(X[i]+(X[i] - Xedges[-1]))
        Xedges = np.array(Xedges)
        ibin = max(np.where(pt>=Xedges)[0])
        y = Y[ibin-1]
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
        print name
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
    file = '/afs/cern.ch/work/m/manzoni/public/rereco2016triggerSF/trigger_sf_et.root'
    
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
        data = reader.getDataWeight(tau[0], tau[1], tau_isocut='VLoose', genuine=genuine, tau_dm=tau[2])
        mc = reader.getMCWeight    (tau[0], tau[1], tau_isocut='VLoose', genuine=genuine, tau_dm=tau[2])
        sf = reader.getSF          (tau[0], tau[1], tau_isocut='VLoose', genuine=genuine, tau_dm=tau[2])
        
        print 'data eff %.4f\t MC eff %.4f\t SF %.4f' %(data, mc, sf)
        
        
