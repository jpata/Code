import ROOT
import numpy as np
ROOT.gSystem.Load("libFWCoreFWLite.so")
ROOT.gSystem.Load("libTTHMEIntegratorStandalone.so")
from ROOT import MEM
from ROOT import TLorentzVector
import math

cfg = ROOT.MEM.MEMConfig()
cfg.defaultCfg()
cfg.save_permutations = True
cfg.transfer_function_method = MEM.TFMethod.Builtin
#cfg.int_code = sum([getattr(MEM.IntegrandType, x) for x in ["Jacobian", "Transfer", "ScattAmpl", "DecayAmpl"]])
#mem = MEM.Integrand(MEM.output+MEM.input+MEM.init+MEM.init_more+MEM.integration, cfg)
mem = MEM.Integrand(MEM.output, cfg)
print mem

def normalize_proba(vec):
    proba_vec = np.array(vec)
    proba_vec[proba_vec <= 1E-50] = 1E-50
    ret = np.array(np.log10(proba_vec), dtype="float64")
    return ret

def save_perm_hists(outfile, res):
    of = ROOT.TFile(outfile, "RECREATE")
    for iperm in range(res.num_perm):
        perm_str = ",".join([str(v) for v in res.permutation_indexes[iperm]])
      
        N = res.num_max_calls
        v_p = normalize_proba([v for v in res.permutation_probas[iperm]])
        v_p_tf = normalize_proba([v for v in res.permutation_probas_transfer[iperm]])
        v_p_me = normalize_proba([v for v in res.permutation_probas_me[iperm]])

        w = np.ones(N, dtype="float64")
        h = ROOT.TH1D("perm_p_{0}".format(iperm), "Permutation {0}: {1} proba".format(iperm, perm_str), 100, 100, 0)
        h.FillN(N, v_p, w)
        h.Write()
        
        h = ROOT.TH1D("perm_p_tf_{0}".format(iperm), "Permutation {0}: {1} transfers".format(iperm, perm_str), 100, 100, 0)
        h.FillN(N, v_p_tf, w)
        h.Write()
        
        h = ROOT.TH1D("perm_p_me_{0}".format(iperm), "Permutation {0}: {1} matrix element".format(iperm, perm_str), 100, 10, 30)
        h.FillN(N, v_p_me, w)
        h.Write()
    #of.Write()
    of.Close()

def add_obj(mem, typ, **kwargs):

    if kwargs.has_key("p4s"):
        pt, eta, phi, mass = kwargs.pop("p4c")
        v = TLorentzVector()
        v.SetPtEtaPhiM(pt, eta, phi, mass);
    elif kwargs.has_key("p4c"):
        v = TLorentzVector(*kwargs.pop("p4c"))
    obsdict = kwargs.pop("obsdict", {})

    o = MEM.Object(v, typ)

    t1 = kwargs.get("tf", None)

    if t1 != None:
        o.addTransferFunction(MEM.TFType.qReco, t1)
        o.addTransferFunction(MEM.TFType.bReco, t1)

        o.addTransferFunction(MEM.TFType.qLost, t1)
        o.addTransferFunction(MEM.TFType.bLost, t1)

    for k, v in obsdict.items():
        o.addObs(k, v)
    mem.push_back_object(o)

t1 = None
# t1 = ROOT.TF1("fb","[0]*exp(-0.5*((x-[1])/[2])**2)",0,500)
# t1.SetParameter(0, 1.0 / 100/math.sqrt(2*3.1415)) #normalization
# t1.SetParameter(1, 100) #mean
# t1.SetParameter(2, math.sqrt(2)*100) #unc

# add_obj(mem,
#     MEM.ObjectType.Jet, p4c=(50, 0, 10, math.sqrt(50*50+10*10)),
#     obsdict={MEM.Observable.BTAG: 0.0}, tf=t1
# )
# add_obj(mem,
#     MEM.ObjectType.Jet, p4c=(30, 30, 40, math.sqrt(30*30+30*30+40*40)),
#     obsdict={MEM.Observable.BTAG: 0.0}, tf=t1
# )


add_obj(mem,
    MEM.ObjectType.Jet, p4c=(0, 50, 20, math.sqrt(50*50+20*20)),
    obsdict={MEM.Observable.BTAG: 1.0}, tf=t1
)

add_obj(mem,
    MEM.ObjectType.Jet, p4c=(70, 20, 10, math.sqrt(70*70+20*20+10*10)),
    obsdict={MEM.Observable.BTAG: 1.0}, tf=t1
)
add_obj(mem,
    MEM.ObjectType.Jet, p4c=(20, 50, 10, math.sqrt(20*20+50*50+10*10)),
    obsdict={MEM.Observable.BTAG: 1.0}, tf=t1
)
add_obj(mem,
    MEM.ObjectType.Jet, p4c=(100, 10, 20, math.sqrt(100*100 + 10*10 + 20*20)),
    obsdict={MEM.Observable.BTAG: 1.0}, tf=t1
)

add_obj(mem,
    MEM.ObjectType.Lepton, p4c=(70, 10, 20, math.sqrt(70*70+20*20+10*10)),
    obsdict={MEM.Observable.CHARGE: 1.0}
)
add_obj(mem,
    MEM.ObjectType.Lepton, p4c=(70, -10, -20, math.sqrt(70*70+20*20+10*10)),
    obsdict={MEM.Observable.CHARGE: -1.0}
)
add_obj(mem,
    MEM.ObjectType.MET, p4c=(30, 0, 0, 30),
)
CvectorPermutations = getattr(ROOT, "std::vector<MEM::Permutations::Permutations>")

pvec = CvectorPermutations()
pvec.push_back(MEM.Permutations.BTagged)
pvec.push_back(MEM.Permutations.QUntagged)
cfg.perm_pruning = pvec

CvectorPSVar = getattr(ROOT, "std::vector<MEM::PSVar::PSVar>")
vars_to_integrate   = CvectorPSVar()
vars_to_marginalize = CvectorPSVar()
r = mem.run(MEM.FinalState.LL, MEM.Hypothesis.TTH, vars_to_integrate, vars_to_marginalize, 1000)

print "tth", r.p
r = mem.run(MEM.FinalState.LL, MEM.Hypothesis.TTBB, vars_to_integrate, vars_to_marginalize, 1000)
save_perm_hists("out_tth.root", r)
print "ttbb", r.p
save_perm_hists("out_ttjets.root", r)
mem.next_event()
