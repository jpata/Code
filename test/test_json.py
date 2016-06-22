import ROOT, json, sys
ROOT.gSystem.Load("libFWCoreFWLite.so")
ROOT.gSystem.Load("libTTHMEIntegratorStandalone.so")
from ROOT import MEM
from ROOT import TLorentzVector
import math

inf = sys.stdin

CvectorPermutations = getattr(ROOT, "std::vector<MEM::Permutations::Permutations>")
CvectorPSVar = getattr(ROOT, "std::vector<MEM::PSVar::PSVar>")

pvec = CvectorPermutations()
pvec.push_back(MEM.Permutations.BTagged)
pvec.push_back(MEM.Permutations.QUntagged)
pvec.push_back(MEM.Permutations.QQbarBBbarSymmetry)

cfg1 = ROOT.MEM.MEMConfig()
cfg1.name = "CASE1"
cfg1.defaultCfg()
cfg1.transfer_function_method = MEM.TFMethod.Builtin
cfg1.perm_pruning = pvec


cfg2 = ROOT.MEM.MEMConfig()
cfg2.name = "CASE2"
cfg2.defaultCfg()
cfg2.transfer_function_method = MEM.TFMethod.Builtin
cfg2.perm_pruning = pvec
cfg2.int_code += MEM.IntegrandType.AdditionalRadiation

confs = [cfg1, cfg2]
#t1 = ROOT.TF1("fb","[0]*exp(-0.5*((x-[1])/[2])**2)",0,500)
#t1.SetParameter(0, 1.0 / 100/math.sqrt(2*3.1415)) #normalization
#t1.SetParameter(1, 100) #mean
#t1.SetParameter(2, math.sqrt(2)*100) #unc

def add_obj(mem, typ, **kwargs):

    if kwargs.has_key("p4s"):
        pt, eta, phi, mass = kwargs.pop("p4s")
        v = TLorentzVector()
        v.SetPtEtaPhiM(pt, eta, phi, mass);
    elif kwargs.has_key("p4c"):
        v = TLorentzVector(*kwargs.pop("p4c"))
    obsdict = kwargs.pop("obsdict", {})

    o = MEM.Object(v, typ)

    #t1 = kwargs.get("tf", None)
    #if t1 != None:
    #    o.addTransferFunction(MEM.TFType.qReco, t1)
    #    o.addTransferFunction(MEM.TFType.bReco, t1)

    #    o.addTransferFunction(MEM.TFType.qLost, t1)
    #    o.addTransferFunction(MEM.TFType.bLost, t1)

    for k, v in obsdict.items():
        o.addObs(k, v)
    mem.push_back_object(o)

events = map(json.loads, inf.readlines())

for cfg in confs:

    mem = MEM.Integrand(
        MEM.output,
        cfg
    )
    for jsev in events:
        print "----"
        jets_p4 = jsev["input"]["selectedJetsP4"]
        jets_csv = jsev["input"]["selectedJetsCSV"]
        jets_btag = jsev["input"]["selectedJetsBTag"]

        njet = 0
        for p4, btag in zip(jets_p4, jets_btag):
            njet += 1
            add_obj(mem,
                MEM.ObjectType.Jet,
                p4s=p4,
                obsdict={MEM.Observable.BTAG: btag},
                #tf=t1
            )
            if njet >= 8:
                print "stopping after {0} jets".format(njet)
                break

        if njet <= 6:
            print "skipping event"
            continue
        leps_p4 = jsev["input"]["selectedLeptonsP4"]
        leps_charge = jsev["input"]["selectedLeptonsCharge"]
        for p4, charge in zip(leps_p4, leps_charge):
            add_obj(mem,
                MEM.ObjectType.Lepton,
                p4s=p4,
                obsdict={MEM.Observable.CHARGE: charge},
            )

        print "event with {0} jets, {1} leptons".format(njet, len(leps_p4))

        add_obj(mem,
            MEM.ObjectType.MET,
            p4s=(jsev["input"]["metP4"][0], 0, jsev["input"]["metP4"][1], 0),
        )
        vars_to_integrate   = CvectorPSVar()
        vars_to_marginalize = CvectorPSVar()
        if jsev["output"]["mem_cfg"] == "SL_2w2h2t":
            r1 = mem.run(MEM.FinalState.LH, MEM.Hypothesis.TTH, vars_to_integrate, vars_to_marginalize, 1000)
            r2 = mem.run(MEM.FinalState.LH, MEM.Hypothesis.TTBB, vars_to_integrate, vars_to_marginalize, 1000)
            print "{0} ev={1} tth={2} ttbb={3} val={4}".format(
                cfg.name,
                jsev["event"]["event"],
                r1.p,
                r2.p,
                r1.p/(r1.p + 0.1*r2.p)
            )
        mem.next_event()
    del mem
