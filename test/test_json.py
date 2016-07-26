import ROOT, json, sys, os
ROOT.gSystem.Load("libFWCoreFWLite.so")
ROOT.gSystem.Load("libTTHMEIntegratorStandalone.so")
from ROOT import MEM
from ROOT import TLorentzVector
import math
import pickle

FILE_NAMES = os.environ["FILE_NAMES"].split()
SKIP_EVENTS = int(os.environ["SKIP_EVENTS"]) 
MAX_EVENTS = int(os.environ["MAX_EVENTS"]) 

inlines = []
nline = 0
for fn in FILE_NAMES:
    fi = open(fn)
    for line in fi.readlines():
        nline += 1
        if nline <= SKIP_EVENTS:
            continue
        if len(inlines) >= MAX_EVENTS:
            break
        inlines += [line]
    fi.close()

print "processing {0} events".format(len(inlines))
CvectorPermutations = getattr(ROOT, "std::vector<MEM::Permutations::Permutations>")
CvectorPSVar = getattr(ROOT, "std::vector<MEM::PSVar::PSVar>")

pvec = CvectorPermutations()
pvec.push_back(MEM.Permutations.BTagged)
pvec.push_back(MEM.Permutations.QUntagged)
pvec.push_back(MEM.Permutations.QQbarBBbarSymmetry)

cfg1 = ROOT.MEM.MEMConfig()
cfg1.name = "CASE1"
cfg1.defaultCfg()
cfg1.points = 2000
cfg1.save_permutations = True
cfg1.vars_to_integrate = CvectorPSVar()
cfg1.vars_to_marginalize = CvectorPSVar()
cfg1.transfer_function_method = MEM.TFMethod.External
cfg1.perm_pruning = pvec
cfg1.do_add_jet = lambda pt, btag: True
cfg1.do_calculate = lambda njet, nlep: njet>=6 and nlep==1

cfg2 = ROOT.MEM.MEMConfig()
cfg2.name = "CASE2"
cfg2.defaultCfg()
cfg2.points = 2000
cfg2.save_permutations = True
cfg2.transfer_function_method = MEM.TFMethod.External
cfg2.perm_pruning = pvec
cfg2.vars_to_integrate = CvectorPSVar()
cfg2.vars_to_marginalize = CvectorPSVar()
cfg2.vars_to_marginalize.push_back(MEM.PSVar.cos_q1)
cfg2.vars_to_marginalize.push_back(MEM.PSVar.phi_q1)
cfg2.vars_to_marginalize.push_back(MEM.PSVar.cos_qbar1)
cfg2.vars_to_marginalize.push_back(MEM.PSVar.phi_qbar1)
cfg2.do_add_jet = lambda pt, btag: btag==True
cfg2.do_calculate = lambda njet, nlep: njet>=4 and nlep==1
#cfg2.int_code += MEM.IntegrandType.AdditionalRadiation

confs = [cfg1, cfg2]
#t1 = ROOT.TF1("fb","[0]*exp(-0.5*((x-[1])/[2])**2)",0,500)
#t1.SetParameter(0, 1.0 / 100/math.sqrt(2*3.1415)) #normalization
#t1.SetParameter(1, 100) #mean
#t1.SetParameter(2, math.sqrt(2)*100) #unc

import TTH.MEAnalysis.TFClasses as TFClasses
sys.modules["TFClasses"] = TFClasses
pi_file = open(os.environ["CMSSW_BASE"]+"/src/TTH/MEAnalysis/data/transfer_functions.pickle" , 'rb')
tf_matrix = pickle.load(pi_file)
pi_file.close()

def configure_transfer_function(cfg):
    for nb in [0, 1]:
        for fl1, fl2 in [('b', MEM.TFType.bLost), ('l', MEM.TFType.qLost)]:
            tf = tf_matrix[fl1][nb].Make_CDF()
            #set pt cut for efficiency function
            tf.SetParameter(0, 30)
            tf.SetNpx(10000)
            tf.SetRange(0, 500)
            cfg.set_tf_global(fl2, nb, tf)

map(configure_transfer_function, confs)
def attach_jet_transfer_function(pt, eta):
    """
    Attaches transfer functions to the supplied jet based on the jet eta bin.
    """
    jet_eta_bin = 0
    if abs(eta)>1.0:
        jet_eta_bin = 1
    tf_b = tf_matrix['b'][jet_eta_bin].Make_Formula(False)
    tf_l = tf_matrix['l'][jet_eta_bin].Make_Formula(False)
    
    tf_b.SetNpx(10000)
    tf_b.SetRange(0, 500)
    tf_l.SetNpx(10000)
    tf_l.SetRange(0, 500)
    
    return tf_b, tf_l

def add_obj(mem, typ, **kwargs):

    if kwargs.has_key("p4s"):
        pt, eta, phi, mass = kwargs.pop("p4s")
        v = TLorentzVector()
        v.SetPtEtaPhiM(pt, eta, phi, mass);
    elif kwargs.has_key("p4c"):
        v = TLorentzVector(*kwargs.pop("p4c"))
    obsdict = kwargs.pop("obsdict", {})

    o = MEM.Object(v, typ)
    if typ == MEM.ObjectType.Jet:
        tb, tl = attach_jet_transfer_function(v.Pt(), v.Eta()) 
        o.addTransferFunction(MEM.TFType.qReco, tl)
        o.addTransferFunction(MEM.TFType.bReco, tb)
    
    for k, v in obsdict.items():
        o.addObs(k, v)
    mem.push_back_object(o)

events = map(json.loads, inlines)

for jsev in events:
    print "---- EVENT"
    results_d = {}
    for cfg in confs:
        mem = MEM.Integrand(
            MEM.output,
            cfg
        )
        jets_p4 = jsev["input"]["selectedJetsP4"]
        jets_csv = jsev["input"]["selectedJetsCSV"]
        jets_btag = jsev["input"]["selectedJetsBTag"]

        njet = 0
        nlep = 0
        for p4, btag in zip(jets_p4, jets_btag):
            print p4, btag
            if cfg.do_add_jet(p4, btag):
                add_obj(mem,
                    MEM.ObjectType.Jet,
                    p4s=p4,
                    obsdict={MEM.Observable.BTAG: btag},
                )
                njet += 1
        
        leps_p4 = jsev["input"]["selectedLeptonsP4"]
        leps_charge = jsev["input"]["selectedLeptonsCharge"]
        for p4, charge in zip(leps_p4, leps_charge):
            nlep += 1
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
        if cfg.do_calculate(njet, nlep):
            r1 = mem.run(MEM.FinalState.LH, MEM.Hypothesis.TTH, cfg.vars_to_integrate, cfg.vars_to_marginalize, cfg.points)
            r2 = mem.run(MEM.FinalState.LH, MEM.Hypothesis.TTBB, cfg.vars_to_integrate, cfg.vars_to_marginalize, cfg.points)
            results_d[cfg.name] = {
                "tth": r1.p,
                "ttbb": r2.p,
                "tth_err": r1.p_err,
                "ttbb_err": r2.p_err,
                "tth_time": r1.time/1000.0, 
                "ttbb_time": r1.time/1000.0, 
                "p": r1.p/(r1.p + 0.1*r2.p) if (r1.p>0 and r2.p>0) else 0.0 
            }
        results_d.update({
            "nMatch_wq_btag": jsev["input"]["nMatch_wq_btag"],
            "nMatch_tb_btag": jsev["input"]["nMatch_tb_btag"],
            "nMatch_hb_btag": jsev["input"]["nMatch_hb_btag"],
        })
        mem.next_event()
        del mem
    print json.dumps(results_d)
