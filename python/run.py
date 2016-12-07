import sys, pickle, json, os

#One cool hack to make TFClasses visible for pickle
import TTH.MEIntegratorStandalone.TFClasses as TFClasses
sys.modules["TFClasses"] = TFClasses

import ROOT
ROOT.gSystem.Load("libFWCoreFWLite")
ROOT.gSystem.Load("libTTHMEIntegratorStandalone")
from ROOT import MEM
from ROOT import TLorentzVector

CvectorPermutations = getattr(ROOT, "std::vector<MEM::Permutations::Permutations>")
CvectorPSVar = getattr(ROOT, "std::vector<MEM::PSVar::PSVar>")

import ConfigParser

INTEGRATION_VARIABLES = {
    "sl_2w2h2t": [],
    "sl_1w2h2t": [
        ROOT.MEM.PSVar.cos_q1, ROOT.MEM.PSVar.phi_q1
    ],
    "sl_0w2h2t": [
        ROOT.MEM.PSVar.cos_q1, ROOT.MEM.PSVar.phi_q1, ROOT.MEM.PSVar.cos_qbar1, ROOT.MEM.PSVar.phi_qbar1
    ]
}

def add_tf(jet_eta, obj, tf_matrix):
    jet_eta_bin = 0
    if abs(jet_eta)>1.0:
        jet_eta_bin = 1
    tf_b = tf_matrix['b'][jet_eta_bin].Make_Formula(False)
    tf_l = tf_matrix['l'][jet_eta_bin].Make_Formula(False)

    tf_b.SetNpx(10000)
    tf_b.SetRange(0,400)

    tf_l.SetNpx(10000)
    tf_l.SetRange(0,400)

    obj.addTransferFunction(MEM.TFType.bReco, tf_b)
    obj.addTransferFunction(MEM.TFType.qReco, tf_l)

def parse_jet(line, tf_matrix):
    pt, eta, phi, m, btag = line.strip().split()

    v = TLorentzVector()
    v.SetPtEtaPhiM(float(pt), float(eta), float(phi), float(m))
    o = MEM.Object(v, MEM.ObjectType.Jet)
    o.addObs(MEM.Observable.BTAG, float(btag))
    add_tf(float(eta), o, tf_matrix)
    return o

def parse_lep(line):
    pt, eta, phi, m, charge = line.strip().split()

    v = TLorentzVector()
    v.SetPtEtaPhiM(float(pt), float(eta), float(phi), float(m))
    o = MEM.Object(v, MEM.ObjectType.Lepton)
    o.addObs(MEM.Observable.CHARGE, float(charge))
    return o

def parse_met(line):
    pt, phi = line.strip().split()

    v = TLorentzVector()
    v.SetPtEtaPhiM(float(pt), 0, float(phi), 0)
    o = MEM.Object(v, MEM.ObjectType.MET)
    return o

def set_tf_global(cfg, tf_matrix):
    for nb in [0, 1]:
        for fl1, fl2 in [('b', MEM.TFType.bLost), ('l', MEM.TFType.qLost)]:
            tf = tf_matrix[fl1][nb].Make_CDF()
            #set pt cut for efficiency function
            tf.SetParameter(0, 30)
            tf.SetNpx(10000)
            tf.SetRange(0,400)
            cfg.set_tf_global(fl2, nb, tf)

def integrate_config_file(conf_fn, tf_matrix):
    config = ConfigParser.ConfigParser()
    config.read(conf_fn)

    objs = []
    jets = config.get("general", "jets").split("\n")
    for jet in jets:
        if len(jet) > 0:
            jet = parse_jet(jet, tf_matrix)
            objs += [jet]

    leps = config.get("general", "leptons").split("\n")
    for lep in leps:
        if len(lep) > 0:
            lep = parse_lep(lep)
            objs += [lep]

    met = config.get("general", "met")
    if len(met) > 0:
        met = parse_met(met)
        objs += [met]

    runs = config.get("general", "runs").split()
    for run in runs:

        cfg = MEM.MEMConfig()
        cfg.defaultCfg()

        if not tf_matrix is None:
            set_tf_global(cfg, tf_matrix)

        tf_method = config.get(run, "tf_method")
        tf_method = getattr(ROOT.MEM.TFMethod, tf_method)
        cfg.transfer_function_method = tf_method

        hypothesis = config.get(run, "hypothesis")
        hypothesis = getattr(ROOT.MEM.Hypothesis, hypothesis)

        fstate = config.get(run, "fstate")
        fstate = getattr(ROOT.MEM.FinalState, fstate)

        integration_hypo = config.get(run, "integration_hypo")
        integration_variables = INTEGRATION_VARIABLES[integration_hypo]


        integrator_type = config.get(run, "integrator_type")
        integrator_type = getattr(ROOT.MEM.IntegratorType, integrator_type)
        cfg.integrator_type = integrator_type
        
        vars_to_integrate   = CvectorPSVar()
        vars_to_marginalize = CvectorPSVar()
        for v in integration_variables:
            vars_to_integrate.push_back(v)

        integration_codes = config.get(run, "integration_codes").split()
        integration_code = None
        if len(integration_codes)>0:
            integration_code = sum([getattr(ROOT.MEM.IntegrandType, ic) for ic in integration_codes])
            cfg.int_code = integration_code

        mem = MEM.Integrand(0, cfg)
        for obj in objs:
            mem.push_back_object(obj)
        ret = mem.run(fstate, hypothesis, vars_to_integrate, vars_to_marginalize)
        
        if config.has_option(run, "assert_result_within"):
            res, res_tol = config.get(run, "assert_result_within").split()
            res = float(res)
            res_tol = float(res_tol)
            assert(abs(ret.p - res)/res < res_tol)

        if config.has_option(run, "assert_time_less"):
            max_time = float(config.get(run, "assert_time_less"))
            assert(ret.time < max_time)

        if config.has_option(run, "assert_num_perm_equal"):
            num_perm = int(config.get(run, "assert_num_perm_equal"))
            assert(ret.num_perm == num_perm)

        if config.has_option(run, "assert_rel_error_less"):
            rel_err = float(config.get(run, "assert_rel_error_less"))
            assert(ret.p_err/ret.p < rel_err)

        print json.dumps({
            "conf": conf_fn,
            "run": run,
            "p": ret.p,
            "p_err": ret.p_err,
            "time": ret.time,
            "num_perm": ret.num_perm
        })

if __name__ == "__main__":

    try:
        pi_file = open(os.environ["CMSSW_BASE"] + "/src/TTH/MEIntegratorStandalone/data/transfers.pickle", 'rb')
        tf_matrix = pickle.load(pi_file)
    except Exception as e:
        print "Could not load transfer functions", e
        tf_matrix = None


    integrate_config_file(sys.argv[1], tf_matrix)

# for inf in infiles:
#     print inf
#     inf = open(inf, "r").readlines()

#     for line in inf:
#         line = line.strip()
#         if line.startswith("#"):
#             continue
#         elif line.startswith("fstate"):
#             fstate = int(line.split()[1])
#         elif line.startswith("integ"):
#             integ = int(line.split()[1])
#         elif line.startswith("hypo"):
#             hypo = int(line.split()[1])
#         elif line.startswith("tf_method"):
#             tf_method = int(line.split()[1])
#             cfg.transfer_function_method = tf_method
#         elif line.startswith("int_code"):
#             print "cfg.int_code", cfg.int_code
#             int_code = int(line.split()[1])
#             cfg.int_code = int_code
#         elif line.startswith("do_minimize"):
#             print "cfg.", cfg.do_minimize
#             do_minimize = int(line.split()[1])
#             cfg.do_minimize = do_minimize
#         elif line.startswith("memBQuark") or line.startswith("lq"):
#             jet_type, pt, eta, phi, m, btagflag, match, match_index = line.split()
#             pt, eta, phi, m, btagflag = tuple(map(float, [pt, eta, phi, m, btagflag]))
#             v = TLorentzVector()
#             v.SetPtEtaPhiM(pt, eta, phi, m)
#             o = MEM.Object(v, MEM.ObjectType.Jet)
#             o.addObs(MEM.Observable.BTAG, btagflag)
#             add_tf(eta, o)
#             mem.push_back_object(o)
#             #print "jet", pt, eta, phi, m, btagflag
#         elif line.startswith("memLepton"):
#             lep_type, pt, eta, phi, m, charge = line.split()
#             pt, eta, phi, m, charge = tuple(map(float, [pt, eta, phi, m, charge]))
#             v = TLorentzVector()
#             v.SetPtEtaPhiM(pt, eta, phi, m)
#             o = MEM.Object(v, MEM.ObjectType.Lepton)
#             o.addObs(MEM.Observable.CHARGE, charge)
#             add_tf(eta, o)
#             mem.push_back_object(o)
#             #print "lep", pt, eta, phi
#         elif line.startswith("memMET"):
#             met_type, pt, phi = line.split()
#             pt, phi = tuple(map(float, [pt, phi]))
#             v = TLorentzVector()
#             v.SetPtEtaPhiM(pt, 0.0, phi, 0.0)
#             o = MEM.Object(v, MEM.ObjectType.MET)
#             mem.push_back_object(o)
#             #print "met", pt, eta, phi

#     pvec = CvectorPermutations()
#     pvec.push_back(MEM.Permutations.BTagged)
#     pvec.push_back(MEM.Permutations.QUntagged)
#     mem.set_permutation_strategy(pvec)

#     psvar_vec = CvectorPSVar()
#     if integ == 1:
#         psvar_vec.push_back(MEM.PSVar.cos_qbar1)
#         psvar_vec.push_back(MEM.PSVar.phi_qbar1)
#     mem.set_cfg(cfg)

#     print "Calling MEM::run", fstate, hypo
#     r = mem.run(fstate, hypo, psvar_vec)
#     mem.next_event()
#     print r.p
