#include "TTH/MEIntegratorStandalone/interface/Integrand.h"
#include "Math/Factory.h"
#include "Math/Functor.h"
#include "Math/GSLMCIntegrator.h"
#include "Math/IOptions.h"
#include "Math/IntegratorOptions.h"
#include "Math/AllIntegrationTypes.h"
#include "TFile.h"
#include <iostream>

using namespace std;
using namespace MEM;

const auto npoints = 10000;

// Returns the transfer function corresponding to a jet flavour and eta
TF1* getTransferFunction(TFile* tffile, const char* flavour, double eta) {
    int etabin = 0;
    if (std::abs(eta) > 1.0) {
        etabin = 1;
    }
    stringstream ss;
    ss << "tf_" << flavour << "_etabin" << etabin;
    const char* fname = ss.str().c_str();
    TF1* tf = (TF1*)(tffile->Get(fname));
    if (tf == 0) {
        cerr << "could not get transfer function " << fname << endl;
        cerr << flush;
        throw exception();
    }
    return tf;
}

int main(){

  //Load the transfer functions
  TFile* tffile = new TFile("root/transfer.root");


  //create a MEM configuration.
  //this needs to be done once per job, not for every event
  MEMConfig cfg;
  cfg.defaultCfg();
  cfg.transfer_function_method = TFMethod::External;
  vector<Permutations::Permutations> pvec({Permutations::BTagged, Permutations::QUntagged, Permutations::QQbarBBbarSymmetry});
  cfg.perm_pruning = pvec;
  //cfg.int_code = cfg.int_code + MEM::IntegrandType::AdditionalRadiation;
  cfg.int_code = (
    MEM::IntegrandType::IntegrandType::Constant |
    MEM::IntegrandType::IntegrandType::ScattAmpl |
    MEM::IntegrandType::IntegrandType::DecayAmpl |
    MEM::IntegrandType::IntegrandType::Jacobian |
    MEM::IntegrandType::IntegrandType::PDF |
    MEM::IntegrandType::IntegrandType::Transfer
  );
  cfg.perm_int = 0;
  cfg.integrator_type = MEM::IntegratorType::Vegas;
  //cfg.cuba_cores = 10;

  //Transfer functions for jet reconstruction efficiency
  cfg.set_tf_global(TFType::bLost, 0, getTransferFunction(tffile, "beff", 0.0));
  cfg.set_tf_global(TFType::bLost, 1, getTransferFunction(tffile, "beff", 2.0));
  cfg.set_tf_global(TFType::qLost, 0, getTransferFunction(tffile, "leff", 0.0));
  cfg.set_tf_global(TFType::qLost, 1, getTransferFunction(tffile, "leff", 2.0));

  //Create the mem integrator, once per job
  Integrand* integrand = new Integrand(0,cfg);
  
  //Add some objects to the MEM
//  TLorentzVector lv_j1;
//  lv_j1.SetPtEtaPhiM(158.72032165527344, -2.042543649673462, 1.3954641819000244, 10.526216506958008);
//  Object j1( lv_j1, ObjectType::Jet );
//  j1.addObs( Observable::BTAG, 0. ); // 0 - jet is assumed to be from a light quark, 1 - a b quark
//  j1.addTransferFunction(TFType::bReco, getTransferFunction(tffile, "b", lv_j1.Eta()));
//  j1.addTransferFunction(TFType::qReco, getTransferFunction(tffile, "l", lv_j1.Eta()));

  TLorentzVector lv_j2;
  lv_j2.SetPtEtaPhiM(154.27105712890625, 0.4824397563934326, -1.7387510538101196, 13.035625457763672);
  Object j2( lv_j2, ObjectType::Jet );
  j2.addObs( Observable::BTAG, 1. );
  j2.addTransferFunction(TFType::bReco, getTransferFunction(tffile, "b", lv_j2.Eta()));
  j2.addTransferFunction(TFType::qReco, getTransferFunction(tffile, "l", lv_j2.Eta()));

  TLorentzVector lv_j3;
  lv_j3.SetPtEtaPhiM(146.58963012695312, 0.4863491654396057, -2.5873682498931885, 16.60299301147461);
  Object j3( lv_j3, ObjectType::Jet );
  j3.addObs( Observable::BTAG, 1. );
  j3.addTransferFunction(TFType::bReco, getTransferFunction(tffile, "b", lv_j3.Eta()));
  j3.addTransferFunction(TFType::qReco, getTransferFunction(tffile, "l", lv_j3.Eta()));

  TLorentzVector lv_j4;
  lv_j4.SetPtEtaPhiM(141.0120849609375, 1.432065725326538, 2.3994698524475098, 11.62883472442627);
  Object j4( lv_j4, ObjectType::Jet );
  j4.addObs( Observable::BTAG, 0. );
  j4.addTransferFunction(TFType::bReco, getTransferFunction(tffile, "b", lv_j4.Eta()));
  j4.addTransferFunction(TFType::qReco, getTransferFunction(tffile, "l", lv_j4.Eta()));

  TLorentzVector lv_j5;
  lv_j5.SetPtEtaPhiM(135.69717407226562, 0.854595422744751, 2.2066357135772705, 10.661768913269043);
  Object j5( lv_j5, ObjectType::Jet );
  j5.addObs( Observable::BTAG, 0. );
  j5.addTransferFunction(TFType::bReco, getTransferFunction(tffile, "b", lv_j5.Eta()));
  j5.addTransferFunction(TFType::qReco, getTransferFunction(tffile, "l", lv_j5.Eta()));

  TLorentzVector lv_j6;
  lv_j6.SetPtEtaPhiM(71.66462707519531, 1.1147513389587402, -2.8738958835601807, 8.748970031738281);
  Object j6( lv_j6, ObjectType::Jet );
  j6.addObs( Observable::BTAG, 1. );
  j6.addTransferFunction(TFType::bReco, getTransferFunction(tffile, "b", lv_j6.Eta()));
  j6.addTransferFunction(TFType::qReco, getTransferFunction(tffile, "l", lv_j6.Eta()));

  TLorentzVector lv_j7;
  lv_j7.SetPtEtaPhiM(41.40342330932617, 1.1593384742736816, -0.756840705871582, 4.991734504699707);
  Object j7( lv_j7, ObjectType::Jet );
  j7.addObs( Observable::BTAG, 1. );
  j7.addTransferFunction(TFType::bReco, getTransferFunction(tffile, "b", lv_j7.Eta()));
  j7.addTransferFunction(TFType::qReco, getTransferFunction(tffile, "l", lv_j7.Eta()));

  
  //create a lepton
  TLorentzVector lv_l1;
  lv_l1.SetPtEtaPhiM(42.008724212646484, 1.87186598777771, 0.4528217911720276, -0.042224351316690445);
  Object l1( lv_l1, ObjectType::Lepton );
  l1.addObs( Observable::CHARGE, +1. );

  //create a MET
  TLorentzVector lv_met;
  lv_met.SetPtEtaPhiM(320.2880554199219, 0.0, -0.37210184335708624, 0.0);
  Object met( lv_met, ObjectType::MET );
  
  for (int i=0; i < 1; i++) {
    //add all objects to the MEM integrator
    //integrand->push_back_object( &j1 );
    integrand->push_back_object( &j2 );
    integrand->push_back_object( &j3 );
    integrand->push_back_object( &j4 );
    integrand->push_back_object( &j5 );
    integrand->push_back_object( &j6 );
    integrand->push_back_object( &j7 );
    integrand->push_back_object( &l1 );
    integrand->push_back_object( &met );

    MEMOutput res;			   
    
    //Evaluate fully reconstructed hypothesis
    //LH - single-leptonic decay channel
    //LL - dileptonic decay channel
    cout << "Fully reconstructed interpretation" << endl;
    cout << "evaluating tth hypo" << endl;
    //third variable is variables to integrate over
    //if nothing is specified, assume that all jets (4b + 2 light) were reconstructed
    res = integrand->run( FinalState::LH, Hypothesis::TTH,  {}, {}, npoints );
    cout << "p = " << res.p << " +- " << res.p_err << endl;
    double p0 = res.p;

    cout << "evaluating ttbb hypo" << endl;
    res = integrand->run( FinalState::LH, Hypothesis::TTBB, {}, {}, npoints  );
    cout << "p = " << res.p << " +- " << res.p_err << endl;
    double p1 = res.p;

    //this is the final discriminator value. the normalization constant is to be optimized
    double mem_w = p0 / (p0 + 0.02*p1);
    cout << "mem 222 discriminator " << mem_w << endl;

    //Evaluate 022 hypothesis. We do not use the information provided by the light quarks.
    cout << "Integrating over light quarks" << endl;
    cout << "evaluating tth hypo" << endl;
    //integrate over the light quark angles
    res = integrand->run( FinalState::LH, Hypothesis::TTH,  {PSVar::cos_q1, PSVar::phi_q1, PSVar::cos_qbar1, PSVar::phi_qbar1}, {}, npoints );
    cout << "p = " << res.p << " +- " << res.p_err << endl;
    p0 = res.p;
    
    cout << "evaluating ttbb hypo" << endl;
    res = integrand->run( FinalState::LH, Hypothesis::TTBB, {PSVar::cos_q1, PSVar::phi_q1, PSVar::cos_qbar1, PSVar::phi_qbar1}, {}, npoints );
    cout << "p = " << res.p << " +- " << res.p_err << endl;
    p1 = res.p;
    
    mem_w = p0 / (p0 + 0.02*p1);
    cout << "mem 022 discriminator " << mem_w << endl;

    integrand->next_event();
  }

  delete integrand;
}
