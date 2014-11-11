#include "TLorentzVector.h"
#include "TVector3.h"
#include "TFormula.h"
#include "TMath.h"

#include<cmath>
#include<cstddef>
#include<iostream>
#include<cstddef>
#include<cstdio>
#include<map> 
#include<string>
#include<algorithm>
#include<assert.h>
#include<memory> 
#include<limits> 

using namespace std;

typedef TLorentzVector LV;

#define VERBOSE 0

namespace Algo {

  constexpr double MTOP = 174.3;
  constexpr double MB   = 4.8;
  constexpr double MW   = 80.19;
  constexpr double DM2  = (MTOP*MTOP-MB*MB-MW*MW)*0.5;
  constexpr double MH   = 125.;
  constexpr double DMH2 = (MH*MH-2*MB*MB)*0.5;

  const string TF_Q = "TMath::Gaus(x, [0] + [1]*y, y*TMath::Sqrt([2]+[3]/y+[4]/y/y))";
  const string TF_B = "TMath::Gaus(x, [0] + [1]*y, y*TMath::Sqrt([2]+[3]/y+[4]/y/y))";
  const double TF_Q_param[2][5] = 
    {  { 0.0e+00, 1.0e+00, 0.0e+00, 1.5e+00, 0.0e+00 },
       { 0.0e+00, 1.0e+00, 1.3e+01, 1.5e+00, 0.0e+00 } 
    };

  const string TF_MET = "TMath::Gaus(x,0.,20)*TMath::Gaus(y,0.,20)";

  size_t eta_to_bin( const LV& );


  enum Decay { TopLep, TopHad, WHad, HiggsHad, Radiation, MET, UNKNOWN };
  string translateDecay(Decay&);


  enum FinalState { TopLep_l=0, TopLep_b, 
		    TopHad_q,   TopHad_qbar,  TopHad_b, 
		    WHad_q,     WHad_qbar,
		    HiggsHad_b, HiggsHad_bbar, 
		    Radiation_q };

  
  
  class TransferFunction  {
    
  public:
    TransferFunction(const string&, const string&);
    ~TransferFunction();
    void init(const double*);
    const string getFormula() const; 
    double eval (const double& , const double&) const ;
  private:
    string formula;
    TFormula* f;
  };
  
  
  class DecayBuilder {
  public:
    virtual ~DecayBuilder() {};
    virtual double eval( const double* , LV&) = 0; 
    virtual void print(ostream&) = 0;     
  };
  
  
  
  class CombBuilder: public DecayBuilder {
  public:
    CombBuilder();
    CombBuilder(vector<DecayBuilder*>&);
    ~CombBuilder();
    void add(DecayBuilder*);
    double eval ( const double* , LV&);
    double eval ( const double* );
    void print(ostream&);   
  private:
    vector<DecayBuilder*> combined;
  };


  class METBuilder: public DecayBuilder {
  public:
    METBuilder();
    ~METBuilder();
    void init(const LV&);
    double eval (  const double* , LV& );
    void print(ostream&);   
  private:
    LV p4_invisible;
    TransferFunction* tf_met;
    Decay decay;
  };
  
  
  class TopHadBuilder: public DecayBuilder {
    
  public:
    TopHadBuilder();
    ~TopHadBuilder();
    void init(const FinalState& , const LV&, const size_t&);
    double eval ( const double* , LV&) ;
    void print(ostream&);     
    
  private:
    
    LV p4_q;
    LV p4_qbar;
    LV p4_b;
    size_t index_q;
    size_t index_qbar;
    size_t index_b;
    TransferFunction* tf_q;
    TransferFunction* tf_qbar;
    TransferFunction* tf_b;
    Decay decay;
    size_t errFlag;
  };

  
  class WHadBuilder: public DecayBuilder {
    
  public:
    WHadBuilder();
    ~WHadBuilder();
    void init(const FinalState& , const LV&, const size_t&);
    double eval ( const double* , LV&) ;
    void print(ostream&);     
    
  private:
  
    LV p4_q;
    LV p4_qbar;
    size_t index_q;
    size_t index_qbar;
    TransferFunction* tf_q;
    TransferFunction* tf_qbar;
    Decay decay;
    size_t errFlag;
  };


 class TopLepBuilder: public DecayBuilder {
    
  public:
    TopLepBuilder();
    ~TopLepBuilder();
    void init(const FinalState& , const LV&, const size_t&);
    double eval ( const double* , LV&) ;
    void print(ostream&);     
    
  private:
    
    LV p4_l;
    LV p4_b;
    size_t index_l;
    size_t index_b;
    TransferFunction* tf_b;
    Decay decay;
    size_t errFlag;
  };



}
