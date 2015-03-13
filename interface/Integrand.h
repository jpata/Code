#ifndef INTEGRAND_H
#define INTEGRAND_H

// user headers
#include "interface/Parameters.h"

extern "C" {
  void pphttxcallme2born_( double*, double[20], double*, double* );
}
extern "C" {
  void ppttxbbxcallme2born_( double*, double[24], double*, double* );
}


namespace LHAPDF {
  void initPDFSet(int nset, const std::string& filename, int member=0);
  int numberPDF (int nset);
  void usePDFMember(int nset, int member);
  double xfx(int nset, double x, double Q, int fl);
  double getXmin(int nset, int member);
  double getXmax(int nset, int member);
  double getQ2min(int nset, int member);
  double getQ2max(int nset, int member);
  void extrapolate(bool extrapolate=true);
}

namespace MEM {
  
  class Integrand {

    /* Public interface */
  public:    

    // constructor (initialise with verbosity)
    Integrand( int =0);
    
    // detstructor
    ~Integrand();

    // add objects (jets, leptons, MET)
    // (ObjectType defined in Parameters.h)
    void push_back_object( const LV& , const ObjectType::ObjectType&);    
    void push_back_object( Object* );

    // filter out permutations
    //void set_permutation_strategy(const std::vector<MEM::Permutations>&);
    void set_permutation_strategy(const std::vector<MEM::Permutations::Permutations>&);

    // choose what to include into the integrand
    void set_integrand(const int =0);
    
    // choose c.o.m. energy
    void set_sqrts(const double&);

    // choose n calls
    void set_ncalls(const std::size_t&);

    // add object information
    // WARNING: to be called just after push_back of related object!
    // (using .back() method of std::vector<>)
    void add_object_observable(const pair<Observable::Observable, double>&, const ObjectType::ObjectType& );

    // main method: call it to have the ME calculated
    //void run( const FinalState =FinalState::LH, const Hypothesis =Hypothesis::TTH, const std::vector<PSVar> ={});
    void run( const FinalState::FinalState=FinalState::LH,
        const Hypothesis::Hypothesis =Hypothesis::TTH,
        const std::vector<PSVar::PSVar> = std::vector<PSVar::PSVar>()
    );

    // clear containers and counters after each event
    void next_event();

    /* Used internally */
  private:

    // initialise (once for event)
    void init( const FinalState::FinalState =FinalState::LH,
        const Hypothesis::Hypothesis=Hypothesis::TTH);

    // create a map between variable names and positions
    //void fill_map( const std::vector<PSVar>& );
    void fill_map( const std::vector<PSVar::PSVar>& );

    // make assumption
    double make_assumption( const std::vector<PSVar::PSVar>& );

    // clear containers before new hypothesis
    void next_hypo();

    // test if given assunption is viable
    bool test_assumption( const std::size_t&, std::size_t& );

    // filter out permutations
    bool accept_perm( const std::vector<int>&, const std::vector<Permutations::Permutations>& ) const;

    // a constanta value for each permutation
    double get_permutation_constants(  const std::vector<int>& ) const;

    // main method. Needed by GSLMCIntegrator
    double Eval(const double*) const;
    
    // create PS point
    int create_PS    (MEM::PS&, const double*, const std::vector<int>&) const;
    int create_PS_LH (MEM::PS&, const double*, const std::vector<int>&) const;
    int create_PS_LL (MEM::PS&, const double*, const std::vector<int>&) const;
    int create_PS_HH (MEM::PS&, const double*, const std::vector<int>&) const;
    int create_PS_TTH(MEM::PS&, const double*, const std::vector<int>&) const;

    void extend_PS(PS&, const PSPart::PSPart&, const double& , const double& , const TVector3&, const int&, const PSVar::PSVar&,const PSVar::PSVar&,const PSVar::PSVar&,const TFType::TFType&, const int =0) const;
    void extend_PS_nodebug(PS&, const PSPart::PSPart&, const double& , const double& , const TVector3&) const;

    double probability(const double*, const vector<int>&) const;

    double constants() const;

    double t_decay_amplitude(const TLorentzVector&, const TLorentzVector&, const TLorentzVector&, const int&) const;

    double H_decay_amplitude(const TLorentzVector&, const TLorentzVector&) const;

    double pdf(const double&, const double&, const double&) const;

    double scattering(const TLorentzVector&, const TLorentzVector&, const TLorentzVector&, const TLorentzVector&, 
		      double&, double&) const;

    // evaluate TF
    double transfer(const PS&, const vector<int>&) const;

    // evaluate ME
    double matrix  (const PS&) const;

    // evaluate ME
    double matrix_nodecay (const PS&) const;

    // solve for energy given the masses and angles
    double solve( const LV&,  const double& ,  const double& , const TVector3&, const double&, int&) const;

    // get integration edges
    void get_edges(double*, const std::vector<PSVar::PSVar>&, const std::size_t&, const std::size_t&);

    // get widths
    double get_width(const double*, const double*, const std::size_t);

    // report an error
    int error_code;

    // count function calls
    int n_calls;

    // set number of calls
    int n_max_calls;

    // count number of invalid phase space points
    int n_skip;

    // debug
    int debug_code;

    // what to include within the integrand
    int int_code;

    // integration type
    Hypothesis::Hypothesis hypo;

    // number of unknowns
    std::size_t num_of_vars;

    // number of original dimensions
    std::size_t ps_dim;

    // keep track of howm many jets one would expect
    std::size_t naive_jet_counting;

    // measured objects
    std::vector<MEM::Object*> obs_jets;
    std::vector<MEM::Object*> obs_leptons;
    std::vector<MEM::Object*> obs_mets;

    // final state
    FinalState::FinalState fs;

    // contain indexes of obs_jets that need permutations
    std::vector<std::vector<int> > perm_indexes;
    std::vector<std::vector<int> > perm_indexes_assumption;
    std::vector< double >      perm_const_assumption;


    // map between parameter names (physical) and positions in
    // VEGAS space
    boost::unordered_map< PSVar::PSVar, std::size_t, PSVarHash, PSVarEqual> map_to_var;

    // map between a particle and the jet position in obs_jets/obs_leptons
    // to which the particle is matched
    boost::unordered_map< PSPart::PSPart, std::size_t, PSPartHash, PSPartEqual> map_to_part;

    // strategy on permutations
    std::vector<MEM::Permutations::Permutations> permutation_strategies;

    // VEGAS integrator
    ROOT::Math::GSLMCIntegrator* ig2;

    // the com energy
    double Sqrt_s;

    //transform energies into adimensional quantitities
    double EMAX;

  };

}

#endif
