#ifndef INTEGRAND_H
#define INTEGRAND_H

// user headers
#include "TTH/MEIntegratorStandalone/interface/Parameters.h"

// Interface
extern "C" {
void ol_setparameter_int(const char *param, int val);
void ol_setparameter_double(const char *param, double val);
void ol_setparameter_string(const char *param, const char *val);
int ol_register_process(const char *process, int amptype);
int ol_n_external(int id);
void ol_phase_space_point(int id, double sqrt_s, double *pp);
void ol_start();
void ol_finish();
void ol_evaluate_tree(int id, double *pp, double *m2_tree);
void ol_evaluate_loop(int id, double *pp, double *m2_tree, double *m2_loop,
                      double *acc);
}

namespace LHAPDF {
void initPDFSet(int nset, const std::string &filename, int member = 0);
int numberPDF(int nset);
void usePDFMember(int nset, int member);
double xfx(int nset, double x, double Q, int fl);
double getXmin(int nset, int member);
double getXmax(int nset, int member);
double getQ2min(int nset, int member);
double getQ2max(int nset, int member);
void extrapolate(bool extrapolate = true);
}

namespace MEM {

class Integrand {
  /* Public interface */
 public:
  // constructor (initialise with verbosity)
  Integrand(int, const MEMConfig &);

  // detstructor
  ~Integrand();

  // add objects (jets, leptons, MET)
  // (ObjectType defined in Parameters.h)
  void push_back_object(const LV &, const ObjectType::ObjectType &);
  void push_back_object(Object *);
  //
  //  // filter out permutations
  //  // void set_permutation_strategy(const std::vector<MEM::Permutations>&);
  //  void set_permutation_strategy(
  //      const std::vector<MEM::Permutations::Permutations> &);
  //
  //  // choose what to include into the integrand
  //  void set_integrand(const int = 0);
  //
  //  // choose c.o.m. energy
  //  void set_sqrts(const double &);
  //
  //  // choose n calls
  //  void set_ncalls(const std::size_t &);
  //
  //  // set main cfg
  void set_cfg(const MEMConfig &);

  // add object information
  // WARNING: to be called just after push_back of related object!
  // (using .back() method of std::vector<>)
  void add_object_observable(const pair<Observable::Observable, double> &,
                             const ObjectType::ObjectType &);

  // main method: call it to have the ME calculated
  // void run( const FinalState =FinalState::LH, const Hypothesis
  // =Hypothesis::TTH, const std::vector<PSVar> ={});
  MEMOutput run(const FinalState::FinalState = FinalState::LH,
                const Hypothesis::Hypothesis = Hypothesis::TTH,
                const std::vector<PSVar::PSVar> = std::vector<PSVar::PSVar>(),
                const std::vector<PSVar::PSVar> = std::vector<PSVar::PSVar>(),
                int ncalls = -1);

  // clear containers and counters after each event
  void next_event();

  const TF1 *get_tf_global(TFType::TFType type, int etabin) const;

  /* Used internally */
 private:
  // initialise (once for event)
  void init(const FinalState::FinalState = FinalState::LH,
            const Hypothesis::Hypothesis = Hypothesis::TTH);

  // create a map between variable names and positions
  // void fill_map( const std::vector<PSVar>& );
  void fill_map(const std::vector<PSVar::PSVar> &);

  // make assumption
  void make_assumption(const std::vector<PSVar::PSVar> &,
                       const std::vector<PSVar::PSVar> &, MEMOutput &);

  void do_integration(unsigned int npar, double *, double *, double &, double &,
                      double &);
  //  void do_minimization(unsigned int npar, double *, double *, double &,
  //                       double &, double &);

  // clear containers before new hypothesis
  void next_hypo();

  // test if given assunption is viable
  bool test_assumption(const std::size_t &);

  // filter out permutations
  bool accept_perm(const std::vector<int> &,
                   const std::vector<Permutations::Permutations> &) const;

  bool accept_perm_btagged(const vector<int> &perm) const;
  bool accept_perm_quntagged(const vector<int> &perm) const;
  bool accept_perm_qqbarsymmetry(const vector<int> &perm) const;
  bool accept_perm_bbbarsymmetry(const vector<int> &perm) const;
  bool accept_perm_qqbarbbbarsymmetry(const vector<int> &perm) const;

  //  // filter out permutations
  //  std::vector<std::size_t> get_permutations_byBTAG(const std::size_t &)
  //  const;
  //
  //  // remove permutations not passing get_permutations_byBTAG()
  //  void filter_permutations_byBTAG();
  //
  // a constanta value for each permutation
  std::map<MEM::PermConstants::PermConstants, double> get_permutation_constants(
      const std::vector<int> &) const;

  // main method. Needed by GSLMCIntegrator
  double Eval(const double *);

  /**
  Creates a physical phase space point ps from an integrable phase space point
  x, given a permutation. Maps variables from the space $x_i \in [0, 1]$ to
  $ps_i \in [E_LOW, E_HIGH]$.
  \param ps the output phase space
  \param x the input phase space
  */
  int create_PS(MEM::PS &ps, const double *x,
                const std::vector<int> &perm) const;
  int create_PS_LH(MEM::PS &ps, const double *x,
                   const std::vector<int> &perm) const;
  int create_PS_LL(MEM::PS &ps, const double *x,
                   const std::vector<int> &perm) const;
  int create_PS_HH(MEM::PS &ps, const double *x,
                   const std::vector<int> &perm) const;
  int create_PS_TTH(MEM::PS &ps, const double *x,
                    const std::vector<int> &perm) const;

  void extend_PS(PS &, const PSPart::PSPart &, const double &, const double &,
                 const TVector3 &, size_t, const PSVar::PSVar &,
                 const PSVar::PSVar &, const PSVar::PSVar &,
                 const TFType::TFType &, int) const;
  void extend_PS_nodebug(PS &, const PSPart::PSPart &, const double &,
                         const double &, const TVector3 &) const;

  /**
  Given a VEGAS phase space point x, creates a MEM phase space point ps,
  evaluates the MEM probablity of the given permutation n_perm. The probability
  is a product of the integrand constants, the transfer functions and the
  scattering matrix:
    ps <- create_PS(x, n_perm)
    p <- constants * transfer(ps) * matrix(ps)

  /param x VEGAS phase space point
  /param n_perm index of permutation
  /return The ME*TF probability corresponding to the current hypothesis
  */
  double probability(const double *x, const std::size_t n_perm);

  double constants() const;

  double t_decay_amplitude(const LV &, const LV &, const LV &,
                           const int &) const;

  double H_decay_amplitude(const LV &, const LV &) const;

  double pdf(const double &, const double &, const double &) const;

  /**
  Evaluates the scattering amplitude (M^2) of a ttbb (+jet) system under the
  current hypothesis. It is possible to specify additional radiation in the form
  of a jet.

  /param top the top quark 4-momentum
  /param atop the antitop quark 4-momentum
  /param b1 the first b-quark 4-momentum
  /param b2 the second b-quark 4-momentum
  /param additional_jet the 4-momentum of additional radiation
  /param x1 the momentum fraction x1, may be modified
  /param x2 the momentum fraction x2, may be modified
  /return The squared amplitude of the current diagram
  */
  double scattering(const LV &top, const LV &atop, const LV &b1, const LV &b2,
                    const LV &additional_jet, double &x1, double &x2) const;

  // evaluate TF
  double transfer(const PS &, const vector<int> &, int &) const;

  // evaluate ME
  double matrix(const PS &) const;

  // evaluate ME
  double matrix_nodecay(const PS &) const;

  // solve for energy given the masses and angles
  double solve(const LV &, const double &, const double &, const TVector3 &,
               const double &, int &) const;

  // get integration edges
  void get_edges(double *, const std::vector<PSVar::PSVar> &, size_t, size_t);

  // get widths
  double get_width(const double *, const double *, const std::size_t);

  // get permutation number n
  vector<int> get_permutation(size_t);

  // setup the minimzer
  void setup_minimizer();

  // improve minimization
  void refine_minimization(std::size_t &, const ROOT::Math::Functor &,
                           const std::size_t &, double *, double *);

  // smear MET
  void smear_met();

  // smear jets
  void smear_jet(MEM::Object *, const bool &);

  // report an error
  int error_code;

  // count function calls
  int n_calls;

  // count maximum function calls
  int n_max_calls;

  // count number of invalid phase space points
  int n_skip;

  // debug
  int debug_code;

  // OpenLoops process ID-s
  map<Process::Process, int> processes;

  // main configurator
  MEMConfig cfg;

  // integration type
  Hypothesis::Hypothesis hypo;

  // number of unknowns
  std::size_t num_of_vars;

  // number of original dimensions
  std::size_t ps_dim;

  // keep track of howm many jets one would expect
  std::size_t naive_jet_counting;
  std::size_t extra_jets;

  // measured objects
  std::vector<MEM::Object *> obs_jets;
  std::vector<MEM::Object *> obs_leptons;
  std::vector<MEM::Object *> obs_mets;

  // final state
  FinalState::FinalState fs;

  // the comparateor between permutations
  CompPerm comparator;

  // contain indexes of obs_jets that need permutations
  // filled in Integrand::init()
  std::vector<int> perm_index;

  std::size_t n_perm_max;

  // the map of permutation index -> actual permutation vector
  std::vector<std::vector<int>> perm_indexes_assumption;
  // the index of the current permutation.
  // needed to keep track of the permutation inside Eval() in case of doing
  // the sum over permutations outside the integral in perm_int = 1.
  size_t this_perm;

  // counts how many times the TF was 0
  size_t tf_zero;

  std::vector<double> perm_const_assumption;
  std::vector<double> perm_btag_assumption;
  std::vector<double> perm_tmpval_assumption;
  std::vector<std::size_t> perm_pruned;

  std::vector<double> perm_btag_bb_assumption;
  std::vector<double> perm_btag_cc_assumption;
  std::vector<double> perm_btag_jj_assumption;

  vector<unsigned long> permutations;
  vector<vector<double>> permutation_probas;
  vector<vector<double>> permutation_probas_constants;
  vector<vector<double>> permutation_probas_transfer;
  vector<vector<double>> permutation_probas_me;

  // status of the prefit (0= not run, 1= run and succesfull, -1= run and
  // unsuccesfull)
  int prefit_code;

  // the btag likelihood (summed over permutations and for the maximal
  // permutation only)
  double btag_weights[3];

  // map between parameter names (physical) and positions in
  // VEGAS space
  std::unordered_map<PSVar::PSVar, std::size_t, PSVarHash, PSVarEqual>
      map_to_var;

  // map between a particle and the jet position in obs_jets/obs_leptons
  // to which the particle is matched
  std::unordered_map<PSPart::PSPart, std::size_t, PSPartHash, PSPartEqual>
      map_to_part;

  // VEGAS integrator
  ROOT::Math::GSLMCIntegrator *ig2;

  // MINUIT2 minimizer
  ROOT::Math::Minimizer *minimizer;

  // an intenal step to perform the pre-fit
  std::size_t prefit_step;

  // Stores the global cumulative transfer functions for jet
  // reconstruction
  // efficiency
  const std::map<std::pair<TFType::TFType, int>, TF1 *> tf_map;

  // the PDF-s for the b-tagger discriminants
  std::map<DistributionType::DistributionType, TH3D *> btag_pdfs;
};
}

#endif
