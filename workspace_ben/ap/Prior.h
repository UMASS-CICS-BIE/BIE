// This is really -*- C++ -*-

#ifndef Prior_h
#define Prior_h

#include <string>

#include <BIEconfig.h>
#include <Distribution.h>
#include <UniformDist.h>
#include <Dirichlet.h>
#include <NormalDist.h>
#include <WeibullDist.h>
#include <InverseGammaDist.h>
#include <TransitionProb.h>
#include <CauchyDist.h>
#include <Chain.h>
#include <MHWidth.h>
#include <Ensemble.h>

#include "Serializable.h"


namespace BIE {
  
  //! @ingroup distribution
  //! @defgroup prior Prior distributions
  //! Classes that define a prior distribution for a Bayesian computation
  //! @{

  //+ CLICLASS Prior SUPER Distribution
  /**
   Combine provided input distributions into a multiple component
   Bayesian prior.

   The constructor uses a vector of distributions, one for each
   component of the parameter vector.

   In addition, InitialState member is provided to provide initial
   states to the MCMC algorithms using an Ensemble

   To provide both initial states and the prior probability from an
   Ensemble, use the PostPrior class.
  */
  class Prior : public Distribution 
  {

  protected:

    //! Ensemble for initialization (if defined)
    Ensemble* _initialdist;

    //! Hold the distribution for variates in each dimension
    vector<Distribution*> _dist;

    //! Beginning index in the state
    vector<unsigned> _pos;

    //! Dimension of each distribution
    vector<unsigned> _len;

    //! StateInfo for temporary states
    vector< boost::shared_ptr<StateInfo> > _tsi;
    //@}

    //! Temporary states
    vector<State> _sta;

    //! One transition probabilty instance for each dimension
    vector<TransProbPtr> _trans;

    //! StateInfo
    StateInfo *_si;

    //! Temporary state storage
    vector<double> v1;

    //! Note whether or not the distribution vector is locally created
    bool local_dist;

    //! Assign the transition probabilities
    virtual void initialize_proposal();

  public:
    
    /**@name Global variables */
    //@{
    /** Width factor for Metropolis-Hastings step using PostMixturePrior

	The input width is assumed to be a variance and the sample
	width is chosen following Gelman et al. 1996, ideal if the
	underlying posterior were Gaussian (see Efficient Metropolis
	jumping rules, Bayesian Statistics V, eds. Bernado, Berger,
	David & Smith, Oxford, pp. 599-608)
    */
    static double width_factor;

    //! Maximum number of iterations for sampling prior
    static unsigned max_iter;
    //@}

    /**@name Constructors */
    //@{

    //+ CLICONSTR
    //! For cloning only: makes an uninitialized prior
    Prior() : _initialdist(0), _si(0) {};

    //! Copy constructor
    Prior(Prior* p);

    //+ CLICONSTR StateInfo* clivectordist*
    //! The originating constructor (from vector [cli version])
    Prior(StateInfo* si, clivectordist *pdist);

    //! The originating constructor (from vector [stl version])
    Prior(StateInfo* si, vector<Distribution*>& dist);

    //! Prior without distributions assigned (for derived classes)
    Prior(StateInfo* si);

    //! Destructor
    ~Prior();

    //! Ensemble for initializing the MCMC states
    //+ CLIMETHOD void InitialStates Ensemble*
    void InitialStates(Ensemble *e) { _initialdist = e; }

    //@}
    
    //! Access to StateInfo
    StateInfo* const SI() { return _si; }

    //! Return the dimension
    virtual unsigned Dim() { return _si->Ndim; }

    /**@name Distribution members */
    //@{
    //! Object factor (clone)
    Prior* New();
    
    //! Differential distribution function P(x)
    virtual double PDF(State&);

    //! Log of differential distribution function P(x)
    virtual double logPDF(State&);

    //! Lower bound on distribution (in each dimension)
    virtual vector<double> lower(void);

    //! Upper bound on distribution (in each dimension)
    virtual vector<double> upper(void);
    
    //! Return mean of distribution (mulitvariate)
    virtual vector<double> Mean(void);

    //! Return standard deviation of distribution (mulitvariate)
    virtual vector<double> StdDev(void);

    //! Return specifided moment of distribution (mulitvariate)
    virtual vector<double> Moments(unsigned);
    
    //! Return random variate from distribution
    virtual State Sample(void);

    //! Return random variate from distribution for a mixture subspace
    virtual State Sample(unsigned M) { return Sample(); }

    //! Return random variate from distribution
    virtual void SamplePrior(State*);

    //! Based on prior type, return random variates
    virtual void SampleProposal(Chain &ch, MHWidth* width);
    
    //! Enforce bounds for chain state
    virtual void EnforceBounds(Chain &ch) { }

    //! Enforce bounds for state
    virtual void EnforceBoundsState(State &s) { }

    //@}
    
    #ifndef SWIG
    // AUTO GENERATED BY ../persistence/autopersist.py
    private:
    friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive & ar, const unsigned int version) {
        this->pre_serialize(ar, version);
         try {                                                         
          ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(Distribution);            
          BIE_CATCH_BOOST_SERIALIZATION_EXCEPTION;                     
         }                                                             
         try {                                                         
          ar & BOOST_SERIALIZATION_NVP(_initialdist);                        
          BIE_CATCH_BOOST_SERIALIZATION_EXCEPTION;                     
         }                                                             
         try {                                                         
          ar & BOOST_SERIALIZATION_NVP(_dist);                        
          BIE_CATCH_BOOST_SERIALIZATION_EXCEPTION;                     
         }                                                             
         try {                                                         
          ar & BOOST_SERIALIZATION_NVP(_pos);                        
          BIE_CATCH_BOOST_SERIALIZATION_EXCEPTION;                     
         }                                                             
         try {                                                         
          ar & BOOST_SERIALIZATION_NVP(_len);                        
          BIE_CATCH_BOOST_SERIALIZATION_EXCEPTION;                     
         }                                                             
         try {                                                         
          ar & BOOST_SERIALIZATION_NVP(_tsi);                        
          BIE_CATCH_BOOST_SERIALIZATION_EXCEPTION;                     
         }                                                             
         try {                                                         
          ar & BOOST_SERIALIZATION_NVP(_sta);                        
          BIE_CATCH_BOOST_SERIALIZATION_EXCEPTION;                     
         }                                                             
         try {                                                         
          ar & BOOST_SERIALIZATION_NVP(_trans);                        
          BIE_CATCH_BOOST_SERIALIZATION_EXCEPTION;                     
         }                                                             
         try {                                                         
          ar & BOOST_SERIALIZATION_NVP(_si);                        
          BIE_CATCH_BOOST_SERIALIZATION_EXCEPTION;                     
         }                                                             
         try {                                                         
          ar & BOOST_SERIALIZATION_NVP(local_dist);                        
          BIE_CATCH_BOOST_SERIALIZATION_EXCEPTION;                     
         }                                                             
        this->post_serialize(ar, version);
    }
    #endif

  };

  //! @}
  //! @}

} // namespace BIE

#ifndef SWIG
BIE_CLASS_TYPE_INFO(BIE::Prior)
BIE_CLASS_EXPORT_KEY(BIE::Prior)
#endif
#endif
