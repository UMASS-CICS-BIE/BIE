// This is really -*- C++ -*-

#ifndef MixturePrior_h
#define MixturePrior_h

#include <vector>
#include <string>

#include <Prior.h>
#include <TransitionProb.h>
#include <MHWidth.h>
#include <Dirichlet.h>
#include <Chain.h>
#include <Poisson.h>

#include "Serializable.h"


namespace BIE {
  
  //+ CLICLASS MixturePrior SUPER Prior
  //! Abstract class.  A prior for a mixture model.
  //! 
  //! The prior for the mixture parameters will be the Dirichlet 
  //! distribution.  See the Prior class for a definition of 
  //! distribution types for each element of the model parameter vector.
  class MixturePrior : public Prior 
  {

  protected:

    //! Return state
    State ret;
    
    //! Dirichlet distribution parameter
    double Alpha;

    //! One distribution for each dimension
    vector<boost::shared_ptr<Dirichlet> > mix;

    //! Uniform random variates in [-1.0, 1.0]
    boost::shared_ptr<Uniform> unit;

    //! Uniform random variates in [0.0, 1.0]
    boost::shared_ptr<Uniform> useg;

    //! Normal with zero mean and unit variance
    boost::shared_ptr<Normal> normal;

    //! Poisson distribution mean
    double pmean;

    //! Poisson distribution
    boost::shared_ptr<Poisson> pois;

    //! Gamma-distributed variates
    double GenGamma(double a);

    //! Beta-distributed variates
    double GenBeta(double a, double b);

    //! Create the M-H proposal distributions
    virtual void initialize_proposal();

    //! Vector of mixture component prior distributions
    vector<Distribution*> _prior0;

    //@{
    //! Beginning position of each distribution
    vector<unsigned> _pos0, _pos1;
    //@}

    //@{
    //! Length of each distribution
    vector<unsigned> _len0, _len1;
    //@}

    //@{
    //! StateInfo for temporary states
    vector< boost::shared_ptr<StateInfo> > _tsi0, _tsi1;
    //@}

    //@{
    //! Temporary states
    vector<State> _sta0, _sta1;
    //@}

    /** Vector of mixture component prior distributions for the
	extended component */
    vector<Distribution*> _prior1;

    //! Sample the mixture model components
    vector<double> priorSampleMix();

    //! Sample the extended components
    vector<double> priorSampleExt();

  public:

    //! Maximum number of tries to find a good state
    static unsigned ITMAX;
    
    //! Very small but non-zero probability "epsilon"
    static double min_prob;

    /**@name Constructors */
    //@{

    //! For cloning only: makes an uninitialized prior
    MixturePrior();

    //! Copy constructor for derived classes
    MixturePrior(MixturePrior*);

    /** The useful constructor
	@param si    is the StateInfo definition
	@param vdist is a vector of distributions for each parameter
    */
    MixturePrior(StateInfo *si, clivectordist* vdist);

    /** The even more useful constructor
	@param si    is the StateInfo definition
	@param alpha is the Dirichlet distribution parameter
	@param pfile is defintion of the prior distribution for parameter vector
    */
    MixturePrior(StateInfo *si, double alpha, string pfile);

    /** The even more useful constructor
	@param si     is the StateInfo definition
	@param alpha  is the Dirichlet distribution parameter
	@param vdist0 is a vector of distributions for each parameter
	@param vdist1 is a vector of distributions for the extended space
    */
    MixturePrior(StateInfo* si, double alpha, clivectordist* vdist0, 
		 clivectordist* vdist1=0);

    /** Constructor for a variable number of components
	@param si    is the StateInfo definition
	@param alpha is the Dirichlet distribution parameter
	@param pmean is Poisson mean for the prior distribution of components
	@param pfile is defintion of the prior distribution for parameter vector
    */
    MixturePrior(StateInfo *si, double alpha, double pmean, string pfile);

    /** Constructor for a variable number of components
	@param si    is the StateInfo definition
	@param alpha is the Dirichlet distribution parameter
	@param pmean is Poisson mean for the prior distribution of components
	@param vdist0 is a vector of distributions for each parameter
	@param vdist1 is a vector of distributions for the extended space
    */
    MixturePrior(StateInfo *si, double alpha, double pmean, 
		 clivectordist* vdist0, clivectordist* vdist1=0);

    //! Destructor
    ~MixturePrior();

    //@}

    //+ CLIMETHOD void setMaxIter int
    /** 
	Set maximum number of iterations in state sampler to generate
	a state within the desired bounds
    */
    void setMaxIter(int i) { ITMAX = i; }

    /**@name Distribution members
       The packing order in vectors is:
       Nmix weights, followed by Nmix vectors of Ndim elements
    */
    //@{
    //! Object factor (clone)
    virtual MixturePrior* New() = 0;
    
    //! Differential distribution function P(x)
    virtual double PDF(State&) = 0;

    //! Log of differential distribution function P(x)
    virtual double logPDF(State&) = 0;

    /** Log of differential distribution function P(x) for subspace @param M
	and component @param n for parameter vector @param V

	This method is used for ReversibleJump classes.
    */
    virtual double logPDFMarginal(unsigned M, unsigned n, 
				  const vector<double>& V) = 0;

    //! Lower bound on distribution
    virtual vector<double> lower(void) = 0;

    //! Upper bound on distribution
    virtual vector<double> upper(void) = 0;
    
    /*
    //! Return mean of distribution
    virtual vector<double> Mean(unsigned m) = 0;

    //! Return standard deviation of distribution
    virtual vector<double> StdDev(unsigned m) = 0;

    //! Return specifided moment of distribution
    virtual vector<double> Moments(unsigned, unsigned) = 0;
    */
    
    //! Sample the number of components
    virtual unsigned SampleM(void) {
      unsigned m = _si->M;
      if (pois) {
	do {
	  m = static_cast<unsigned>(floor((*pois)()+1.0e-6));
	} while (m>_si->M || m<1);
      }
      return m;
    }

    //! Return random variate for variable dimension model
    virtual State Sample(void) { return Sample(SampleM()); }

    //! Return random variate from distribution
    virtual State Sample(unsigned m) = 0;

    //@}

    //! Sample default
    virtual void SamplePrior(State *s);

    //! Sample: a parameter vector for a single component in a subspace
    virtual void SamplePrior(unsigned, unsigned, vector<double>&) = 0;

    //! Based on prior type, return random variates for mixture model
    virtual void SampleProposal(Chain &ch, MHWidth* width);
    
    //! Choose a new component weight conditional on old weights
    virtual double BirthWeight(unsigned M);

    //! Probability of new weight and old weights
    virtual double BirthWeightPDF(unsigned M, double p);

    //! Check existence of subspace for m
    virtual bool Exists(int m) { return m>0 && m<static_cast<int>(_si->M); }

    #ifndef SWIG
    // AUTO GENERATED BY ../persistence/autopersist.py
    private:
    friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive & ar, const unsigned int version) {
        this->pre_serialize(ar, version);
         try {                                                         
          ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(Prior);            
          BIE_CATCH_BOOST_SERIALIZATION_EXCEPTION;                     
         }                                                             
         try {                                                         
          ar & BOOST_SERIALIZATION_NVP(ret);                        
          BIE_CATCH_BOOST_SERIALIZATION_EXCEPTION;                     
         }                                                             
         try {                                                         
          ar & BOOST_SERIALIZATION_NVP(Alpha);                        
          BIE_CATCH_BOOST_SERIALIZATION_EXCEPTION;                     
         }                                                             
         try {                                                         
          ar & BOOST_SERIALIZATION_NVP(mix);                        
          BIE_CATCH_BOOST_SERIALIZATION_EXCEPTION;                     
         }                                                             
         try {                                                         
          ar & BOOST_SERIALIZATION_NVP(unit);                        
          BIE_CATCH_BOOST_SERIALIZATION_EXCEPTION;                     
         }                                                             
         try {                                                         
          ar & BOOST_SERIALIZATION_NVP(useg);                        
          BIE_CATCH_BOOST_SERIALIZATION_EXCEPTION;                     
         }                                                             
         try {                                                         
          ar & BOOST_SERIALIZATION_NVP(normal);                        
          BIE_CATCH_BOOST_SERIALIZATION_EXCEPTION;                     
         }                                                             
         try {                                                         
          ar & BOOST_SERIALIZATION_NVP(pmean);                        
          BIE_CATCH_BOOST_SERIALIZATION_EXCEPTION;                     
         }                                                             
         try {                                                         
          ar & BOOST_SERIALIZATION_NVP(pois);                        
          BIE_CATCH_BOOST_SERIALIZATION_EXCEPTION;                     
         }                                                             
         try {                                                         
          ar & BOOST_SERIALIZATION_NVP(_prior0);                        
          BIE_CATCH_BOOST_SERIALIZATION_EXCEPTION;                     
         }                                                             
         try {                                                         
          ar & BOOST_SERIALIZATION_NVP(_pos0);                        
          BIE_CATCH_BOOST_SERIALIZATION_EXCEPTION;                     
         }                                                             
         try {                                                         
          ar & BOOST_SERIALIZATION_NVP(_pos1);                        
          BIE_CATCH_BOOST_SERIALIZATION_EXCEPTION;                     
         }                                                             
         try {                                                         
          ar & BOOST_SERIALIZATION_NVP(_len0);                        
          BIE_CATCH_BOOST_SERIALIZATION_EXCEPTION;                     
         }                                                             
         try {                                                         
          ar & BOOST_SERIALIZATION_NVP(_len1);                        
          BIE_CATCH_BOOST_SERIALIZATION_EXCEPTION;                     
         }                                                             
         try {                                                         
          ar & BOOST_SERIALIZATION_NVP(_tsi0);                        
          BIE_CATCH_BOOST_SERIALIZATION_EXCEPTION;                     
         }                                                             
         try {                                                         
          ar & BOOST_SERIALIZATION_NVP(_tsi1);                        
          BIE_CATCH_BOOST_SERIALIZATION_EXCEPTION;                     
         }                                                             
         try {                                                         
          ar & BOOST_SERIALIZATION_NVP(_sta0);                        
          BIE_CATCH_BOOST_SERIALIZATION_EXCEPTION;                     
         }                                                             
         try {                                                         
          ar & BOOST_SERIALIZATION_NVP(_sta1);                        
          BIE_CATCH_BOOST_SERIALIZATION_EXCEPTION;                     
         }                                                             
         try {                                                         
          ar & BOOST_SERIALIZATION_NVP(_prior1);                        
          BIE_CATCH_BOOST_SERIALIZATION_EXCEPTION;                     
         }                                                             
         try {                                                         
          ar & BOOST_SERIALIZATION_NVP(ITMAX);                        
          BIE_CATCH_BOOST_SERIALIZATION_EXCEPTION;                     
         }                                                             
         try {                                                         
          ar & BOOST_SERIALIZATION_NVP(min_prob);                        
          BIE_CATCH_BOOST_SERIALIZATION_EXCEPTION;                     
         }                                                             
        this->post_serialize(ar, version);
    }
    #endif

  };
} // namespace BIE

#ifndef SWIG
BIE_CLASS_ABSTRACT(BIE::MixturePrior)

BIE_CLASS_EXPORT_KEY(BIE::MixturePrior)
#endif
#endif
