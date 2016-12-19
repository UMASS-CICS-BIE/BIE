// This is really -*- C++ -*-

#ifndef InitialMixturePrior_h
#define InitialMixturePrior_h

#include <ACG.h>
#include <Uniform.h>
#include <Normal.h>

#include <MixturePrior.h>

#include "Serializable.h"


namespace BIE {
  
  //! @addtogroup prior Prior distributions
  //! @{

  //+ CLICLASS InitialMixturePrior SUPER MixturePrior
  /**
     Prior for mixture model
  */
  class InitialMixturePrior : public MixturePrior 
  {

  protected:
    
    //! Temporary state vector
    State ttwght;

  public:
    
    /**@name Constructors */
    //@{

    //+ CLICONSTR
    //! For cloning only: makes an uninitialized prior
    InitialMixturePrior();

    //+ CLICONSTR StateInfo* double clivectordist*
    /** The useful constructor
	@param si    is the StateInfo definition
	@param alpha is the Dirichlet parameter
	@param vdist is the vector of distributions, one for each each dimension of the parameter vector
    */
    InitialMixturePrior(StateInfo *si, double alpha, clivectordist *vdist);
    
    
    //+ CLICONSTR StateInfo* double double clivectordist*
    //+ CLICONSTR StateInfo* double double clivectordist* clivectordist*
    /** The useful constructor
	@param si    is the StateInfo definition
	@param alpha is the Dirichlet parameter
	@param pmean is the Poisson prior mean for the number of components
	@param vdist0 is the vector of distributions, one for each each dimension of the parameter vector
	@param vdist1 is the vector of distributions, one for each each dimension of the extended vector
    */
    InitialMixturePrior(StateInfo *si, double alpha, double pmean, clivectordist *vdist0, clivectordist* vdist1=0);

    /**@name Distribution members
       The packing order in vectors is:
       Nmix weights, followed by Nmix vectors of Ndim elements
     */
    //@{
    //! Object factor (clone)
    virtual InitialMixturePrior* New();
    
    //! Differential distribution function P(x)
    virtual double PDF(State&);

    //! Log of differential distribution function P(x)
    virtual double logPDF(State&);

    //! Log of differential distribution function
    //! @param M is the mixture subspace
    //! @param n is the component index
    //! @param V is the input parameter vector for component n
    virtual double logPDFMarginal(unsigned M, unsigned n, 
				  const vector<double>& V);

    //! Differential distribution function for single component
    // virtual double PDF(int, State&);

    //! Lower bound on distribution
    virtual vector<double> lower(void);

    //! Upper bound on distribution
    virtual vector<double> upper(void);
    
    //! Return mean of distribution
    virtual vector<double> Mean(unsigned m);

    //! Return standard deviation of distribution
    virtual vector<double> StdDev(unsigned m);

    //! Return specifided moment of distribution
    virtual vector<double> Moments(unsigned m, unsigned n);
    
    //! Return random variate from distribution
    virtual State Sample(unsigned m);
    //@}

    //! Sample based on StateInfo
    virtual void SamplePrior(State* s);

    //! Sample: return a single component @param V for subspace @param
    //! M and component @param n
    virtual void SamplePrior(unsigned M, unsigned n, vector<double>& V);

  private:
    //! Temporary storage for mixture evaluation
    vector<double> p;

    #ifndef SWIG
    // AUTO GENERATED BY ../persistence/autopersist.py
    private:
    friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive & ar, const unsigned int version) {
        this->pre_serialize(ar, version);
         try {                                                         
          ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(MixturePrior);            
          BIE_CATCH_BOOST_SERIALIZATION_EXCEPTION;                     
         }                                                             
         try {                                                         
          ar & BOOST_SERIALIZATION_NVP(ttwght);                        
          BIE_CATCH_BOOST_SERIALIZATION_EXCEPTION;                     
         }                                                             
        this->post_serialize(ar, version);
    }
    #endif

  };

  //! @}
}

#ifndef SWIG
BIE_CLASS_TYPE_INFO(BIE::InitialMixturePrior)
BIE_CLASS_EXPORT_KEY(BIE::InitialMixturePrior)
#endif
#endif