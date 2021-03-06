// This is really -*- C++ -*-

#ifndef InitialMixturePriorPoisson_h
#define InitialMixturePriorPoisson_h

#include <ACG.h>
#include <Uniform.h>
#include <Normal.h>

#include <InitialMixturePrior.h>

#include "Serializable.h"


namespace BIE {
  
  //+ CLICLASS InitialMixturePriorPoisson SUPER InitialMixturePrior
  /**
     Prior for mixture model
  */
  class InitialMixturePriorPoisson : public InitialMixturePrior
  {

  protected:

    //! The Poisson frequency
    double Lambda;

    //! Vector containing the prefactors for each value of \f$k\in[0, nmix]\f$
    vector<double> Poisson;

    //! Uniform random variates
    Uniform * unit;

    //! Generate a variate by table look up
    int GenPoisson();

  public:
    
    /**@name Constructors */
    //@{

    //+ CLICONSTR
    //! For cloning only: makes an uninitialized prior
    InitialMixturePriorPoisson();

    //+ CLICONSTR StateInfo* double double clivectordist*
    /** The useful constructor
	@param si    is the StateInfo definition
	@param lambda is the Poisson mean
	@param alpha is the Dirichlet parameter
	@param vdist is the vector of distributions, one for each each dimension of the parameter vector
    */
    InitialMixturePriorPoisson(StateInfo *si, double lambda, 
			       double alpha, clivectordist *vdist);

    //! Destructor
    ~InitialMixturePriorPoisson();


    /**@name Distribution members
       The packing order in vectors is:
       Nmix weights, followed by Nmix vectors of Ndim elements
     */
    //@{
    //! Object factor (clone)
    InitialMixturePriorPoisson* New();
    //@}
    
    //! Sample: return weights and components separately
    void SamplePrior(unsigned& M, vector<double>& wght, 
		     vector< vector<double> >& phi);

    //! Sample: return a single parameter vector for subspace M
    void SamplePrior(unsigned M, vector<double>& V);

    #ifndef SWIG
    // AUTO GENERATED BY ../persistence/autopersist.py
    private:
    friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive & ar, const unsigned int version) {
        this->pre_serialize(ar, version);
         try {                                                         
          ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(InitialMixturePrior);            
          BIE_CATCH_BOOST_SERIALIZATION_EXCEPTION;                     
         }                                                             
         try {                                                         
          ar & BOOST_SERIALIZATION_NVP(Lambda);                        
          BIE_CATCH_BOOST_SERIALIZATION_EXCEPTION;                     
         }                                                             
         try {                                                         
          ar & BOOST_SERIALIZATION_NVP(Poisson);                        
          BIE_CATCH_BOOST_SERIALIZATION_EXCEPTION;                     
         }                                                             
         try {                                                         
          ar & BOOST_SERIALIZATION_NVP(unit);                        
          BIE_CATCH_BOOST_SERIALIZATION_EXCEPTION;                     
         }                                                             
        this->post_serialize(ar, version);
    }
    #endif

  };

}

#ifndef SWIG
BIE_CLASS_TYPE_INFO(BIE::InitialMixturePriorPoisson)
BIE_CLASS_EXPORT_KEY(BIE::InitialMixturePriorPoisson)
#endif
#endif
