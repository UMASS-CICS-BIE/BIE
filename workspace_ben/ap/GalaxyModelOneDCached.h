// This is really -*- C++ -*-

#ifndef GalaxyModelOneDCached_h
#define GalaxyModelOneDCached_h

#include <gaussQ.h>
#include <Model.h>
#include <Histogram1D.h>

#include <GalaxyModelOneD.h>

#include "Serializable.h"


namespace BIE {
  
  //+ CLICLASS GalaxyModelOneDCached SUPER GalaxyModelOneD
  //! Simple test galaxy model with one flux and caching.
  //! 
  //! This is model differs from GalaxyModelOneD by precomputing
  //! the line-of-sight quanties over a grid of parameter space
  //! a priori.  This may use too much memory to be practical, so
  //! be careful.
  class GalaxyModelOneDCached : public GalaxyModelOneD
  {
  public:

    //+ CLICONSTR int int SampleDistribution*
    //! Constructor 
    GalaxyModelOneDCached(int ndim, int mdim, SampleDistribution *_dist)
      : GalaxyModelOneD(ndim, mdim, _dist), current2(0) { }
      

    //! Destructor 
    ~GalaxyModelOneDCached();

    //! Compute normalization of tiles
    double NormEval(double x, double y, SampleDistribution *d);

    //! Main method returning source density
    vector<double> Evaluate(double x, double y, SampleDistribution *d);

    /** @name Global parameters */
    //@{

    //! Size of A grid
    static int numA;
    
    //! Size of H grid
    static int numH;
    
    //@}

  private:
    double dA, dH;

    mmapGalCMG cache2;
    mmapGalCMG::iterator mit2;
    CacheGalaxyModelGrid *current2;

    void generate(double L, double B, SampleDistribution *sd=NULL);
    void compute_bins();
    void manageCache(coordPair&);

    #ifndef SWIG
    // AUTO GENERATED BY ../persistence/autopersist.py
    protected:
    GalaxyModelOneDCached() {}
    private:
    friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive & ar, const unsigned int version) {
        this->pre_serialize(ar, version);
         try {                                                         
          ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(GalaxyModelOneD);            
          BIE_CATCH_BOOST_SERIALIZATION_EXCEPTION;                     
         }                                                             
         try {                                                         
          ar & BOOST_SERIALIZATION_NVP(numA);                        
          BIE_CATCH_BOOST_SERIALIZATION_EXCEPTION;                     
         }                                                             
         try {                                                         
          ar & BOOST_SERIALIZATION_NVP(numH);                        
          BIE_CATCH_BOOST_SERIALIZATION_EXCEPTION;                     
         }                                                             
         try {                                                         
          ar & BOOST_SERIALIZATION_NVP(dA);                        
          BIE_CATCH_BOOST_SERIALIZATION_EXCEPTION;                     
         }                                                             
         try {                                                         
          ar & BOOST_SERIALIZATION_NVP(dH);                        
          BIE_CATCH_BOOST_SERIALIZATION_EXCEPTION;                     
         }                                                             
         try {                                                         
          ar & BOOST_SERIALIZATION_NVP(cache2);                        
          BIE_CATCH_BOOST_SERIALIZATION_EXCEPTION;                     
         }                                                             
         try {                                                         
          ar & BOOST_SERIALIZATION_NVP(current2);                        
          BIE_CATCH_BOOST_SERIALIZATION_EXCEPTION;                     
         }                                                             
        this->post_serialize(ar, version);
    }
    #endif


  };

}

#ifndef SWIG
BIE_CLASS_TYPE_INFO(BIE::GalaxyModelOneDCached)
BIE_CLASS_EXPORT_KEY(BIE::GalaxyModelOneDCached)
#endif
#endif
