// This is really -*- C++ -*-



#ifndef GalaxyModelNDCached_h
#define GalaxyModelNDCached_h

#include <gaussQ.h>
#include <Model.h>
#include <Histogram1D.h>

#include <GalaxyModelND.h>

#include "Serializable.h"


namespace BIE {
  
  //+ CLICLASS GalaxyModelNDCached SUPER GalaxyModelND
  /**
     Simple test galaxy model with one flux and no caching
  */
  class GalaxyModelNDCached : public GalaxyModelND
  {
  public:

    //+ CLICONSTR int int SampleDistribution*
    //! Constructor 
    GalaxyModelNDCached(int ndim, int mdim, SampleDistribution *_dist);

    //! Destructor 
    ~GalaxyModelNDCached();

    //! Compute normalization of tiles
    double NormEval(double x, double y, SampleDistribution *d);

    //! Main method returning source density
    vector<double> EvaluateBinned(double x, double y, SampleDistribution *d);

    /** @name Global parameters */
    //@{

    //! Size of A grid
    static int numA;
    
    //! Size of H grid
    static int numH;
    
    //@}

  private:
    double dA, dH;
    vector<CacheGalaxyModelGrid*> cache;

    void generate(double L, double B, SampleDistribution *d);
    void compute_bins();

    mmapGalCMG cache2;
    mmapGalCMG::iterator mit2;
    CacheGalaxyModelGrid *current2;

    #ifndef SWIG
    // AUTO GENERATED BY ../persistence/autopersist.py
    protected:
    GalaxyModelNDCached() {}
    private:
    friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive & ar, const unsigned int version) {
        this->pre_serialize(ar, version);
         try {                                                         
          ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(GalaxyModelND);            
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
          ar & BOOST_SERIALIZATION_NVP(cache);                        
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
BIE_CLASS_TYPE_INFO(BIE::GalaxyModelNDCached)
BIE_CLASS_EXPORT_KEY(BIE::GalaxyModelNDCached)
#endif
#endif
