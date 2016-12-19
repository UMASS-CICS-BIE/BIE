// This is really -*- C++ -*-

#ifndef GalaxyModelOneD_h
#define GalaxyModelOneD_h

#include <deque>

#include <gaussQ.h>
#include <Model.h>
#include <Histogram1D.h>

#include <CacheGalaxyModel.h>

#include "Serializable.h"


namespace BIE {
  
  //+ CLICLASS GalaxyModelOneD SUPER Model
  //! Simple test galaxy model with one flux and line-of-sight caching 
  //! 
  //! The line-of-sight quantities that are parameter independent
  //! are precomputed and cached for each line of sight.
  //! 
  //! This is implemented with a hash map and a FIFO queue that is
  //! set by default to an infinite number of elements.  This may be
  //! changed using the CacheLimit() member function.
  class GalaxyModelOneD : public Model
  {
  public:

    //+ CLICONSTR int int SampleDistribution*
    //! Constructor 
    GalaxyModelOneD(int ndim, int mdim, SampleDistribution *_dist);

    //! Null constructor 
    GalaxyModelOneD();

    //! Destructor 
    ~GalaxyModelOneD();

    //+ CLIMETHOD void SetKnots int
    //! Reset number of integration points
    void SetKnots(int num) 
    { 
      NUM = num;
      delete intgr;
      intgr = new JacoQuad(NUM, ALPHA, BETA);
    }

    //+ CLIMETHOD void LogLineOfSight double
    //! Logarithmic spacing of integration knots along the line of sight
    void LogLineOfSight(double Smin)
    { 
      logs = true;
      smin = Smin;
    }

    //+ CLIMETHOD void LinearLineOfSight
    //! Linear spacing of integration knots along the line of sight
    void LinearLineOfSight()
    { 
      logs = false;
    }

    //+ CLIMETHOD void SetExtinction double double
    //! Radial and Vertical size of extinction slab
    void SetExtinction(double A, double Z) {A1 = A; Z1 = Z;}

    //+ CLIMETHOD void ExtinctionCoefficient double
    //! Set value of extinction coefficient (nominally k-band)
    void ExtinctionCoefficient(double ak) { AK = ak; }

    //+ CLIMETHOD void CacheLimit int
    //! Maximum number of cache elements (default=0, unlimited)
    void CacheLimit(int num) { cache_limit=num; }

    //+ CLIMETHOD void ModelLimits double double double double
    //! Set limiting values for model parameters
    void ModelLimits(double amin, double amax, double hmin, double hmax) {
      AMIN = amin;
      AMAX = amax;
      HMIN = hmin;
      HMAX = hmax;
    }

    //+ CLIMETHOD void MaximumDistance double
    //! Set line-of-sight extent
    void MaximumDistance(double rmax) { RMAX=rmax; }

    //! Initialize state dependent part of calculation
    //@{
    //! From State vector
    void Initialize(State&);
    //! From component weights and component parameter vectors
    void Initialize(vector<double>& w, vector< vector<double> >& p);
    //@}

    //! Reset model cache
    void ResetCache() {}

    //! Compute normalization of tiles
    virtual double NormEval(double x, double y, SampleDistribution *d=NULL);

    //! Compute normalization for point likelihood
    virtual double NormEval(double x, double y);

    //! Integration measure
    virtual double NormEvalMeasure(double x, double y) { return cos(y); }

    //! Main method returning source density
    vector<double> Evaluate(double x, double y, SampleDistribution *d=NULL);

    //! Identify labels for StateInfo
    vector<string> ParameterLabels();

    /** @name Global parameters */
    //@{
    
    //! Minimum magnitude (default 6);
    static double LMAG;

    //! Maximum magnitude (default 16);
    static double HMAG;

    //! Extinction scale length (default 20 [kpc])
    static double A1;

    //! Extinction scale height (default 100 [pc])
    static double Z1;

    //! Standard candle magnitude (default -4.0)
    static double K0;

    //! Std dev in standard candle magnitudes (default 0.25)
    static double SIGK;

    //! Number of integration knots for Jacobi quadrature (default 200)
    static int NUM;

    //! Jacobi quadrature parameters (default: 0, 0)
    static double ALPHA;
    static double BETA;

    //! Observers position (kpc) (default 8.0)
    static double R0;

    // Maximum radius (kpc) (default 20.0)
    static double RMAX;

    //! Extinction (mags) (default 0.1)
    static double AK;

    //! Limits for model parameters
    //@{
    //! Minimum scale length
    static double AMIN;		// Defaults: 0.2
    //! Maximum scale length
    static double AMAX;		//           8.0
    //! Minimum scale height
    static double HMIN;		//          50.0
    //! Maximum scale height
    static double HMAX;		//        1200.0
    //@}

    //! Log line of sight
    //@{
    //! Flag
    static bool logs;		// Default: false
    //! Minimum los
    static double smin;		// Default: 0.01
    //@}

    //@}

  protected:

    //! Constants
    //@{
    //! \f$\ln(10)\f$
    static double Log10;
    //! Number of radians per degree
    static double onedeg;
    //@}

    //! Cache limit exceeded
    bool cache_full;

    //! Maximum number of sight lines kept in cache
    unsigned int cache_limit;

    //! Maintain memory store for cache elements
    virtual void manageCache(coordPair&);

    //! Maximum number of components in mixture
    int M;

    //! Current number of components in mixture
    int Mcur;

    //! Number of model dimensions
    int Ndim;

    //! Constants related to record type.
    //@{
    //! Descriptor string for scale length
    static const char* LENGTH_FIELDNAME;
    //! Descriptor string for scale height
    static const char* HEIGHT_FIELDNAME;
    //! Descriptor strings for all other parameters
    static const char* PARAM_NAMES[];
    //@}

    //! Histogram components
    //@{
    //! Number of bins
    int nbins;
    //! Low value for bin edges
    vector<double> lowb;
    //! Low value for bin edges
    vector<double> highb;
    //@}

    //! Integrator
    JacoQuad *intgr;

    //! Component weights
    vector<double> wt;
    //! Parameter vectors for each component
    vector< vector<double> > pt;

    //! Compute values along the line-of-sight
    virtual void generate(double L, double B, SampleDistribution *sd=NULL);

    //! Compute the bin predictions for the given model paraemeters
    virtual void compute_bins();

    //! Flag is true if parameter values are in bounds
    bool good_bounds;

    //! Check for in-bounds parameters and set flag
    void check_bounds();

    //! The line-of-sight element cache
    mmapGalCM cache;
    
    //! Cache iterator
    mmapGalCM::iterator mit;

    //! List of cache coordinates (for management)
    deque<coordPair> cacheList;

    //! Current cache element
    CacheGalaxyModel *current;

    //! Evaluation type (e.g. binned or point)
    EvalType type;

    #ifndef SWIG
    // AUTO GENERATED BY ../persistence/autopersist.py
    private:
    friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive & ar, const unsigned int version) {
        this->pre_serialize(ar, version);
         try {                                                         
          ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(Model);            
          BIE_CATCH_BOOST_SERIALIZATION_EXCEPTION;                     
         }                                                             
         try {                                                         
          ar & BOOST_SERIALIZATION_NVP(LMAG);                        
          BIE_CATCH_BOOST_SERIALIZATION_EXCEPTION;                     
         }                                                             
         try {                                                         
          ar & BOOST_SERIALIZATION_NVP(HMAG);                        
          BIE_CATCH_BOOST_SERIALIZATION_EXCEPTION;                     
         }                                                             
         try {                                                         
          ar & BOOST_SERIALIZATION_NVP(A1);                        
          BIE_CATCH_BOOST_SERIALIZATION_EXCEPTION;                     
         }                                                             
         try {                                                         
          ar & BOOST_SERIALIZATION_NVP(Z1);                        
          BIE_CATCH_BOOST_SERIALIZATION_EXCEPTION;                     
         }                                                             
         try {                                                         
          ar & BOOST_SERIALIZATION_NVP(K0);                        
          BIE_CATCH_BOOST_SERIALIZATION_EXCEPTION;                     
         }                                                             
         try {                                                         
          ar & BOOST_SERIALIZATION_NVP(SIGK);                        
          BIE_CATCH_BOOST_SERIALIZATION_EXCEPTION;                     
         }                                                             
         try {                                                         
          ar & BOOST_SERIALIZATION_NVP(NUM);                        
          BIE_CATCH_BOOST_SERIALIZATION_EXCEPTION;                     
         }                                                             
         try {                                                         
          ar & BOOST_SERIALIZATION_NVP(ALPHA);                        
          BIE_CATCH_BOOST_SERIALIZATION_EXCEPTION;                     
         }                                                             
         try {                                                         
          ar & BOOST_SERIALIZATION_NVP(BETA);                        
          BIE_CATCH_BOOST_SERIALIZATION_EXCEPTION;                     
         }                                                             
         try {                                                         
          ar & BOOST_SERIALIZATION_NVP(R0);                        
          BIE_CATCH_BOOST_SERIALIZATION_EXCEPTION;                     
         }                                                             
         try {                                                         
          ar & BOOST_SERIALIZATION_NVP(RMAX);                        
          BIE_CATCH_BOOST_SERIALIZATION_EXCEPTION;                     
         }                                                             
         try {                                                         
          ar & BOOST_SERIALIZATION_NVP(AK);                        
          BIE_CATCH_BOOST_SERIALIZATION_EXCEPTION;                     
         }                                                             
         try {                                                         
          ar & BOOST_SERIALIZATION_NVP(AMIN);                        
          BIE_CATCH_BOOST_SERIALIZATION_EXCEPTION;                     
         }                                                             
         try {                                                         
          ar & BOOST_SERIALIZATION_NVP(AMAX);                        
          BIE_CATCH_BOOST_SERIALIZATION_EXCEPTION;                     
         }                                                             
         try {                                                         
          ar & BOOST_SERIALIZATION_NVP(HMIN);                        
          BIE_CATCH_BOOST_SERIALIZATION_EXCEPTION;                     
         }                                                             
         try {                                                         
          ar & BOOST_SERIALIZATION_NVP(HMAX);                        
          BIE_CATCH_BOOST_SERIALIZATION_EXCEPTION;                     
         }                                                             
         try {                                                         
          ar & BOOST_SERIALIZATION_NVP(logs);                        
          BIE_CATCH_BOOST_SERIALIZATION_EXCEPTION;                     
         }                                                             
         try {                                                         
          ar & BOOST_SERIALIZATION_NVP(smin);                        
          BIE_CATCH_BOOST_SERIALIZATION_EXCEPTION;                     
         }                                                             
         try {                                                         
          ar & BOOST_SERIALIZATION_NVP(Log10);                        
          BIE_CATCH_BOOST_SERIALIZATION_EXCEPTION;                     
         }                                                             
         try {                                                         
          ar & BOOST_SERIALIZATION_NVP(onedeg);                        
          BIE_CATCH_BOOST_SERIALIZATION_EXCEPTION;                     
         }                                                             
         try {                                                         
          ar & BOOST_SERIALIZATION_NVP(cache_limit);                        
          BIE_CATCH_BOOST_SERIALIZATION_EXCEPTION;                     
         }                                                             
         try {                                                         
          ar & BOOST_SERIALIZATION_NVP(M);                        
          BIE_CATCH_BOOST_SERIALIZATION_EXCEPTION;                     
         }                                                             
         try {                                                         
          ar & BOOST_SERIALIZATION_NVP(Mcur);                        
          BIE_CATCH_BOOST_SERIALIZATION_EXCEPTION;                     
         }                                                             
         try {                                                         
          ar & BOOST_SERIALIZATION_NVP(Ndim);                        
          BIE_CATCH_BOOST_SERIALIZATION_EXCEPTION;                     
         }                                                             
         try {                                                         
          ar & BOOST_SERIALIZATION_NVP(nbins);                        
          BIE_CATCH_BOOST_SERIALIZATION_EXCEPTION;                     
         }                                                             
         try {                                                         
          ar & BOOST_SERIALIZATION_NVP(lowb);                        
          BIE_CATCH_BOOST_SERIALIZATION_EXCEPTION;                     
         }                                                             
         try {                                                         
          ar & BOOST_SERIALIZATION_NVP(highb);                        
          BIE_CATCH_BOOST_SERIALIZATION_EXCEPTION;                     
         }                                                             
         try {                                                         
          ar & BOOST_SERIALIZATION_NVP(intgr);                        
          BIE_CATCH_BOOST_SERIALIZATION_EXCEPTION;                     
         }                                                             
         try {                                                         
          ar & BOOST_SERIALIZATION_NVP(wt);                        
          BIE_CATCH_BOOST_SERIALIZATION_EXCEPTION;                     
         }                                                             
         try {                                                         
          ar & BOOST_SERIALIZATION_NVP(pt);                        
          BIE_CATCH_BOOST_SERIALIZATION_EXCEPTION;                     
         }                                                             
         try {                                                         
          ar & BOOST_SERIALIZATION_NVP(good_bounds);                        
          BIE_CATCH_BOOST_SERIALIZATION_EXCEPTION;                     
         }                                                             
         try {                                                         
          ar & BOOST_SERIALIZATION_NVP(type);                        
          BIE_CATCH_BOOST_SERIALIZATION_EXCEPTION;                     
         }                                                             
        this->post_serialize(ar, version);
    }
    #endif

  };

}

#ifndef SWIG
BIE_CLASS_TYPE_INFO(BIE::GalaxyModelOneD)
BIE_CLASS_EXPORT_KEY(BIE::GalaxyModelOneD)
#endif
#endif
