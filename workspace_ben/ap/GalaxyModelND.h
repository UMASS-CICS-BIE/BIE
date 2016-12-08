// This is really -*- C++ -*-


#ifndef GalaxyModelND_h
#define GalaxyModelND_h

#include <gaussQ.h>
#include <Model.h>
#include <Histogram1D.h>
#include <BIEconfig.h>

#include <cmath>
#include <deque>
#include <utility>

#include <CacheGalaxyModel.h>

#include "Serializable.h"


namespace BIE {
  
  //+ CLICLASS GalaxyModelND SUPER Model
  //! Simple test galaxy model with one flux, one color and no caching
  class GalaxyModelND : public Model
  {
  public:

    //+ CLICONSTR int int SampleDistribution*
    //! Constructor 
    GalaxyModelND(int ndim, int mdim, SampleDistribution *histo);

    //! Destructor 
    ~GalaxyModelND();

    //+ CLIMETHOD void SetKnots int
    //! Reset number of integration points
    void SetKnots(int num) 
    { 
      NUM = num;
      delete intgr;
      intgr = new JacoQuad(NUM, ALPHA, BETA);
    }

    //+ CLIMETHOD void CacheLimit int
    //! Maximum number of cache elements (default=0, unlimited)
    void CacheLimit(int num) { cache_limit=num; }

    //+ CLIMETHOD void SetExtinction double double
    //! Radial and Vertical size if extinction slab
    void SetExtinction(double A, double Z) {A1 = A; Z1 = Z;}

    //! Initialize state dependent part of calculation
    //@{
    //@! From state vector
    void Initialize(State&);
    //@! From component weights and component parameter vectors
    void Initialize(vector<double>& w, vector< vector<double> >& p);
    //@}

    //! Reset model cache
    void ResetCache() {}

    //! Compute normalization of tiles (binned)
    virtual double NormEval(double x, double y, SampleDistribution *d);

    //! Line of sight eval for point normalization
    virtual double NormEval(double x, double y);

    //! Integration measure
    virtual double NormEvalMeasure(double x, double y) { return cos(y); }

    //! Main method returning source density
    vector<double> Evaluate(double x, double y, SampleDistribution *d);

    //! Identify labels for StateInfo
    vector<string> ParameterLabels();

    /** @name Global parameters */
    //@{

    //! Minimum apparent magnitude (default 6)
    static double LMAG;

    //! Minimum apparent magnitude (default 16)
    static double HMAG;

    //! Extinction scale length (default 20 [kpc])
    static double A1;

    //! Extinction scale height (default 100 [pc])
    static double Z1;

    //! Standard candle magnitude (default -4.0)
    static double K0;

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
    static double AMIN;		// Defaults: 0.2
    static double AMAX;		//           8.0
    static double HMIN;		//          50.0
    static double HMAX;		//        2000.0

    //! File defining stellar classes
    static string BASISDATA;	// Defaults: "basisdata.2d"

    //@}

  protected:

    //! Constants
    //@{
    //! \f$\ln(10)\f$
    static double Log10;
    //! Number of radian per degree
    static double onedeg;
    //@}

    //! Maximum number of components in mixture
    int M;

    //! Current number of components in mixture
    int Mcur;

    //! Number of model dimensions
    int Ndim;

    //! Maximum number of sight lines kept in cache
    unsigned int cache_limit;

    //! Maintain memory store for cache elements
    virtual void manageCache(coordPair&);

    //! Constants related to record type.
    //@{
    //! Length name
    static const char* LENGTH_FIELDNAME;
    //! Height name
    static const char* HEIGHT_FIELDNAME;
    //! Other paramter names
    static const char* PARAM_NAMES[];
    //@}

    //! Histogram components
    //@{
    //! Number of bins
    int nbins;
    //! Low values for bin boundaries
    vector<dvector> lowb;
    //! High values for bin boundaries
    vector<dvector> highb;
    //@}

    //! Point components
    vector<double> flux;

    //! Number of data dimensions
    int Nflux;

    //! Low flux limit
    vector<double> lolim;

    //! High flux limit
    vector<double> hilim;

    //! Stellar components
    //@{
    //! Number of components
    int nparam;
    //! Weights for each component
    vector<double> w;
    //! Centers for each component
    vector<dvector> pos;
    //! Widths for each component
    vector<dvector> sig;
    //@}

    //! Integrator
    JacoQuad *intgr;

    //! Component weights
    vector<double> wt;
    //! Parameter vectors for each comonent
    vector< vector<double> > pt;
    //! Work vector for bin computation
    vector<double> work;

    //! Compute line-of-site quantities
    virtual void generate(double L, double B, SampleDistribution *sd=NULL);
    //! Compute bins for current parameter vector
    virtual void compute_bins();

    //! Flag is true if parameter value are in bounds
    bool good_bounds;

    //! Check parameter bounds and set flag
    void check_bounds();

    //! The line-of-sight element cache
    mmapGalCM cache;
    //! Line-of-sight cache iterator
    mmapGalCM::iterator mit;
    //! List of cche keys (for cache maintenance)
    deque<coordPair> cacheList;
    //! Current line of sight
    CacheGalaxyModel *current;
    //! Number of cache misses
    int missed;

    //! Evaluation type (e.g. binned or point)
    EvalType type;

    #ifndef SWIG
    // AUTO GENERATED BY ../persistence/autopersist.py
    protected:
    GalaxyModelND() {}
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
          ar & BOOST_SERIALIZATION_NVP(BASISDATA);                        
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
          ar & BOOST_SERIALIZATION_NVP(cache_limit);                        
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
          ar & BOOST_SERIALIZATION_NVP(flux);                        
          BIE_CATCH_BOOST_SERIALIZATION_EXCEPTION;                     
         }                                                             
         try {                                                         
          ar & BOOST_SERIALIZATION_NVP(Nflux);                        
          BIE_CATCH_BOOST_SERIALIZATION_EXCEPTION;                     
         }                                                             
         try {                                                         
          ar & BOOST_SERIALIZATION_NVP(lolim);                        
          BIE_CATCH_BOOST_SERIALIZATION_EXCEPTION;                     
         }                                                             
         try {                                                         
          ar & BOOST_SERIALIZATION_NVP(hilim);                        
          BIE_CATCH_BOOST_SERIALIZATION_EXCEPTION;                     
         }                                                             
         try {                                                         
          ar & BOOST_SERIALIZATION_NVP(nparam);                        
          BIE_CATCH_BOOST_SERIALIZATION_EXCEPTION;                     
         }                                                             
         try {                                                         
          ar & BOOST_SERIALIZATION_NVP(w);                        
          BIE_CATCH_BOOST_SERIALIZATION_EXCEPTION;                     
         }                                                             
         try {                                                         
          ar & BOOST_SERIALIZATION_NVP(pos);                        
          BIE_CATCH_BOOST_SERIALIZATION_EXCEPTION;                     
         }                                                             
         try {                                                         
          ar & BOOST_SERIALIZATION_NVP(sig);                        
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
          ar & BOOST_SERIALIZATION_NVP(work);                        
          BIE_CATCH_BOOST_SERIALIZATION_EXCEPTION;                     
         }                                                             
         try {                                                         
          ar & BOOST_SERIALIZATION_NVP(good_bounds);                        
          BIE_CATCH_BOOST_SERIALIZATION_EXCEPTION;                     
         }                                                             
         try {                                                         
          ar & BOOST_SERIALIZATION_NVP(cache);                        
          BIE_CATCH_BOOST_SERIALIZATION_EXCEPTION;                     
         }                                                             
         try {                                                         
          ar & BOOST_SERIALIZATION_NVP(cacheList);                        
          BIE_CATCH_BOOST_SERIALIZATION_EXCEPTION;                     
         }                                                             
         try {                                                         
          ar & BOOST_SERIALIZATION_NVP(current);                        
          BIE_CATCH_BOOST_SERIALIZATION_EXCEPTION;                     
         }                                                             
         try {                                                         
          ar & BOOST_SERIALIZATION_NVP(missed);                        
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
BIE_CLASS_TYPE_INFO(BIE::GalaxyModelND)
BIE_CLASS_EXPORT_KEY(BIE::GalaxyModelND)
#endif
#endif
