// This is really -*- C++ -*-

#ifndef GalaxyModelOneDHashed_h
#define GalaxyModelOneDHashed_h

#include <gaussQ.h>
#include <Model.h>
#include <Histogram1D.h>
#include <CacheGalaxyModel.h>

#include "Serializable.h"


namespace BIE {
  
  //+ CLICLASS GalaxyModelOneDHashed SUPER Model
  //! Simple test galaxy model with one flux and hash map caching 
  class GalaxyModelOneDHashed : public Model
  {
  public:

    //+ CLICONSTR int int SampleDistribution*
    //! Constructor 
    GalaxyModelOneDHashed(int ndim, int mdim, SampleDistribution *_dist);

    //! Destructor 
    ~GalaxyModelOneDHashed();

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
    //! Radial and Vertical size if extinction slab
    void SetExtinction(double A, double Z) {A1 = A; Z1 = Z;}

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

    //! Main method returning source density
    vector<double> Evaluate(double x, double y, SampleDistribution *d=NULL);

    //! Identify labels for StateInfo
    vector<string> ParameterLabels();

    /** @name Global parameters */
    //@{

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

    /** Jacobi quadrature parameters (default: 0, 0)
	@see GalaxyModelND for more details
    */
    //@{
    //! Alpha exponent
    static double ALPHA;
    //! Beta exponent
    static double BETA;
    //@}

    //! Observers position (kpc) (default 8.0)
    static double R0;

    //! Maximum radius (kpc) (default 20.0)
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

    //! Maximum number of components in the mixture
    int M;

    //! Current number of components in the mixture
    int Mcur;

    //! Number of model dimensions
    int Ndim;

    //! Constants related to record type.
    //@{
    //! Description of scale-length field
    static const char* LENGTH_FIELDNAME;
    //! Description of scale-height field
    static const char* HEIGHT_FIELDNAME;
    //! Description of all other fields
    static const char* PARAM_NAMES[];
    //@}

    //! Histogram components
    //@{
    //! Number of bins
    int nbins;
    //! Low values of bin edges
    vector<double> lowb;
    //! High values of bin edges
    vector<double> highb;
    //@}

    //! Integrator
    JacoQuad *intgr;

    //! Component weights
    vector<double> wt;

    //! Parameter vectors for each component
    vector< vector<double> > pt;

    //! Compute line-of-sight values 
    virtual CacheGalaxyModel* 
      generate(const coordPair &, SampleDistribution *sd=NULL);

    //! Compute bin values for current parameter vector
    virtual CacheGalaxyModel* 
      compute_bins(const coordPair &, SampleDistribution *sd);

    //! True if paramter vector values are in bounds
    bool good_bounds;

    //! Check that parameter values are in bounds and set flag
    void check_bounds();

    //! Cache element hash map
    mmapGalCM cache;

    //! Last coordinate evaluated
    coordPair lastP;

    //! Tally hash misses
    int missed;

    //! Type of evaluation (e.g. point or bined)
    EvalType type;

    #ifndef SWIG
    // AUTO GENERATED BY ../persistence/autopersist.py
    protected:
    GalaxyModelOneDHashed() {}
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
          ar & BOOST_SERIALIZATION_NVP(cache);                        
          BIE_CATCH_BOOST_SERIALIZATION_EXCEPTION;                     
         }                                                             
         try {                                                         
          ar & BOOST_SERIALIZATION_NVP(lastP);                        
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
BIE_CLASS_TYPE_INFO(BIE::GalaxyModelOneDHashed)
BIE_CLASS_EXPORT_KEY(BIE::GalaxyModelOneDHashed)
#endif
#endif
