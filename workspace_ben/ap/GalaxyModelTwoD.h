// This is really -*- C++ -*-



#ifndef GalaxyModelTwoD_h
#define GalaxyModelTwoD_h

#include <gaussQ.h>
#include <Model.h>
#include <Histogram1D.h>
#include <CacheGalaxyModel.h>

#include "Serializable.h"


namespace BIE {

  //+ CLICLASS GalaxyModelTwoD SUPER Model
  //! Simple test galaxy model with one flux, one color and no caching
  class GalaxyModelTwoD : public Model
  {
  public:

    //+ CLICONSTR int int SampleDistribution*
    //! Constructor 
    GalaxyModelTwoD(int ndim, int mdim, SampleDistribution *_dist);

    //! Destructor 
    ~GalaxyModelTwoD();

    //+ CLIMETHOD void SetKnots int
    //! Reset number of integration points
    void SetKnots(int num) 
    { 
      NUM = num;
      delete intgr;
      intgr = new JacoQuad(NUM, ALPHA, BETA);
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

    //! Compute normalization for point likelihood
    virtual double NormEval(double x, double y);

    //! Compute normalization of tiles
    virtual double NormEval(double x, double y, SampleDistribution *d);

    //! Integration measure
    virtual double NormEvalMeasure(double x, double y) { return cos(y); }

    //! Main method returning source density
    vector<double> Evaluate(double x, double y, SampleDistribution *d);

    /// Identify labels for StateInfo
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

    /** Jacobi quadrature parameters (default: 0, 0)
	@see GalaxyModelND for more details
    */
    //@{
    //! Alpha exponent value
    static double ALPHA;
    //! Beta exponent value
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
    //! Maxiimum scale height
    static double HMAX;		//        1200.0
    //@}

    //! File defining stellar classes
    static string BASISDATA;	// Defaults: "basisdata.2d"

    //@}

  protected:

    //! Constants
    //@{
    //! \f$\ln(10)\f$
    static double Log10;
    //! Radians per degree
    static double onedeg;
    //@}
    
    //! Total number of components in mixture
    int M;

    //! Current number of components in mixture
    int Mcur;

    //! Number of model dimensions
    int Ndim;
     
    //! Constants related to record type.
    //@{
    //! String descriptor for scale length
    static const char* LENGTH_FIELDNAME;
    //! String descriptor for scale height
    static const char* HEIGHT_FIELDNAME;
    //! String descriptors for remaining parameters
    static const char* PARAM_NAMES[];
    //@}

    //! Histogram components
    //@{
    //! Number of bins
    int nbins;
    //! Low values for first flux dimension
    vector<double> lowb1;
    //! High values for first flux dimension
    vector<double> highb1;
    //! Low values for second flux dimension
    vector<double> lowb2;
    //! High values for second flux dimension
    vector<double> highb2;
    //@}

    //! Point components
    vector<double> flux;

    //! Stellar components
    //@{
    //! Number of parameters for the model
    int nparam;
    //! Centers for first dimension
    vector<double> x;
    //! Centers for second dimension
    vector<double> y;
    //! Widths for first dimension
    vector<double> sx;
    //! Widths for second dimension
    vector<double> sy;
    //! Weights for each component
    vector<double> w;
    //@}

    //! Integrator
    JacoQuad *intgr;

    //! Weights for each component in the mixture
    vector<double> wt;

    //! Parameter vector for each component in the mixture
    vector< vector<double> >pt;

    //! Compute line-of-sight quantities for these coordinates
    virtual void generate(double L, double B, SampleDistribution *d);

    //! Compute predicted bin values for current parameter vector
    virtual void compute_bins();

    //! Flag is true if parameter values are in bounds
    bool good_bounds;

    //! Check parameter bounds and set flag
    void check_bounds();

    //! Line-of-sight element cache
    mmapGalCM cache;

    //! Cache iterator
    mmapGalCM::iterator mit;

    //! Current cache element
    CacheGalaxyModel *current;

    //! Type of evaluation (e.g. point or binned)
    EvalType type;

    #ifndef SWIG
    // AUTO GENERATED BY ../persistence/autopersist.py
    protected:
    GalaxyModelTwoD() {}
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
          ar & BOOST_SERIALIZATION_NVP(nbins);                        
          BIE_CATCH_BOOST_SERIALIZATION_EXCEPTION;                     
         }                                                             
         try {                                                         
          ar & BOOST_SERIALIZATION_NVP(lowb1);                        
          BIE_CATCH_BOOST_SERIALIZATION_EXCEPTION;                     
         }                                                             
         try {                                                         
          ar & BOOST_SERIALIZATION_NVP(highb1);                        
          BIE_CATCH_BOOST_SERIALIZATION_EXCEPTION;                     
         }                                                             
         try {                                                         
          ar & BOOST_SERIALIZATION_NVP(lowb2);                        
          BIE_CATCH_BOOST_SERIALIZATION_EXCEPTION;                     
         }                                                             
         try {                                                         
          ar & BOOST_SERIALIZATION_NVP(highb2);                        
          BIE_CATCH_BOOST_SERIALIZATION_EXCEPTION;                     
         }                                                             
         try {                                                         
          ar & BOOST_SERIALIZATION_NVP(flux);                        
          BIE_CATCH_BOOST_SERIALIZATION_EXCEPTION;                     
         }                                                             
         try {                                                         
          ar & BOOST_SERIALIZATION_NVP(x);                        
          BIE_CATCH_BOOST_SERIALIZATION_EXCEPTION;                     
         }                                                             
         try {                                                         
          ar & BOOST_SERIALIZATION_NVP(y);                        
          BIE_CATCH_BOOST_SERIALIZATION_EXCEPTION;                     
         }                                                             
         try {                                                         
          ar & BOOST_SERIALIZATION_NVP(sx);                        
          BIE_CATCH_BOOST_SERIALIZATION_EXCEPTION;                     
         }                                                             
         try {                                                         
          ar & BOOST_SERIALIZATION_NVP(sy);                        
          BIE_CATCH_BOOST_SERIALIZATION_EXCEPTION;                     
         }                                                             
         try {                                                         
          ar & BOOST_SERIALIZATION_NVP(w);                        
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
          ar & BOOST_SERIALIZATION_NVP(current);                        
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
BIE_CLASS_TYPE_INFO(BIE::GalaxyModelTwoD)
BIE_CLASS_EXPORT_KEY(BIE::GalaxyModelTwoD)
#endif
#endif
