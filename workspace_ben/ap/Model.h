// This is really -*- C++ -*-

#ifndef Model_h
#define Model_h

#include <new>
#include <typeinfo>

#include <gaussQ.h>
#include <BIEconfig.h>
#include <RecordType.h>
#include <Distribution.h>

#include "Serializable.h"


namespace BIE {
  
  //+ CLICLASS Model
  //! Abstract: determines partial distribution of source density for
  //! given cell.  The ListIterator interface allows the model to
  //! cache configuration information.
  //! 
  //! For each likelihood evaluation, the model is called first to get
  //! compute the overall normalizaton and then to compute the
  //! likelihood for each datum.
  //! 
  //! For tiled distributions, the measure is set by the tile and the
  //! tessellation itself is used to perform the normalization integral.
  //! The intermediate values of the model prediction for each (x, y)
  //! may be cached during the normalization computation and reused by
  //! the likelihood computation for efficiency.
  //! 
  //! For point distributions, there is no natural measure and the
  //! normalization must be performed with an independent quadrature.
  //! Since this is done multiple times, the implementor may consider
  //! creating two caches, one for the normalization and one for
  //! likelihood evaluations.  Two quadrature schemes are currently
  //! available: two-dimensional Gauss-Legendre integration and a
  //! recursive, adaptive cubature scheme with error control.
  //! 
  //! @see QuadTreeIntegrator for more details on the cubature scheme
  class Model: public Serializable
  {
  public:

    //! Use quadtree integrator for point norm in serial mode (default: true)
    static bool quadtree;

    //! Maximum number of recursion levels for quadtree (default: 16)
    static int maxlevels;

    //! Tolerance for quadtree (default: 0.0001)
    static double qeps;

    //! Maximum x grid size for quadtree (default: 0.5)
    static double dX0;

    //! Maximum y grid size for quadtree (default: 0.1)
    static double dY0;

    //! Minimum x grid size for quadtree (default: 0.005)
    static double dX1;

    //! Minimum y grid size for quadtree (default: 0.001)
    static double dY1;

    //! Number of x-direction integration points for point norm (default: 20)
    static int numX;

    //! Number of y-direction integration points for point norm (default: 20)
    static int numY;

    virtual ~Model() 
    { 
      delete iX; delete iY;
      if (parametertype) {delete parametertype;} 
    }
    
    //+ CLIMETHOD RecordType* getParameterType
    //! Returns a record type describing the parameters of the model.
    RecordType * getParameterType() { return parametertype; }  
    
    //+ CLIMETHOD void setNormKnots double double
    /** Change the default Gauss-Legendre knot numbers for the point-type
	normalization integral */
    void setNormKnots(int nx, int ny) { numX = nx; numY = ny; }
    
    //+ CLIMETHOD void setQuadTreeParams double double double double double int
    /** Change the default QuadTree parameters
	@param dx0 is the maximum allowed grid space in the x direction
	@param dy0 is the maximum allowed grid space in the x direction
	@param dx1 is the minimum allowed grid space in the x direction
	@param dy1 is the minimum allowed grid space in the x direction
	@param eps is the desired tolerance
	@param mlev is the limit on recursive bisection of each dimension
    */
    void setQuadTreeParams(double dx0, double dy0, double dx1, double dy1,
			   double eps, int mlev);
    
    //! Initialize for particular state
    virtual void Initialize(State&) = 0;

    //! Hook for implementation dependent caching
    virtual void ResetCache() {};

    //! Integrated norm (e.g. for point likelihood)
    virtual double NormEval(double xmin, double xmax, 
			    double ymin, double ymax, bool serial);

    //! Contribution to norm at this line of sight
    virtual double NormEval(double x, double y, SampleDistribution *d) = 0;

    /** Integrated norm (e.g. for point likelihood). This should be
	overridden by the inheriting class or never called (e.g. if
	the normalization integral is overridden above). Send a
	message to the programmer this is overlooked.  This could be
	made pure virtual to enforce this behavior at compile time,
	but this approach saves the inheriting class from defining a
	dummy implementation that will never be used.
    */
    virtual double NormEval(double X, double Y) {
      string msg  = 
	"Model classes must implement a NormEval(double, double) method for\n"\
	"handling point type data or override the normalization calculation\n"\
	"NormEval(double, double, double, double)";
      throw InternalError(msg, __FILE__, __LINE__);
      return 0.0; 
    }

    //! Integration measure
    virtual double NormEvalMeasure(double X, double Y) { return 1.0; }

    //! Main method for returning source density
    virtual vector<double> Evaluate(double x, double y, SampleDistribution* d);

    /** Specific evaulate functions to be supplied by derived classes
	Either Evaluate or should be overloaded by the derived class or
	the following two classes need to be overloaded.
    */
    //@{

    //! Contribution to each bin from this line of sight
    virtual vector<double> EvaluateBinned(double x, double y, 
					  BinnedDistribution* d)
    { return vector<double>(1, 0.0); }

    //! Contribution to each point from this line of sight
    virtual vector<double> EvaluatePoint(double x, double y, 
					 PointDistribution* d)
    { return vector<double>(1, 0.0); }

    //@}

    //! Identify labels for StateInfo (mixture or non-mixture component)
    virtual vector<string> ParameterLabels() { return vector<string>(0); }

    //! Identify labels for StateInfo (extended component for mixture model)
    virtual vector<string> ExtendedLabels() { return vector<string>(0); }

    //! Enum for distribution types
    enum EvalType {binned, point};

  protected:

    //! Quadrature for point tessellation norm
    //@{
    LegeQuad *iX, *iY;
    //@}

    /** Constructor called implicitly by subclasses - makes sure the
	record type is 0, so that destructor knows not to delete the type
	if construction of the subclass fails before typoe creation .
    */
    Model() { iX = 0; iY = 0; parametertype = 0; }

    //! Record type of a model's parameters.
    RecordType * parametertype;

    //@{
    //! Create a record type for a mixture model from a list of parameters.
    RecordType* createMixParameterType
      (const char ** paramnames, int numparams, int mixturedim);
    RecordType* createMixParameterType
      (const vector<string>& paramnames, int numparams, int mixturedim);
    //@}
    
    #ifndef SWIG
    // AUTO GENERATED BY ../persistence/autopersist.py
    private:
    friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive & ar, const unsigned int version) {
        this->pre_serialize(ar, version);
         try {                                                         
          ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(Serializable);            
          BIE_CATCH_BOOST_SERIALIZATION_EXCEPTION;                     
         }                                                             
         try {                                                         
          ar & BOOST_SERIALIZATION_NVP(iX);                        
          BIE_CATCH_BOOST_SERIALIZATION_EXCEPTION;                     
         }                                                             
         try {                                                         
          ar & BOOST_SERIALIZATION_NVP(iY);                        
          BIE_CATCH_BOOST_SERIALIZATION_EXCEPTION;                     
         }                                                             
         try {                                                         
          ar & BOOST_SERIALIZATION_NVP(parametertype);                        
          BIE_CATCH_BOOST_SERIALIZATION_EXCEPTION;                     
         }                                                             
        this->post_serialize(ar, version);
    }
    #endif

  };
} // namespace BIE

#ifndef SWIG
BIE_CLASS_ABSTRACT(BIE::Model)

BIE_CLASS_EXPORT_KEY(BIE::Model)
#endif
#endif
