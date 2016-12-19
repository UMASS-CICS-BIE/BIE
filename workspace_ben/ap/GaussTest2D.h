// This is -*- C++ -*-

#ifndef GaussTest2D_h
#define GaussTest2D_h

#include "Serializable.h"


#include <LikelihoodFunction.h>

namespace BIE {

  //+ CLICLASS GaussTest2D SUPER LikelihoodFunction
  /** A "user-defined" likelihood function for testing a image object
      similar to a galaxy image but without any intrisic magnitude.
      The singal-to-noise ratio is determined by the number of points
      used to define the image.  For example, if we define the signal
      to noise ratio as the counts over the root variance of the counts
      for the inner half-light, then a signal to noise ratio of 100
      would obtain for
      \f[
      100 = \frac{N/2}{\sqrt{N/2}} = \sqrt(N/2)
      \f]
      which yields \f$N=200000\f$.
   
      By default, the image "data" is the combination of a single
      2-dimensional Gaussian distribution with center at 0.5 and
      variance of 0.03 in each dimension.  The default model may be
      changed at construction.
   
      The variance may be modeled (dim=2) or fixed (dim=1) and one may
      include an arbitrary position angle (dim=3) using the
      SetDim(dim) member.  dim=3 implies variable center, variance and
      position angle.

      Finally, two components packed into a single component, in order
      to keep the centers coincident and fix the shape of the second
      component.  The first component is the relative weight of the
      two components followed by \f$x_c\f$ and \f$y_c\f$ if dim=1,
      \f$x_c, \sigma_x^2, y_c, \sigma_y^2\f$ if dim=2 and \f$x_c,
      \sigma_x^2, y_c, \sigma_y^2, \phi_1, \phi_2\f$ if dim=3.

      @ingroup likefunc
  */
  class GaussTest2D : public LikelihoodFunction 
  {

  private:

    vector < vector<unsigned> > fdata;
    vector < vector<unsigned> > ldata;

    unsigned ncomp, dim, mdim, N;
    unsigned nbins, nbins1;
    unsigned levels, total;
    int llevel;

    vector<double> weights, var2, angles;
    vector< vector<double> > centers, variance;
    bool fixed2;

    unsigned lastl;

    void compute_level();

 public:


    //+ CLICONSTR
    //! Null constructor
    GaussTest2D();

    //+ CLICONSTR int int int
    //+ CLICONSTR int int int clivectord* clivectord*
    //+ CLICONSTR int int int clivectord* clivectord* clivectord*
    //+ CLICONSTR int int int clivectord* clivectord* clivectord* clivectord*

    /** Constructor with
	@param N0 number of points
	@param Nbins number of bins 
	@param Levels number of partitions

	If Levels is equal to smaller than 1, a single partition is used.
    */
    GaussTest2D(int N0, int Nbins, int Levels);


    //@{
    /** Constructor with
	@param N0 number of points
	@param Nbins number of bins 
	@param Levels number of partitions
	@param cen0 are the centers of the components
	@param var0 are the variance values of the components

	If Levels is equal to smaller than 1, a single partition is
	used.  The dimension of the vectors (which must agree)
	implicitly specifies the dimension.
    */
    GaussTest2D(int N0, int Nbins, int Levels, 
		clivectord* cen0, clivectord* var0);

    GaussTest2D(int N0, int Nbins, int Levels, 
		vector<double>* cen0, vector<double>* var0);
    //@}

    //@{
    /** Constructor with
	@param N0 number of points
	@param Nbins number of bins 
	@param Levels number of partitions
	@param wght is a multicomponent weight vector
	@param cen0 are the centers of the components
	@param var0 are the variance values of the components

	If Levels is equal to smaller than 1, a single partition is
	used.  The dimension of the vectors (which must agree)
	implicitly specifies the dimension.  If dim(wght) > 1, the
	vectors are assumed to be appended serially.
    */
    GaussTest2D(int N0, int Nbins, int Levels, clivectord* wght,
		    clivectord* cen0, clivectord* var0);

    GaussTest2D(int N0, int Nbins, int Levels, vector<double>* wght,
		    vector<double>* cen0, vector<double>* var0);
    //@}

    //@{
    /** Constructor with
	@param N0 number of points
	@param Nbins number of bins 
	@param Levels number of partitions
	@param wght is a multicomponent weight vector
	@param cen0 are the centers of the components
	@param var0 are the variance values of the components
	@param phi0 are the inclination angles of the x-axis

	If Levels is equal to smaller than 1, a single partition is
	used.  The dimension of the vectors (which must agree)
	implicitly specifies the dimension.  If dim(wght) > 1, the
	vectors are assumed to be appended serially.
    */
    GaussTest2D(int N0, int Nbins, int Levels, clivectord* wght,
		clivectord* cen0, clivectord* var0, 
		clivectord* phi0);

    GaussTest2D(int N0, int Nbins, int Levels, vector<double>* wght,
		vector<double>* cen0, vector<double>* var0, 
		vector<double>* phi0);
    //@}
    
    //+ CLIMETHOD void SetDim int
    //! Set the model dimension (currently either 1, 2 or 3)
    void SetDim(int n);
    
    //+ CLIMETHOD void SetFix2
    /** Set the centers of all components to be coincident.  The
	center will be carried by the extended component in the
	mixture structure.
    */
    void SetFix2() { fixed2 = true; }

    //+ CLIMETHOD void PrintData
    //! Print data to a file for the current level
    void PrintData();

    //! This is likelihood function
    double LikeProb(std::vector<double> &z, SampleDistribution* sd, 
		    double norm, Tile *t, State *s, int indx)
    { return LocalLikelihood(s); }

    //! Parameter labels
    vector<string> ParameterLabels();

    //! Position labels if fixed center model
    vector<string> ExtendedLabels();

 protected:

    //! Make the synthetic data
    void makeSyntheticData();

    //! This is likelihood function
    double LocalLikelihood(State* s);

    #ifndef SWIG
    // AUTO GENERATED BY ../persistence/autopersist.py
    private:
    friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive & ar, const unsigned int version) {
        this->pre_serialize(ar, version);
         try {                                                         
          ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(LikelihoodFunction);            
          BIE_CATCH_BOOST_SERIALIZATION_EXCEPTION;                     
         }                                                             
         try {                                                         
          ar & BOOST_SERIALIZATION_NVP(fdata);                        
          BIE_CATCH_BOOST_SERIALIZATION_EXCEPTION;                     
         }                                                             
         try {                                                         
          ar & BOOST_SERIALIZATION_NVP(ldata);                        
          BIE_CATCH_BOOST_SERIALIZATION_EXCEPTION;                     
         }                                                             
         try {                                                         
          ar & BOOST_SERIALIZATION_NVP(ncomp);                        
          BIE_CATCH_BOOST_SERIALIZATION_EXCEPTION;                     
         }                                                             
         try {                                                         
          ar & BOOST_SERIALIZATION_NVP(dim);                        
          BIE_CATCH_BOOST_SERIALIZATION_EXCEPTION;                     
         }                                                             
         try {                                                         
          ar & BOOST_SERIALIZATION_NVP(mdim);                        
          BIE_CATCH_BOOST_SERIALIZATION_EXCEPTION;                     
         }                                                             
         try {                                                         
          ar & BOOST_SERIALIZATION_NVP(N);                        
          BIE_CATCH_BOOST_SERIALIZATION_EXCEPTION;                     
         }                                                             
         try {                                                         
          ar & BOOST_SERIALIZATION_NVP(nbins);                        
          BIE_CATCH_BOOST_SERIALIZATION_EXCEPTION;                     
         }                                                             
         try {                                                         
          ar & BOOST_SERIALIZATION_NVP(nbins1);                        
          BIE_CATCH_BOOST_SERIALIZATION_EXCEPTION;                     
         }                                                             
         try {                                                         
          ar & BOOST_SERIALIZATION_NVP(levels);                        
          BIE_CATCH_BOOST_SERIALIZATION_EXCEPTION;                     
         }                                                             
         try {                                                         
          ar & BOOST_SERIALIZATION_NVP(total);                        
          BIE_CATCH_BOOST_SERIALIZATION_EXCEPTION;                     
         }                                                             
         try {                                                         
          ar & BOOST_SERIALIZATION_NVP(llevel);                        
          BIE_CATCH_BOOST_SERIALIZATION_EXCEPTION;                     
         }                                                             
         try {                                                         
          ar & BOOST_SERIALIZATION_NVP(weights);                        
          BIE_CATCH_BOOST_SERIALIZATION_EXCEPTION;                     
         }                                                             
         try {                                                         
          ar & BOOST_SERIALIZATION_NVP(var2);                        
          BIE_CATCH_BOOST_SERIALIZATION_EXCEPTION;                     
         }                                                             
         try {                                                         
          ar & BOOST_SERIALIZATION_NVP(angles);                        
          BIE_CATCH_BOOST_SERIALIZATION_EXCEPTION;                     
         }                                                             
         try {                                                         
          ar & BOOST_SERIALIZATION_NVP(centers);                        
          BIE_CATCH_BOOST_SERIALIZATION_EXCEPTION;                     
         }                                                             
         try {                                                         
          ar & BOOST_SERIALIZATION_NVP(variance);                        
          BIE_CATCH_BOOST_SERIALIZATION_EXCEPTION;                     
         }                                                             
         try {                                                         
          ar & BOOST_SERIALIZATION_NVP(fixed2);                        
          BIE_CATCH_BOOST_SERIALIZATION_EXCEPTION;                     
         }                                                             
         try {                                                         
          ar & BOOST_SERIALIZATION_NVP(lastl);                        
          BIE_CATCH_BOOST_SERIALIZATION_EXCEPTION;                     
         }                                                             
        this->post_serialize(ar, version);
    }
    #endif

  };

} // namespace BIE

#ifndef SWIG
BIE_CLASS_TYPE_INFO(BIE::GaussTest2D)
BIE_CLASS_EXPORT_KEY(BIE::GaussTest2D)
#endif
#endif
