// This is -*- C++ -*-

#ifndef GaussTestMultiD_h
#define GaussTestMultiD_h

#include "Serializable.h"


#include <Uniform.h>
#include <Normal.h>
#include <Gamma.h>
#include <BetaRDist.h>
#include <LikelihoodFunction.h>
#include <CDFGenerator.h>

namespace BIE {

  //+ CLICLASS GaussTestMultiD SUPER LikelihoodFunction
  //! A "user-defined" likelihood function for testing
  //! 
  //! By default, the "data" is the combination of a single multidimensional
  //! Gaussian distribution in 10 dimensions, with center at 0.5 and variance 
  //! of 0.03 in each dimension
  //! 
  //! The default model may be changed at construction.
  //! 
  //! The variance may be modeled or fixed (using the SetDim member).
  //! 
  //! @ingroup likefunc
  class GaussTestMultiD : public LikelihoodFunction, public FunctDV
  {

  private:

    vector < vector<double> > fdata;
    unsigned ncomp, dim, mdim, N, levels;

    vector<double> weights;
    vector< vector<double> > centers, variance;

    UniformPtr unit;
    NormalPtr  norm;
    GammaPtr   gama;

    boost::shared_ptr<CDFGenerator> cdf;
    boost::shared_ptr<BetaRDist>    beta;
    
    bool varBeta, useBeta;
    double varBetaPDF(double, vector<double>&);
    double varBetaCDF(double, vector<double>&);

 public:

    /** Use the analytic CDF rather than the sampled CDF if true
	(default: false) */
    static bool useAnalytic; 

    //! Number of samples for CDF generation (default: 100000)
    static int cdfSamples; 

    //! Bucket size in KD tree for CDF generation (default: 16)
    static int ncut; 

    //+ CLICONSTR
    /** Default constructor

	Model with 1000 points by default, with variable centers, 10
	dimenions by default, one component (i.e. not a mixture)
    */
    GaussTestMultiD();

    //+ CLICONSTR int int
    /** Constructor with 
	@param N0 number of points
	@param Levels partitions.  

	If Levels is equal to smaller than 1, a single parition is
	used.

	As in the default model but with a user-specified number of
	points, user-specified number of levels.
     */
    GaussTestMultiD(int N0, int Levels);

    //+ CLICONSTR int int clivectord* clivectord*
    //+ CLICONSTR int int clivectord* clivectord* double
    /** Constructor with
	@param N0 number of points
	@param Levels partitions.  

	If Levels is equal to smaller than 1, a single parition is
	used.  The input vectors

	@param cen0 are the component centers
	@param var0 are the variance values

	User-specified number of points and levels, user-specified
	centers, user-specified variance.  The rank of these vectors
	(which must agree) implicitly specifies the dimension.

	@param beta is the shape of the beta function model to be used 
	(BetaRDist) instead of a Normal if beta>=0.0
     */
    GaussTestMultiD(int N0, int Levels, 
		    clivectord* cen0, clivectord* var0, double beta=-1.0);

    //+ CLICONSTR string int clivectord* clivectord*
    //+ CLICONSTR string int clivectord* clivectord* double
    /** Constructor with
	@param file is the data file name
	@param Levels partitions.  

	If Levels is equal to smaller than 1, a single parition is
	used.  The input vectors

	@param cen0 are the component centers
	@param var0 are the variance values

	User-specified data file and levels, user-specified
	centers, user-specified variance.  The rank of these vectors
	(which must agree) implicitly specifies the dimension.

	@param beta is the shape of the beta function model to be used 
	(BetaRDist) instead of a Normal if beta>=0.0
     */
    GaussTestMultiD(string file, int Levels, 
		    clivectord* cen0, clivectord* var0, double beta=-1.0);

    //+ CLICONSTR int int clivectord* clivectord* clivectord*
    //+ CLICONSTR int int clivectord* clivectord* clivectord* double
    /** Constructor with
	@param N0 number of points
	@param Levels partitions.  

	If Levels is equal to smaller than 1, a single parition is
	used.  The input vectors

	@param wght is a multicomponent weight vector
	@param cen0 are the component centers
	@param var0 are the variance values
	@param beta if beta>=0.0, we use a BetaRDist model rather 
	than Normal model


	User-specified number of points and levels, multiple
	components with user-specified weights, user-specified
	centers, user-specified variance.  The rank of these vectors
	(which must agree) implicitly specifies the dimension.  The
	rank of these vectors (which must agree) implicitly specifies
	the dimension.  If dim(wght) > 1, the vectors are assumed to
	be appended serially.
     */
    GaussTestMultiD(int N0, int Levels, clivectord* wght,
		    clivectord* cen0, clivectord* var0, double beta=-1.0);


    //+ CLICONSTR string int clivectord* clivectord* clivectord*
    //+ CLICONSTR string int clivectord* clivectord* clivectord* double
    /** Constructor with
	@param file contains the data set
	@param Levels partitions.  

	If Levels is equal to smaller than 1, a single parition is
	used.  The input vectors

	@param wght is a multicomponent weight vector
	@param cen0 are the component centers
	@param var0 are the variance values
	@param beta if beta>=0.0, we use a BetaRDist model rather 
	than Normal model


	User-specified data file and levels, multiple components with
	user-specified weights, user-specified centers, user-specified
	variance.  The rank of these vectors (which must agree)
	implicitly specifies the dimension.  The rank of these vectors
	(which must agree) implicitly specifies the dimension.  If
	dim(wght) > 1, the vectors are assumed to be appended
	serially.
     */
    GaussTestMultiD(string file, int Levels, clivectord* wght,
		    clivectord* cen0, clivectord* var0, double beta=-1.0);


    //+ CLIMETHOD void SetDim int
    //! Set the model dimension (currently either 1 or 2)
    void SetDim(int n);
    
    //+ CLIMETHOD void dumpData string
    //! Dump the entire data set to a file
    void dumpData(string s);

    //+ CLIMETHOD void readData string
    //! Read the entire data set from a file
    void readData(string s);

    //+ CLIMETHOD void printData
    //! Print data to a file for the current level
    void printData();

    /** This may be called after instantiation to create a model
	different from the synthetic data (for testing).  Obviously,
	this will have no effect if the variance is a parameter (model
	dimension=2)
    */
    //+ CLIMETHOD void newModel clivectord*
    void newModel(clivectord* var0);

    //! Use analytic rather than sampled CDF
    //+ CLIMETHOD void Analytic
    void Analytic() { useAnalytic = true; }

    /** Use a Gaussian model on Beta data.  This assumes that you have
	already instantiated using alpha>0 */
    //+ CLIMETHOD void Gaussian
    void Gaussian();

    /** Use a Beta model on Gaussian data.  This assumes that you have
	already instantiated using a Gaussian distribution (alpha<=0,
	the default constructor */
    //+ CLIMETHOD void BetaR double
    void BetaR(double alpha);

    /** Use a Beta model on Gaussian data with variable shape.  The
	range of the shape exponent is controlled by the prior
	distribution. This also assumes that you have already
	instantiated using a Gaussian distribution (alpha<=0, the
	default constructor */
    //+ CLIMETHOD void VarBetaR
    void VarBetaR() { 
      useBeta     = true;	//< Use beta-R distribution
      varBeta     = true;	//< Assume that shape (alpha) varies
      useAnalytic = true;	//< Analytic CDF rather than sampled
    }

    //! Reset the number of samples for the CDF generation (default: 100000)
    //+ CLIMETHOD void setSamples int int
    void setSamples(int val, int cut) { cdfSamples = val; ncut = cut; }

    //! This is likelihood function
    double LikeProb(std::vector<double> &z, SampleDistribution* sd, 
		    double norm, Tile *t, State *s, int indx)
    {
      if (s->Type()==StateInfo::None)
	return LikeProbSingle(s);
      else
	return LikeProbMixture(s);
    }


    //! The joint cumulative function
    double CumuProb(std::vector<double> &z, SampleDistribution* sd, 
		    double norm, Tile *t, State *s, int indx, Fct1dPtr f)
    {
      if (s->Type()==StateInfo::None)
	return CumuProb_single(s, f);
      else
	return CumuProb_mixture(s, f);
    }

    //! Label parameters.  Scheduled for removal.
    const std::string ParameterDescription(int i);

 protected:

    /// @name Likelihood functions
    //@{
    double LikeProbSingle(State *s);
    double LikeProbMixture(State *s);
    //@}

    /// @name Joint cumulative function
    //@{
    double CumuProb_single(State *s, Fct1dPtr f);
    double CumuProb_mixture(State *s, Fct1dPtr f);
    //@}

    //! Make the synthetic data
    void makeSyntheticData();

    //! Make the sampled CDF
    void makeCDF();

    //@{
    //! The samplers
    PairDV Sampler();
    vector<double> makeData();
    //! This is only used for initializing the CDF
    PairDV dd;
    //@}

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
          ar & BOOST_SERIALIZATION_NVP(levels);                        
          BIE_CATCH_BOOST_SERIALIZATION_EXCEPTION;                     
         }                                                             
         try {                                                         
          ar & BOOST_SERIALIZATION_NVP(weights);                        
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
          ar & BOOST_SERIALIZATION_NVP(unit);                        
          BIE_CATCH_BOOST_SERIALIZATION_EXCEPTION;                     
         }                                                             
         try {                                                         
          ar & BOOST_SERIALIZATION_NVP(norm);                        
          BIE_CATCH_BOOST_SERIALIZATION_EXCEPTION;                     
         }                                                             
         try {                                                         
          ar & BOOST_SERIALIZATION_NVP(gama);                        
          BIE_CATCH_BOOST_SERIALIZATION_EXCEPTION;                     
         }                                                             
         try {                                                         
          ar & BOOST_SERIALIZATION_NVP(cdf);                        
          BIE_CATCH_BOOST_SERIALIZATION_EXCEPTION;                     
         }                                                             
         try {                                                         
          ar & BOOST_SERIALIZATION_NVP(beta);                        
          BIE_CATCH_BOOST_SERIALIZATION_EXCEPTION;                     
         }                                                             
         try {                                                         
          ar & BOOST_SERIALIZATION_NVP(varBeta);                        
          BIE_CATCH_BOOST_SERIALIZATION_EXCEPTION;                     
         }                                                             
         try {                                                         
          ar & BOOST_SERIALIZATION_NVP(useBeta);                        
          BIE_CATCH_BOOST_SERIALIZATION_EXCEPTION;                     
         }                                                             
         try {                                                         
          ar & BOOST_SERIALIZATION_NVP(useAnalytic);                        
          BIE_CATCH_BOOST_SERIALIZATION_EXCEPTION;                     
         }                                                             
         try {                                                         
          ar & BOOST_SERIALIZATION_NVP(cdfSamples);                        
          BIE_CATCH_BOOST_SERIALIZATION_EXCEPTION;                     
         }                                                             
         try {                                                         
          ar & BOOST_SERIALIZATION_NVP(ncut);                        
          BIE_CATCH_BOOST_SERIALIZATION_EXCEPTION;                     
         }                                                             
        this->post_serialize(ar, version);
    }
    #endif

  };

} // namespace BIE

#ifndef SWIG
BIE_CLASS_TYPE_INFO(BIE::GaussTestMultiD)
BIE_CLASS_EXPORT_KEY(BIE::GaussTestMultiD)
#endif
#endif
