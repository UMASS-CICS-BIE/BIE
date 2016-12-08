// This is really -*- C++ -*-

#ifndef GaussTestLikelihoodFunctionRJ_h
#define GaussTestLikelihoodFunctionRJ_h

#include <LikelihoodFunction.h>
#include <gaussQ.h>

#include "Serializable.h"


namespace BIE {
  
  //+ CLICLASS GaussTestLikelihoodFunctionRJ SUPER LikelihoodFunction
  //! A "user-defined" likelihood function for testing
  //! 
  //! The "data" is the combination of two one-dimensional Gaussians,
  //! one at 0.2 with variance of 0.03 and one at 0.9 with variance of 0.03
  //! with 50% weights each.
  //! 
  //! The variance may be modeled or fixed (using the SetDim member).
  //! One may choose either point or binned data.  The default is point.
  //!
  //! @ingroup likefunc
  class GaussTestLikelihoodFunctionRJ : public LikelihoodFunction 
  {
    
  private:
    
    double         xmin, xmax, dx;
    vector<double> fdata, pdata;
    vector<double> centers, variance, weights;
    bool           point, cauchyData, cauchyModel;

    LegeQPtr lq;
    uint32   number;

    int  dim;
    int  nbins;
    int  N;
    
  public:

    //! Number of knots for bin average (default: 6)
    static int nint;

    //+ CLICONSTR
    //! Null constructor
    GaussTestLikelihoodFunctionRJ();
    
    //+ CLICONSTR int int
    //! Constructor with number of bins and number of sample points
    GaussTestLikelihoodFunctionRJ(int Nbins, int N0);
    
    //+ CLICONSTR int
    //! Constructor with number of sample points for point model
    GaussTestLikelihoodFunctionRJ(int N0);
    
    //! Destructor
    virtual ~GaussTestLikelihoodFunctionRJ() {}
    
    //+ CLIMETHOD void useCauchyModel bool
    //! Set the likelihood model to a Cauchy distribution
    void useCauchyModel(bool b) { cauchyModel = b; }

    //+ CLIMETHOD void useCauchyData bool
    //! Generate the data from a Cauchy distribution
    void useCauchyData(bool b) { cauchyData = b; }

    //+ CLIMETHOD void newModel string
    //! Read the model from a file with name @param file
    void newModel(string file);
    
    //+ CLIMETHOD void SetDim int
    //! Set the model dimension (currently either 1 or 2)
    void SetDim(int n);
    
    //! This is likelihood function
    double LikeProb(std::vector<double> &z, SampleDistribution* sd, 
		    double norm, Tile *t, State *s, int indx);
    
    //! The joint cumulative function
    double CumuProb(std::vector<double> &z, SampleDistribution* sd, 
		    double norm, Tile *t, State *s, int indx, Fct1dPtr f);

    //! Label parameters.  Scheduled for removal.
    const std::string ParameterDescription(int i);
    
  protected:
    
    //! Likelihood function
    double LikeProbTwo(State *s);
    
    //! The joint cumulative function
    double CumuProbTwo(State *s, Fct1dPtr f);

    //! Make a model with default parameters
    void defaultModel();
    
    //! Make the synthetic binned data
    void makeSyntheticData();
    
    //! Make the synthetic point data
    void makeSyntheticPointData();
    
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
          ar & BOOST_SERIALIZATION_NVP(xmin);                        
          BIE_CATCH_BOOST_SERIALIZATION_EXCEPTION;                     
         }                                                             
         try {                                                         
          ar & BOOST_SERIALIZATION_NVP(xmax);                        
          BIE_CATCH_BOOST_SERIALIZATION_EXCEPTION;                     
         }                                                             
         try {                                                         
          ar & BOOST_SERIALIZATION_NVP(dx);                        
          BIE_CATCH_BOOST_SERIALIZATION_EXCEPTION;                     
         }                                                             
         try {                                                         
          ar & BOOST_SERIALIZATION_NVP(fdata);                        
          BIE_CATCH_BOOST_SERIALIZATION_EXCEPTION;                     
         }                                                             
         try {                                                         
          ar & BOOST_SERIALIZATION_NVP(pdata);                        
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
          ar & BOOST_SERIALIZATION_NVP(weights);                        
          BIE_CATCH_BOOST_SERIALIZATION_EXCEPTION;                     
         }                                                             
         try {                                                         
          ar & BOOST_SERIALIZATION_NVP(point);                        
          BIE_CATCH_BOOST_SERIALIZATION_EXCEPTION;                     
         }                                                             
         try {                                                         
          ar & BOOST_SERIALIZATION_NVP(cauchyData);                        
          BIE_CATCH_BOOST_SERIALIZATION_EXCEPTION;                     
         }                                                             
         try {                                                         
          ar & BOOST_SERIALIZATION_NVP(cauchyModel);                        
          BIE_CATCH_BOOST_SERIALIZATION_EXCEPTION;                     
         }                                                             
         try {                                                         
          ar & BOOST_SERIALIZATION_NVP(lq);                        
          BIE_CATCH_BOOST_SERIALIZATION_EXCEPTION;                     
         }                                                             
         try {                                                         
          ar & BOOST_SERIALIZATION_NVP(number);                        
          BIE_CATCH_BOOST_SERIALIZATION_EXCEPTION;                     
         }                                                             
         try {                                                         
          ar & BOOST_SERIALIZATION_NVP(dim);                        
          BIE_CATCH_BOOST_SERIALIZATION_EXCEPTION;                     
         }                                                             
         try {                                                         
          ar & BOOST_SERIALIZATION_NVP(nbins);                        
          BIE_CATCH_BOOST_SERIALIZATION_EXCEPTION;                     
         }                                                             
         try {                                                         
          ar & BOOST_SERIALIZATION_NVP(N);                        
          BIE_CATCH_BOOST_SERIALIZATION_EXCEPTION;                     
         }                                                             
        this->post_serialize(ar, version);
    }
    #endif

  };

} //namespace BIE

#ifndef SWIG
BIE_CLASS_TYPE_INFO(BIE::GaussTestLikelihoodFunctionRJ)
BIE_CLASS_EXPORT_KEY(BIE::GaussTestLikelihoodFunctionRJ)
#endif
#endif
