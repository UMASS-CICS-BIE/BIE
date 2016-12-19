// This is really -*- C++ -*-

#ifndef EnsembleStat_h
#define EnsembleStat_h

#include <Ensemble.h>

#include "Serializable.h"


using namespace std;

namespace BIE {
  
  //+ CLICLASS EnsembleStat SUPER Ensemble
  //! Caches an group of posterior states in a vector and provides
  //! member functions to characterize the posterior from this
  //! distribution and produce statistical diagnostics on the
  //! convergence of the simulation.
  //! 
  //! The posterior state is characterized from the ensemble variance
  //! by a principal component analysis.
  //!
  //! @ingroup ensembles
  //!
  class EnsembleStat : public Ensemble
  {
  protected:
    //! Trimmed dimension for each subspace (that is, remove the null space)
    map<int, int> icut;

    //! Trimmed dimension for each subspace for each component
    map<int, vector<int> > mcut;

    //! Fraction of states in each subspace per component
    map<int, int> mcnt;

    //! Fraction of states in each subspace per component
    map<int, double> mcum;

    //! True of mean and variance has been computed
    bool var_computed;

    //! Mean for each dimension
    //! @{
    map< int, vector<double> > vmean;
    map< int, vector< vector<double> > > mmean;
    //! @}

    //! Variance for each dimension
    //! @{
    map< int, vector<double> > vvar;
    map< int, vector< vector<double> > > mvar;
    //! @}
    
    //! @{
    //! Lower and upper limits
    vector<double> vlower;
    vector<double> vupper;
    //! @}

    //! Normal distribution generator
    boost::shared_ptr<Normal> normal;

    /** Scaling factor for realizing variates from the
	multidimensional Gaussian defined by the mean and covariance
	of the ensemble */
    double factor;

    //! Accumulate subspace counts
    virtual void addCount(int m);

  public:

    //! Minimum number of states per subspace to retain subspace (default 3)
    static int minsub;

    //! Constructor
    EnsembleStat();

    //+ CLICONSTR StateInfo*
    /** Mixture constructor
	\param si is the state metadata structure
    */
    EnsembleStat(StateInfo *si);

    //+ CLICONSTR StateInfo* int int string int
    /** General constructor
	\param si is the state metadata structure
	\param level is the update level
	\param nburn is the desired number of converged states
	\param filename is the logfile for ensemble statistics
	\param keypos is the parameter index for computing PDF and CDF 
	evaluations 
    */
    EnsembleStat(StateInfo *si,
		 int level, int nburn, string filename, int keypos);

    //! Copy constructor
    EnsembleStat(const EnsembleStat&);

    //! Copy constructor for the CLI
    //+ CLICONSTR Ensemble*
    EnsembleStat(Ensemble*);

    //! Destructor
    virtual ~EnsembleStat() {}
    
    //+ CLIMETHOD void Reset StateInfo*
    //! Reinitialize the EnsembleStat
    void Reset(StateInfo *si);

    //+ CLIMETHOD void Reset StateInfo* int int string int
    //! Reinitialize the EnsembleStat
    void Reset(StateInfo *si,
	       int level, int nburn, string filename, int keypos);

    //+ CLIMETHOD void setDimensions StateInfo*
    //! Set total size of arrays from the metadata structure
    void setDimensions(StateInfo *si);


    /**
       The covariance matrix for all parameters and weights in the
       deque.  This and Teigen and eigen are only defined after
       EnsembleStat::ComputeDistribution() is called.  
    */
    map<int, MatrixM> covar;

    /**
       The covariance matrix for all parameters and weights in the
       deque by component.  This and TeigenM and eigenM are only
       defined after EnsembleStat::ComputeDistributionMarginal() is
       called.
    */
    map<int, vector<MatrixM> > covarM;

    //! The transpose of the eigenvector matrix
    map<int, MatrixM> Teigen;

    //! The transpose of the eigenvector matrix for the components
    map<int, vector<MatrixM> > TeigeM;

    /** The eigenvector matrix.  Eigenvectors are ranked in decreasing
	order by their eigenvalues.
    */
    //@{
    //! The eigenvectors for the subspace
    map<int, MatrixM> eigen;

    //! The eigenvectors for each component in the subspace
    map<int, vector<MatrixM> > eigenM;

    //! The mean for each dimension
    map<int, VectorM> mean;

    //! The mean for each dimension for each component
    map<int, vector<VectorM> > meanM;

    //! The eigenvalues
    map<int, VectorM> evals;

    //! The eigenvalues by component
    map<int, vector<VectorM> > evalsM;

    //+ CLIMETHOD void ComputeDistribution
    //+ CLIMETHOD void ComputeDistribution int
    /** 
	@name ComputeDistribution
	Computes the covariance matrix and its eigenvalues/vectors
	
	Return = 1 (variance computed), Return = 0 (too few values)
    */
    void ComputeDistribution();
    void ComputeDistribution(unsigned n) { burnIn=n; ComputeDistribution(); }
    void ComputeDistributionMarginal(ostream& out);
    
    //! Index of non-null space for subspace m
    int trim_location(unsigned m) { return icut[m]; }

    /** Return a single value of the marginal distribution where 
	@param m is the index of non-null space for subspace 
	@param n is the desired component  index
    */
    int trim_location_marginal(unsigned m, unsigned n) { return mcut[m][n]; }
    
    //! Object factory (clone)
    EnsembleStat* New() {
      return new EnsembleStat(_si);
    }

    //! Differential distribution function P(x)
    double PDF(State&);

    //! Log of differential distribution function P(x)
    double logPDF(State&);

    /** Log of differential distribution function where
	@param m is the mixture component subspace
	@param n is the component index
	@param V is parameter vector
    */
    double logPDFMarginal(unsigned m, unsigned n, const vector<double>& V);

    //! Differential distribution function for a single component
    double PDF(unsigned, State&);

    //! Lower bound on distribution (in each dimension)
    vector<double> lower(void);

    //! Upper bound on distribution (in each dimension)
    vector<double> upper(void);
    
    //! Return mean of distribution (mulitvariate)
    vector<double> Mean(unsigned m);

    //! Return standard deviation of distribution (mulitvariate)
    vector<double> StdDev(unsigned m);

    //! Return specifided moment of distribution (mulitvariate)
    vector<double> Moments(unsigned m, unsigned k);
    
    //! Return random variate from distribution for the given subspace m
    virtual State Sample(unsigned m);
    
    //! Return random variate from fully marginalized distribution
    double SampleOne(unsigned m, unsigned j);
    
    //! Retrieve the std dev for the component marginal
    vector<double> StdDevMarginal(unsigned m, unsigned n);

    //! Retrieve the mean for the component marginal
    vector<double> MeanMarginal(unsigned m, unsigned n);

    //! Sample the distribution for the component marginal
    vector<double> SampleMarginal(unsigned m, unsigned n);

    //! Width for component marginal
    vector<double> WidthMarginal(unsigned m, unsigned n) 
    { return StdDevMarginal(m, n); }

    //! Print current covariance matrix (debug and diag)
    //+ CLIMETHOD void PrintDiag
    //+ CLIMETHOD void PrintDiag string
    //@{
    //! To a given ostream
    virtual void PrintDiag(ostream& out);
    //! To a given file
    virtual void PrintDiag(string& outfile);
    //! Print to console
    virtual void PrintDiag();

    //+ CLIMETHOD void PrintDensity int int int int string
    //! Density image (for debugging)
    void PrintDensity(int dim1, int dim2, int num1, int num2,
		      string file);

    #ifndef SWIG
    // AUTO GENERATED BY ../persistence/autopersist.py
    private:
    friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive & ar, const unsigned int version) {
        this->pre_serialize(ar, version);
         try {                                                         
          ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(Ensemble);            
          BIE_CATCH_BOOST_SERIALIZATION_EXCEPTION;                     
         }                                                             
         try {                                                         
          ar & BOOST_SERIALIZATION_NVP(icut);                        
          BIE_CATCH_BOOST_SERIALIZATION_EXCEPTION;                     
         }                                                             
         try {                                                         
          ar & BOOST_SERIALIZATION_NVP(mcut);                        
          BIE_CATCH_BOOST_SERIALIZATION_EXCEPTION;                     
         }                                                             
         try {                                                         
          ar & BOOST_SERIALIZATION_NVP(mcnt);                        
          BIE_CATCH_BOOST_SERIALIZATION_EXCEPTION;                     
         }                                                             
         try {                                                         
          ar & BOOST_SERIALIZATION_NVP(mcum);                        
          BIE_CATCH_BOOST_SERIALIZATION_EXCEPTION;                     
         }                                                             
         try {                                                         
          ar & BOOST_SERIALIZATION_NVP(var_computed);                        
          BIE_CATCH_BOOST_SERIALIZATION_EXCEPTION;                     
         }                                                             
         try {                                                         
          ar & BOOST_SERIALIZATION_NVP(vmean);                        
          BIE_CATCH_BOOST_SERIALIZATION_EXCEPTION;                     
         }                                                             
         try {                                                         
          ar & BOOST_SERIALIZATION_NVP(mmean);                        
          BIE_CATCH_BOOST_SERIALIZATION_EXCEPTION;                     
         }                                                             
         try {                                                         
          ar & BOOST_SERIALIZATION_NVP(vvar);                        
          BIE_CATCH_BOOST_SERIALIZATION_EXCEPTION;                     
         }                                                             
         try {                                                         
          ar & BOOST_SERIALIZATION_NVP(mvar);                        
          BIE_CATCH_BOOST_SERIALIZATION_EXCEPTION;                     
         }                                                             
         try {                                                         
          ar & BOOST_SERIALIZATION_NVP(vlower);                        
          BIE_CATCH_BOOST_SERIALIZATION_EXCEPTION;                     
         }                                                             
         try {                                                         
          ar & BOOST_SERIALIZATION_NVP(vupper);                        
          BIE_CATCH_BOOST_SERIALIZATION_EXCEPTION;                     
         }                                                             
         try {                                                         
          ar & BOOST_SERIALIZATION_NVP(normal);                        
          BIE_CATCH_BOOST_SERIALIZATION_EXCEPTION;                     
         }                                                             
         try {                                                         
          ar & BOOST_SERIALIZATION_NVP(factor);                        
          BIE_CATCH_BOOST_SERIALIZATION_EXCEPTION;                     
         }                                                             
         try {                                                         
          ar & BOOST_SERIALIZATION_NVP(minsub);                        
          BIE_CATCH_BOOST_SERIALIZATION_EXCEPTION;                     
         }                                                             
         try {                                                         
          ar & BOOST_SERIALIZATION_NVP(covar);                        
          BIE_CATCH_BOOST_SERIALIZATION_EXCEPTION;                     
         }                                                             
         try {                                                         
          ar & BOOST_SERIALIZATION_NVP(covarM);                        
          BIE_CATCH_BOOST_SERIALIZATION_EXCEPTION;                     
         }                                                             
         try {                                                         
          ar & BOOST_SERIALIZATION_NVP(Teigen);                        
          BIE_CATCH_BOOST_SERIALIZATION_EXCEPTION;                     
         }                                                             
         try {                                                         
          ar & BOOST_SERIALIZATION_NVP(TeigeM);                        
          BIE_CATCH_BOOST_SERIALIZATION_EXCEPTION;                     
         }                                                             
         try {                                                         
          ar & BOOST_SERIALIZATION_NVP(eigen);                        
          BIE_CATCH_BOOST_SERIALIZATION_EXCEPTION;                     
         }                                                             
         try {                                                         
          ar & BOOST_SERIALIZATION_NVP(eigenM);                        
          BIE_CATCH_BOOST_SERIALIZATION_EXCEPTION;                     
         }                                                             
         try {                                                         
          ar & BOOST_SERIALIZATION_NVP(mean);                        
          BIE_CATCH_BOOST_SERIALIZATION_EXCEPTION;                     
         }                                                             
         try {                                                         
          ar & BOOST_SERIALIZATION_NVP(meanM);                        
          BIE_CATCH_BOOST_SERIALIZATION_EXCEPTION;                     
         }                                                             
         try {                                                         
          ar & BOOST_SERIALIZATION_NVP(evals);                        
          BIE_CATCH_BOOST_SERIALIZATION_EXCEPTION;                     
         }                                                             
         try {                                                         
          ar & BOOST_SERIALIZATION_NVP(evalsM);                        
          BIE_CATCH_BOOST_SERIALIZATION_EXCEPTION;                     
         }                                                             
        this->post_serialize(ar, version);
    }
    #endif

  };
  #ifndef SWIG
  BOOST_SERIALIZATION_SHARED_PTR(EnsembleStat)
  #endif
  //! @}

} // namespace BIE

#ifndef SWIG
BIE_CLASS_TYPE_INFO(BIE::EnsembleStat)
BIE_CLASS_EXPORT_KEY(BIE::EnsembleStat)
#endif
#endif
  
  
  
  
  
  
  
  
