// This is really -*- C++ -*-

#ifndef EnsembleKD_h
#define EnsembleKD_h

#include <map>

#include <KDTree.h>
#include <Ensemble.h>

#include "Serializable.h"


using namespace std;

namespace BIE {
  
  //+ CLICLASS EnsembleKD SUPER Ensemble
  //! Caches an group of posterior states in a vector and provides
  //! member functions to characterize the posterior from this
  //! distribution and produce statistical diagnostics on the
  //! convergence of the simulation.
  //!
  //! The posterior state is sampled from the input chain and
  //! characterized by a KD tessellation
  //!
  //! @ingroup ensembles
  //!
  class EnsembleKD : public Ensemble
  {
  protected:

    //! True of the distribution has already been computed
    bool dist_computed;
    
    //! Number of states in each subspace
    map<int, int> mcnt;

    //! Fraction of states in each subspace
    map<int, double> mcum;

    //! Probability scaling
    vector<double> dataScale;

    //! Holds the tree for each mixture cardinality
    //! @{
    map<int, boost::shared_ptr<KDTree> > kde;
    map<int, vector<boost::shared_ptr<KDTree> > > kdeM;
    //! @}

    //! The mean for each component
    //! @{
    map<int, vector<double> > vmean;
    map<int, vector< vector<double> > > mmean;
    //! @}

    //! The peak for each component (used for density plotting)
    map<int, vector<double> > peak;

    //! The variance for each component
    //! @{
    map< int, vector<double> > vvar;
    map< int, vector< vector<double> > > mvar;
    //! @}
    
    //! The mean for each mixture subspace
    //! @{
    map< int, VectorM > mean;
    map< int, vector<VectorM> > meanM;
    //! @}

    //! The covariance for each mixture subspace
    //! @{
    map< int, MatrixM > covar;
    map< int, vector<MatrixM> > covarM;
    //! @}

    //! Cache of all states used for sampling
    boost::shared_ptr<StateCache> cache;

    //! The points used in the tessellation
    //! @{
    map<int, vector< vector<double> > > points;
    map<int, vector< vector< vector<double> > > > pointsM;
    //! @}

    //! The computed weight for each point
    //! @{
    map<int, vector<double> > weight;
    map<int, vector< vector<double> > > weightM;
    //! @}

    //! The compute lower bound
    std::vector<double> vlower;

    //! The compute upper bound
    std::vector<double> vupper;

    //! Accumulate subspace counts
    virtual void addCount(int m);

  public:

    //! Target bucket size for KD-tree (default: 32)
    static int ncut;

    //! Minimum number of states per subspace to retain subspace (default: 3)
    static int minsub;

    //+ CLICONSTR
    //! Default constructor
    EnsembleKD();

    //+ CLICONSTR StateInfo*
    //! Constructor from StateInfo
    EnsembleKD(StateInfo *si);

    //+ CLICONSTR StateInfo* int int string int
    /** General constructor
	\param si is the state metadata instance
	\param level is the update level
	\param nburn is the desired number of converged states
	\param filename is the logfile for ensemble statistics
	\param keypos is the parameter index for computing PDF and CDF 
	evaluations 
    */
    EnsembleKD(StateInfo* si,
	       int level, int nburn, string filename, int keypos);

    //! Copy constructor
    EnsembleKD(const EnsembleKD&);

    //! Copy constructor (base class instance for CLI)
    //+ CLICONSTR Ensemble*
    EnsembleKD(Ensemble*);

    //! Destructor
    virtual ~EnsembleKD() {}
    
    //+ CLIMETHOD void Reset StateInfo*
    //! Reinitialize the EnsembleKD
    void Reset(StateInfo*);

    //+ CLIMETHOD void Reset StateInfo* int int string int
    //! Reinitialize the EnsembleKD
    void Reset(StateInfo* si,
	       int level, int nburn, string filename, int keypos);

    //+ CLIMETHOD void setBucketSize int
    //! Set arget bucket size for KD-tree
    void setBucketSize (int n) { ncut = n; }

    //+ CLIMETHOD void ComputeDistribution
    //+ CLIMETHOD void ComputeDistribution int
    /** 
	@name ComputeDistribution
	Computes the covariance matrix and its eigenvalues/vectors
    */
    void ComputeDistribution();
    void ComputeDistribution(int n) { burnIn=n; ComputeDistribution(); }
    void ComputeDistributionMarginal(ostream& out);
    
    //! Object factory (clone)
    EnsembleKD* New() {
      return new EnsembleKD(_si);
    }

    //! Get dimension of state vector from cache
    unsigned getDim() { return cache->stateSize(); }

    //! Differential distribution function P(x)
    double PDF(State&);

    //! Log of differential distribution function P(x)
    double logPDF(State&);

    /** 
	Log of differential distribution function for a component in
	@param m subspace 
	@param n component
	@param V parameter vector
    */
    double logPDFMarginal(unsigned m, unsigned n, vector<double>& V);

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
    State Sample(unsigned m);
    
    /**
       Return random variate for a single component from distribution
       for the given 

       @param m subspace 
       @param n component
    */
    vector<double> SampleMarginal(unsigned m, unsigned n);

    //! Return random variate from fully marginalized distribution
    double SampleOne(unsigned m, unsigned j);
    
    //! Print current covariance matrix (debug and diag)
    //@{
    //! Use a supplied ostream for diag info
    void PrintDiag(ostream& out);

    //! Open a new file for diag into
    void PrintDiag(string& outfile);

    //! Print to console
    void PrintDiag();
    //@}

    //+ CLIMETHOD void PrintDensity int int int int string
    //! Density image (for debugging)
    void PrintDensity(int dim1, int dim2, int num1, int num2,
		      string file);

  private:
    //! Hold over from previous implementation that did not persist the tree
    void RestoreDistribution();

  private:
    //
    // Work around for problem with map restore.  Not sure what is
    // going here.
    //
    template<class Archive>
      void post_save(Archive &ar, const unsigned int file_version) const
    {
      SAVEDMDV(points,  vector<vector<double> >         );
      SAVEDMDV(pointsM, vector<vector<vector<double> > >);
      SAVEDMDV(weight,  vector<double>                  );
      SAVEDMDV(weightM, vector<vector<double> >         );
    }

    template<class Archive>
    void post_load(Archive &ar, const unsigned int file_version) 
    {
      RESTRMDV(points,  vector<vector<double> >         );
      RESTRMDV(pointsM, vector<vector<vector<double> > >);
      RESTRMDV(weight,  vector<double>                  );
      RESTRMDV(weightM, vector<vector<double> >         );
    }
    
    #ifndef SWIG
    // AUTO GENERATED BY ../persistence/autopersist.py
    private:
    friend class boost::serialization::access;
    BOOST_SERIALIZATION_SPLIT_MEMBER();

    template<class Archive>
    void save(Archive & ar, const unsigned int version) const {
        this->pre_save(ar, version);
         try {                                                         
          ar << BOOST_SERIALIZATION_BASE_OBJECT_NVP(Ensemble);            
          BIE_CATCH_BOOST_SERIALIZATION_EXCEPTION;                     
         }                                                             
         try {                                                         
          ar << BOOST_SERIALIZATION_NVP(dist_computed);                        
          BIE_CATCH_BOOST_SERIALIZATION_EXCEPTION;                     
         }                                                             
         try {                                                         
          ar << BOOST_SERIALIZATION_NVP(mcnt);                        
          BIE_CATCH_BOOST_SERIALIZATION_EXCEPTION;                     
         }                                                             
         try {                                                         
          ar << BOOST_SERIALIZATION_NVP(mcum);                        
          BIE_CATCH_BOOST_SERIALIZATION_EXCEPTION;                     
         }                                                             
         try {                                                         
          ar << BOOST_SERIALIZATION_NVP(dataScale);                        
          BIE_CATCH_BOOST_SERIALIZATION_EXCEPTION;                     
         }                                                             
         try {                                                         
          ar << BOOST_SERIALIZATION_NVP(kde);                        
          BIE_CATCH_BOOST_SERIALIZATION_EXCEPTION;                     
         }                                                             
         try {                                                         
          ar << BOOST_SERIALIZATION_NVP(kdeM);                        
          BIE_CATCH_BOOST_SERIALIZATION_EXCEPTION;                     
         }                                                             
         try {                                                         
          ar << BOOST_SERIALIZATION_NVP(vmean);                        
          BIE_CATCH_BOOST_SERIALIZATION_EXCEPTION;                     
         }                                                             
         try {                                                         
          ar << BOOST_SERIALIZATION_NVP(mmean);                        
          BIE_CATCH_BOOST_SERIALIZATION_EXCEPTION;                     
         }                                                             
         try {                                                         
          ar << BOOST_SERIALIZATION_NVP(peak);                        
          BIE_CATCH_BOOST_SERIALIZATION_EXCEPTION;                     
         }                                                             
         try {                                                         
          ar << BOOST_SERIALIZATION_NVP(vvar);                        
          BIE_CATCH_BOOST_SERIALIZATION_EXCEPTION;                     
         }                                                             
         try {                                                         
          ar << BOOST_SERIALIZATION_NVP(mvar);                        
          BIE_CATCH_BOOST_SERIALIZATION_EXCEPTION;                     
         }                                                             
         try {                                                         
          ar << BOOST_SERIALIZATION_NVP(mean);                        
          BIE_CATCH_BOOST_SERIALIZATION_EXCEPTION;                     
         }                                                             
         try {                                                         
          ar << BOOST_SERIALIZATION_NVP(meanM);                        
          BIE_CATCH_BOOST_SERIALIZATION_EXCEPTION;                     
         }                                                             
         try {                                                         
          ar << BOOST_SERIALIZATION_NVP(covar);                        
          BIE_CATCH_BOOST_SERIALIZATION_EXCEPTION;                     
         }                                                             
         try {                                                         
          ar << BOOST_SERIALIZATION_NVP(covarM);                        
          BIE_CATCH_BOOST_SERIALIZATION_EXCEPTION;                     
         }                                                             
         try {                                                         
          ar << BOOST_SERIALIZATION_NVP(cache);                        
          BIE_CATCH_BOOST_SERIALIZATION_EXCEPTION;                     
         }                                                             
         try {                                                         
          ar << BOOST_SERIALIZATION_NVP(vlower);                        
          BIE_CATCH_BOOST_SERIALIZATION_EXCEPTION;                     
         }                                                             
         try {                                                         
          ar << BOOST_SERIALIZATION_NVP(vupper);                        
          BIE_CATCH_BOOST_SERIALIZATION_EXCEPTION;                     
         }                                                             
         try {                                                         
          ar << BOOST_SERIALIZATION_NVP(ncut);                        
          BIE_CATCH_BOOST_SERIALIZATION_EXCEPTION;                     
         }                                                             
         try {                                                         
          ar << BOOST_SERIALIZATION_NVP(minsub);                        
          BIE_CATCH_BOOST_SERIALIZATION_EXCEPTION;                     
         }                                                             
        this->post_save(ar, version);
    }

    template<class Archive>
    void load(Archive & ar, const unsigned int version) {
        this->pre_load(ar, version);
         try {                                                         
          ar >> BOOST_SERIALIZATION_BASE_OBJECT_NVP(Ensemble);            
          BIE_CATCH_BOOST_SERIALIZATION_EXCEPTION;                     
         }                                                             
         try {                                                         
          ar >> BOOST_SERIALIZATION_NVP(dist_computed);                        
          BIE_CATCH_BOOST_SERIALIZATION_EXCEPTION;                     
         }                                                             
         try {                                                         
          ar >> BOOST_SERIALIZATION_NVP(mcnt);                        
          BIE_CATCH_BOOST_SERIALIZATION_EXCEPTION;                     
         }                                                             
         try {                                                         
          ar >> BOOST_SERIALIZATION_NVP(mcum);                        
          BIE_CATCH_BOOST_SERIALIZATION_EXCEPTION;                     
         }                                                             
         try {                                                         
          ar >> BOOST_SERIALIZATION_NVP(dataScale);                        
          BIE_CATCH_BOOST_SERIALIZATION_EXCEPTION;                     
         }                                                             
         try {                                                         
          ar >> BOOST_SERIALIZATION_NVP(kde);                        
          BIE_CATCH_BOOST_SERIALIZATION_EXCEPTION;                     
         }                                                             
         try {                                                         
          ar >> BOOST_SERIALIZATION_NVP(kdeM);                        
          BIE_CATCH_BOOST_SERIALIZATION_EXCEPTION;                     
         }                                                             
         try {                                                         
          ar >> BOOST_SERIALIZATION_NVP(vmean);                        
          BIE_CATCH_BOOST_SERIALIZATION_EXCEPTION;                     
         }                                                             
         try {                                                         
          ar >> BOOST_SERIALIZATION_NVP(mmean);                        
          BIE_CATCH_BOOST_SERIALIZATION_EXCEPTION;                     
         }                                                             
         try {                                                         
          ar >> BOOST_SERIALIZATION_NVP(peak);                        
          BIE_CATCH_BOOST_SERIALIZATION_EXCEPTION;                     
         }                                                             
         try {                                                         
          ar >> BOOST_SERIALIZATION_NVP(vvar);                        
          BIE_CATCH_BOOST_SERIALIZATION_EXCEPTION;                     
         }                                                             
         try {                                                         
          ar >> BOOST_SERIALIZATION_NVP(mvar);                        
          BIE_CATCH_BOOST_SERIALIZATION_EXCEPTION;                     
         }                                                             
         try {                                                         
          ar >> BOOST_SERIALIZATION_NVP(mean);                        
          BIE_CATCH_BOOST_SERIALIZATION_EXCEPTION;                     
         }                                                             
         try {                                                         
          ar >> BOOST_SERIALIZATION_NVP(meanM);                        
          BIE_CATCH_BOOST_SERIALIZATION_EXCEPTION;                     
         }                                                             
         try {                                                         
          ar >> BOOST_SERIALIZATION_NVP(covar);                        
          BIE_CATCH_BOOST_SERIALIZATION_EXCEPTION;                     
         }                                                             
         try {                                                         
          ar >> BOOST_SERIALIZATION_NVP(covarM);                        
          BIE_CATCH_BOOST_SERIALIZATION_EXCEPTION;                     
         }                                                             
         try {                                                         
          ar >> BOOST_SERIALIZATION_NVP(cache);                        
          BIE_CATCH_BOOST_SERIALIZATION_EXCEPTION;                     
         }                                                             
         try {                                                         
          ar >> BOOST_SERIALIZATION_NVP(vlower);                        
          BIE_CATCH_BOOST_SERIALIZATION_EXCEPTION;                     
         }                                                             
         try {                                                         
          ar >> BOOST_SERIALIZATION_NVP(vupper);                        
          BIE_CATCH_BOOST_SERIALIZATION_EXCEPTION;                     
         }                                                             
         try {                                                         
          ar >> BOOST_SERIALIZATION_NVP(ncut);                        
          BIE_CATCH_BOOST_SERIALIZATION_EXCEPTION;                     
         }                                                             
         try {                                                         
          ar >> BOOST_SERIALIZATION_NVP(minsub);                        
          BIE_CATCH_BOOST_SERIALIZATION_EXCEPTION;                     
         }                                                             
        this->post_load(ar, version);
    }
    #endif

  };
  #ifndef SWIG
  BOOST_SERIALIZATION_SHARED_PTR(EnsembleKD)
  #endif
}
#ifndef SWIG
BIE_CLASS_TYPE_INFO(BIE::EnsembleKD)
BIE_CLASS_EXPORT_KEY(BIE::EnsembleKD)
#endif
#endif
