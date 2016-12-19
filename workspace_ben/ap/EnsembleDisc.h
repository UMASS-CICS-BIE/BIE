// This is really -*- C++ -*-

#ifndef EnsembleDisc_h
#define EnsembleDisc_h

#include <map>

#include <MetricTreeDensity.h>
#include <Ensemble.h>

#include "Serializable.h"


using namespace std;

namespace BIE {
  
  //+ CLICLASS EnsembleDisc SUPER Ensemble
  /** Caches a group of posterior states in a vector.

      Provides member functions to characterize the posterior from
      this distribution and produce statistical diagnostics on the
      convergence of the simulation.
  
      The posterior state is sampled from the input chain and
      characterized kernel density estimation
  
      @ingroup ensembles
  */
  class EnsembleDisc : public Ensemble
  {

  protected:

    //! Debug counter
    static int dbg_ctr;

    //! True of the distribution has already been computed
    bool dist_computed;
    
    //@{
    //! KDE has already been computed
    bool kde_computed;
    bool kdeM_computed;
    //@}
    
    //! Target number of points in metric tree bucket (default: 16)
    unsigned bucketSize;

    //! Number of states in each subspace
    map<int, int> mcnt;

    //! Fraction of states in each subspace
    map<int, double> mcum;

    //! @{
    //! Holds the kernel density estimation for each mixture cardinality
    map<int, MetricDensPtr > kde;
    map<int, vector<MetricDensPtr > > kdeM;
    //! @}

    //! @{
    //! The mean for each component
    map<int, vector<double> > vmean;
    map<int, vector< vector<double> > > mmean;
    //! @}

    //! The peak for each component (used for density plotting)
    map<int, vector<double> > peak;

    //! @{
    //! Lower and upper limits
    vector<double> vlower;
    vector<double> vupper;
    //! @}

    //! @{
    //! The variance for each component
    map< int,  vector<double> > vvar;
    map< int, vector< vector<double> > > mvar;
    //! @}
    
    //! @{
    //! The mean for each mixture subspace
    map<int, VectorM> mean;
    map<int, vector<VectorM> > meanM;
    //! @}


    //! @{
    /** Scale each dimension in the KDE distance metric computation by
	user specified scaling vector */
    std::vector<double> fscale;
    //! @}

    //! @{
    //! The covariance for each mixture subspace
    map<int, MatrixM> covar;
    map<int, vector<MatrixM> > covarM;
    //! @}

    //! Cache of all states used for sampling
    CachePtr cache;

    //! @{
    //! The points used in the density estimation
    // map<int, vector< vector<double> > > @autopersist(points);
    map<int, vector< vector<double> > > points;
    map<int, vector< vector< vector<double> > > > pointsM;
    //! @}

    //! @{
    //! The computed width for each point
    map<int, vector< vector<double> > > widths;
    map<int, vector< vector< vector<double> > > > widthsM;
    //! @}

    //! @{
    //! The computed weight for each point
    map<int, vector<double> > weight;
    map<int, vector< vector<double> > > weightM;
    //! @}
    
    //! The kernel type used for density estimation
    KernelPtr kernel;

    //! @{
    /// Compute density using the metric tree
    void ComputeDensity ();
    void ComputeDensityM();
    //! @}

    /** For the metric tree, sort the points and remove duplicates.
	Duplicates will break the ball tree construction algorithm.
    */
    //! @{
    void sortPoints(map<int, vector< vector<double> > >& points,
		    map<int, vector< vector<double> > >& widths,
		    map<int, vector<double> >&           weight );

    void sortPointsM(map<int, vector< vector< vector<double> > > >& points,
		     map<int, vector< vector< vector<double> > > >& widths,
		     map<int, vector< vector<double> > >&           weight );
    //! @}

    //! The maximum number of states to use in creating the density estimate
    //! (default: 0, which means all of the states)
    int nMax;

    //! Accumulate subspace counts
    virtual void addCount(int m);

    //! Scale factor for nearest-neighbor determined kernel width
    double bbfac;
    
    //! Target point list size for density construction
    int targetSize;
    
  public:

    /**  Minimum number of nearest neighbors for computing the kernel
	 density estimate (default: 32) */
    static int nnear;

    //! Target tolerance on density estimation
    static double errortol;

    //! Square of the distance threshold on individual state separation
    static double mindist;

    //! Minimum width for the smoothing kernel
    static double minwidth;

    /** Scale each dimension in parameter space seperately (default:
	true), otherwise use the parameter values in their natural
	metric (don't do this unless you are abosolutely sure it makes
	sense */
    static bool scaled;

    /** Scale each dimension in the KDE distance metric computation by
	the variance (default: true).  Don't change this unless you
	are abosolutely sure it makes sense */
    static bool mscale;

    //! Default constructor
    EnsembleDisc();

    //+ CLICONSTR StateInfo*
    /** Construct state from metadata info
	\param si is a pointer to the metadata structure
    */
    EnsembleDisc(StateInfo *si);

    //+ CLICONSTR StateInfo* int int string int
    /** General constructor
	\param si is a pointer to the metadata structure
	\param level is the update level
	\param nburn is the desired number of converged states
	\param filename is the logfile for ensemble statistics
	\param keypos is the parameter index for computing PDF and CDF 
    */
    EnsembleDisc(StateInfo *si,
		 int level, int nburn, string filename, int keypos);

    //! Copy constructor
    EnsembleDisc(const EnsembleDisc&);

    //! Copy constructor for CLI from base class
    //+ CLICONSTR Ensemble*
    EnsembleDisc(Ensemble*);

    //! Destructor
    virtual ~EnsembleDisc() {}
    
    //+ CLIMETHOD void Reset StateInfo*
    //! Reinitialize the EnsembleDisc from new state metadata
    void Reset(StateInfo *);

    //+ CLIMETHOD void Reset StateInfo* int int string int
    //! Reinitialize the EnsembleDisc
    void Reset(StateInfo* si,
	       int level, int nburn, string filename, int keypos);

    //+ CLIMETHOD void setDimensions StateInfo*
    //! Set total size of arrays from the state metadata structure
    void setDimensions(StateInfo *si);

    //+ CLIMETHOD void setTarget int
    /** Set the target number of points for density construction
	(default: 0 which means use the full list) */
    void setTarget(int n) { targetSize = n; }

    //+ CLIMETHOD void setBucket int
    //! Set the target number of points per ball-tree bucket
    void setBucket(int n) { bucketSize = n; }

    //+ CLIMETHOD void setNnear int
    //! Set the target number of nearest neighbors for kernel smoothing
    void setNnear(int n) {
      if (nnear == n) return;
      nnear = n; 
      dist_computed = kde_computed = kdeM_computed = false;
    }

    //+ CLIMETHOD void setRelative bool
    //! Use scaled relative distances in deduplication
    void setRelative (bool b) { scaled = b; }

    //+ CLIMETHOD void setDensityScale clivectord*
    //! User defined scaling for metric distance
    /** Sets a constant vector \f$s^2_j\f$ for the inverse metric coeffients.
	That is, the distance between points with indices \f$\mu\f$ and 
	\f$\nu\f$ given by \f$x^{[\mu]})\f$ and \f$x^{[\nu]}\f$, is:
	\f[
	ds^2 = \sum_j (x^{[\mu]}_j - x^{[\nu]}_j)^2/s^2_j
	\f]1
	This may provide a better representation of the kernel shape than 
	the variance for multimodal posterior distributions.
    */
    void setDensityScale(clivectord* f); 

    //+ CLIMETHOD void setKernelScale double
    //! Set the kernel width scale factor for nearest neighbor width estimate
    void setKernelScale (double s) {
      if (fabs(bbfac - s)<1.0e-16) return;
      bbfac = s; 
      dist_computed = kde_computed = kdeM_computed = false;
    }

    //+ CLIMETHOD void ComputeDistribution
    //+ CLIMETHOD void ComputeDistribution int
    //+ CLIMETHOD void ComputeDistribution int int
    /** 
	@name ComputeDistribution
	Computes the covariance matrix and its eigenvalues/vectors
    */
    void ComputeDistribution();
    void ComputeDistribution(int n) { burnIn=n; ComputeDistribution(); }
    void ComputeDistribution(int n, int nmax) { burnIn=n; nMax=nmax; ComputeDistribution(); }
    void ComputeDistributionMarginal(ostream& out);
    
    //+ CLIMETHOD void preComputeDensity
    //+ CLIMETHOD void preComputeDensity int
    //+ CLIMETHOD void preComputeDensity int int
    //@{
    /** 
	@name preComputeDensity
	Computes the KDE structures for later use (and archiving)
    */
    void preComputeDensity() {
      ComputeDistribution();
      ComputeDensity ();
      ComputeDensityM();
    }
    void preComputeDensity(int n) {
      ComputeDistribution(n);
      ComputeDensity ();
      ComputeDensityM();
    }
    void preComputeDensity(int n, int nmax) {
      ComputeDistribution(n, nmax);
      ComputeDensity ();
      ComputeDensityM();
    }
    //@}

    //! Object factory (clone)
    EnsembleDisc* New() {
      return new EnsembleDisc(_si);
    }

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
    double logPDFMarginal(unsigned m, unsigned n, const vector<double>& V);

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
    //+ CLIMETHOD void PrintDiag
    //+ CLIMETHOD void PrintDiag string
    //@{
    //! Use a supplied ostream for diag info
    void PrintDiag(ostream& out);

    //! Open a new file for diag into
    void PrintDiag(string& outfile);

    //! Print to console
    void PrintDiag();
    //@}

    //+ CLIMETHOD void PrintDensity int int int int string
    //! Density image (for debugging) using state vector with peak
    //! probability as the default
    void PrintDensity(int dim1, int dim2, int num1, int num2,
		      string file);

    //+ CLIMETHOD void PrintDensity int int int int clivectord* string
    //! Density image (for debugging) with user supplied default state
    //! vector
    void PrintDensity(int dim1, int dim2, int num1, int num2,
		      clivectord* def, string file);

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
      SAVEDMDV(widths,  vector<vector<double> >         );
      SAVEDMDV(widthsM, vector<vector<vector<double> > >);
      SAVEDMDV(weight,  vector<double>                  );
      SAVEDMDV(weightM, vector<vector<double> >         );
    }

    template<class Archive>
    void post_load(Archive &ar, const unsigned int file_version) 
    {
      RESTRMDV(points,  vector<vector<double> >         );
      RESTRMDV(pointsM, vector<vector<vector<double> > >);
      RESTRMDV(widths,  vector<vector<double> >         );
      RESTRMDV(widthsM, vector<vector<vector<double> > >);
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
          ar << BOOST_SERIALIZATION_NVP(kde_computed);                        
          BIE_CATCH_BOOST_SERIALIZATION_EXCEPTION;                     
         }                                                             
         try {                                                         
          ar << BOOST_SERIALIZATION_NVP(kdeM_computed);                        
          BIE_CATCH_BOOST_SERIALIZATION_EXCEPTION;                     
         }                                                             
         try {                                                         
          ar << BOOST_SERIALIZATION_NVP(bucketSize);                        
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
          ar << BOOST_SERIALIZATION_NVP(vlower);                        
          BIE_CATCH_BOOST_SERIALIZATION_EXCEPTION;                     
         }                                                             
         try {                                                         
          ar << BOOST_SERIALIZATION_NVP(vupper);                        
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
          ar << BOOST_SERIALIZATION_NVP(fscale);                        
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
          ar << BOOST_SERIALIZATION_NVP(kernel);                        
          BIE_CATCH_BOOST_SERIALIZATION_EXCEPTION;                     
         }                                                             
         try {                                                         
          ar << BOOST_SERIALIZATION_NVP(nMax);                        
          BIE_CATCH_BOOST_SERIALIZATION_EXCEPTION;                     
         }                                                             
         try {                                                         
          ar << BOOST_SERIALIZATION_NVP(bbfac);                        
          BIE_CATCH_BOOST_SERIALIZATION_EXCEPTION;                     
         }                                                             
         try {                                                         
          ar << BOOST_SERIALIZATION_NVP(targetSize);                        
          BIE_CATCH_BOOST_SERIALIZATION_EXCEPTION;                     
         }                                                             
         try {                                                         
          ar << BOOST_SERIALIZATION_NVP(nnear);                        
          BIE_CATCH_BOOST_SERIALIZATION_EXCEPTION;                     
         }                                                             
         try {                                                         
          ar << BOOST_SERIALIZATION_NVP(errortol);                        
          BIE_CATCH_BOOST_SERIALIZATION_EXCEPTION;                     
         }                                                             
         try {                                                         
          ar << BOOST_SERIALIZATION_NVP(mindist);                        
          BIE_CATCH_BOOST_SERIALIZATION_EXCEPTION;                     
         }                                                             
         try {                                                         
          ar << BOOST_SERIALIZATION_NVP(minwidth);                        
          BIE_CATCH_BOOST_SERIALIZATION_EXCEPTION;                     
         }                                                             
         try {                                                         
          ar << BOOST_SERIALIZATION_NVP(scaled);                        
          BIE_CATCH_BOOST_SERIALIZATION_EXCEPTION;                     
         }                                                             
         try {                                                         
          ar << BOOST_SERIALIZATION_NVP(mscale);                        
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
          ar >> BOOST_SERIALIZATION_NVP(kde_computed);                        
          BIE_CATCH_BOOST_SERIALIZATION_EXCEPTION;                     
         }                                                             
         try {                                                         
          ar >> BOOST_SERIALIZATION_NVP(kdeM_computed);                        
          BIE_CATCH_BOOST_SERIALIZATION_EXCEPTION;                     
         }                                                             
         try {                                                         
          ar >> BOOST_SERIALIZATION_NVP(bucketSize);                        
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
          ar >> BOOST_SERIALIZATION_NVP(vlower);                        
          BIE_CATCH_BOOST_SERIALIZATION_EXCEPTION;                     
         }                                                             
         try {                                                         
          ar >> BOOST_SERIALIZATION_NVP(vupper);                        
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
          ar >> BOOST_SERIALIZATION_NVP(fscale);                        
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
          ar >> BOOST_SERIALIZATION_NVP(kernel);                        
          BIE_CATCH_BOOST_SERIALIZATION_EXCEPTION;                     
         }                                                             
         try {                                                         
          ar >> BOOST_SERIALIZATION_NVP(nMax);                        
          BIE_CATCH_BOOST_SERIALIZATION_EXCEPTION;                     
         }                                                             
         try {                                                         
          ar >> BOOST_SERIALIZATION_NVP(bbfac);                        
          BIE_CATCH_BOOST_SERIALIZATION_EXCEPTION;                     
         }                                                             
         try {                                                         
          ar >> BOOST_SERIALIZATION_NVP(targetSize);                        
          BIE_CATCH_BOOST_SERIALIZATION_EXCEPTION;                     
         }                                                             
         try {                                                         
          ar >> BOOST_SERIALIZATION_NVP(nnear);                        
          BIE_CATCH_BOOST_SERIALIZATION_EXCEPTION;                     
         }                                                             
         try {                                                         
          ar >> BOOST_SERIALIZATION_NVP(errortol);                        
          BIE_CATCH_BOOST_SERIALIZATION_EXCEPTION;                     
         }                                                             
         try {                                                         
          ar >> BOOST_SERIALIZATION_NVP(mindist);                        
          BIE_CATCH_BOOST_SERIALIZATION_EXCEPTION;                     
         }                                                             
         try {                                                         
          ar >> BOOST_SERIALIZATION_NVP(minwidth);                        
          BIE_CATCH_BOOST_SERIALIZATION_EXCEPTION;                     
         }                                                             
         try {                                                         
          ar >> BOOST_SERIALIZATION_NVP(scaled);                        
          BIE_CATCH_BOOST_SERIALIZATION_EXCEPTION;                     
         }                                                             
         try {                                                         
          ar >> BOOST_SERIALIZATION_NVP(mscale);                        
          BIE_CATCH_BOOST_SERIALIZATION_EXCEPTION;                     
         }                                                             
        this->post_load(ar, version);
    }
    #endif

  };
  #ifndef SWIG
  BOOST_SERIALIZATION_SHARED_PTR(EnsembleDisc)
  #endif 
}
#ifndef SWIG
BIE_CLASS_TYPE_INFO(BIE::EnsembleDisc)
BIE_CLASS_EXPORT_KEY(BIE::EnsembleDisc)
#endif
#endif
