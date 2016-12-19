// This is really -*- C++ -*-


#ifndef GelmanRubinConverge_h
#define GelmanRubinConverge_h

#include <mpi.h>

#include <deque>

#include <Ensemble.h>
#include <Converge.h>

extern double inv_student_t_1sided(double alpha, double df);

#include "Serializable.h"


namespace BIE {
  
  //+ CLICLASS GelmanRubinConverge SUPER Converge
  //! Convergence class based on multiple chain analyses as desribed
  //! in Gelman and Rubin (1992)
  //!
  //! \par Intro and brief description
  //! This implements the convergence metric developed by Gelman and
  //! Rubin for measuring convergence of multiple chains.  The measure
  //! uses both the within- and between-chain variance estimate the
  //! total variance of the pooled chain variance.  From this they
  //! compute the ratio of the pooled variance to in-chain variance as
  //! an estimate of the reduction in variance that might be expected
  //! from a convergenced simulation.  They call this the "scale
  //! reduction".  A score of 1 indicates convergence.  Since a score
  //! of 1 is difficult to achieve, Gelman and Rubin recommend using
  //! 1.2 or 1.1 to declare convergence.  We assume 1.2 by default but
  //! this value may be changed using the #setRhatMax method.
  //!
  //! \par Outlier detection
  //! We use Grubbs' test (Grubbs 1969 and Stefansky 1972) to identify
  //! outliers in a univariate data set. It assumes that the
  //! underlying distribution is normal and this is appropriate for
  //! converged MCMC output.  Grubbs' test detects one outlier at a
  //! time, and successive outliers are iteratively removed from the
  //! dataset.
  //!
  //! \par Grubbs' test
  //! The Grubbs' test defined the following: 
  //! \f[ 
  //! G = {\max|x_i - {\bar x}|\over \sigma_s} 
  //! \f] 
  //! where \f$x_i\f$ are the \f$i=1,\ldots,N\f$ chain samples for a
  //! particular parameter, \f${\bar x}\f$ is the sample mean and
  //! \f$\sigma_s\f$ is the root of sample variance (standard
  //! deviation), respectively.  In other words, the Grubbs test
  //! statistic is the largest absolute deviation from the sample mean
  //! in units of the sample standard deviation.
  //! Let the confidence level be \f$\alpha\f$.  Then the hypothesis of no
  //! outliers is rejected if 
  //! \f[ G > {N-1\over\sqrt{N}}
  //! \sqrt{t^2_{\alpha/(2N),N-2}\over N - 2 + t^2_{\alpha/(2N),N-2}}
  //! \f]
  //! where \f$t_{\alpha/(2N),N-2}\f$ denoting the critical value of the
  //! Student t-distribution with \f$N-2\f$ degrees of freedom and a
  //! significance level of \f$\alpha/(2N)\f$.
  //!
  //! \par How to enable Grubbs' test
  //! To diagnose aberrant chains, Grubbs' test may be used by setting
  //! the outlier confidence \f$\alpha\f$ using the #setAlpha method.
  //! A value \f$\alpha=0\f$ implies no outlier testing.  Grubbs' test
  //! will only be applied after a specified number of steps have been
  //! computed; this is set by the #setNoutlier method.  The default
  //! value is 500.
  //!
  //! \par Other parameters
  //! The remaining parameters are similar to the SubsampleConverge
  //! routine.  In particular #setNgood chooses the number of states
  //! to compute after convergence is obtained and and #setNskip
  //! chooses the number of steps to skip between convergence checks.
  //!
  //! \ingroup converge
  //!
  class GelmanRubinConverge : public Converge
    {
    private:
      int next;
      int maxit;
      int nconv;
      string _id;
				// For notational convenience
      typedef vector<double> dvec;
				// Double-ended queue of vectors for parameters
      deque< vector<State> > data;
				// Double-ended queue for posterior probability
      deque< vector<double> > prob;
				// Currently valid chains
      vector<unsigned char> chain_mask;

				// Looks for abberant chains using Grubbs' test
      int ComputeMaskAll();
      int ComputeMaskRJ();
				// Two forms of the convergence test
      bool ConvergedAll();
      bool ConvergedRJ();

				// Keep count of instances for debugging
      static unsigned ninstance;

    public:
      /**@name Global variables */
      //@{
      //! Number of iterations after convergence (default: int=1000)
      static unsigned ngood;

      //! Number of iterations to skip betwen tests (default: int=100)
      static unsigned nskip;

      //! Outlier confidence value 1 - alpha (default: double=0.05)
      static double alpha;

      //! Maximum number of allowed outliers (default: int=6)
      static int maxoutlier;

      //! Number of steps before outlier test (default: int=500)
      static int noutlier;

      //! Tolerance for Rhat (default: double=1.2)
      static double rtol;

      //! Diagnostic files (default: bool=true)
      static bool verbose;

      //! Log of offset in peak probability for outliers (default: -30.0)
      static double poffset;
      
      //! Additional debug output (default: bool=false)
      static bool debug;

      //! Maximum number of states to keep if maxit=0 (default: 100000)
      static int maxkp;

      //@}

      //+ CLICONSTR int Ensemble* string
      /** 

      The parameter <code>m</code> determines the "burn-in" handling
      procedure.  If \f$m>0\f$, the most recently computed \f$m\f$
      steps are used to determine convergence.  If \f$m=0\f$, the most
      recent half of the steps are used. In all cases, \f$N_{skip}\f$
      retained steps are required before a convergence test is
      attempted.  Therefore if \f$m>0\f$, the first test will be
      computed after \f$max(m,N_{skip})\f$ steps.  <code>d</code> is
      an Ensemble instance used to accumulate statistics about the
      posterior simulation. <code>id</code> is used to tag diagnostic
      output files

      Note: this is a change from older versions that used \f$m<0\f$
      to request that the first \f$|m|\f$ steps be discarded.  This
      option resulted in nearly the same behavior as \f$m>0\f$ for
      long simulations and and otherwise caused confusion.  \f$m>=0\f$
      is now enforced on construction for compatibility with older
      applications.
      */
      GelmanRubinConverge(int m, Ensemble* d, string id);
      
      //! Factory method
      virtual GelmanRubinConverge* New(int m, Ensemble* d, string id);

      //! True if converged
      virtual bool Converged();
      
      //! This is a parallel chain convergence test
      virtual bool IsParallel() { return true; }
      
      //+ CLIMETHOD void setMax int
      //! Set number of states to retainfor next analysis.
      //! Set to zero to retain last half.
      void setMax(int n) { maxit = abs(n); }

      //+ CLIMETHOD void setNskip int
      //! Set number of iterations to skip between convergence tests
      void setNskip(int n) { nskip = n; next = count + nskip; }

      //+ CLIMETHOD void setNoutlier int
      //! Set number of iterations before performing an outlier test
      void setNoutlier(int n) { noutlier = n; }

      //+ CLIMETHOD void setMaxout int
      //! Set number maximum number of outliers that can be masked
      void setMaxout(int n) { maxoutlier = n; }

      //+ CLIMETHOD void setNgood int
      //! Set number of states to sample after convergence
      void setNgood(int n) { ngood = n; nburn = -1; next = count + nskip; }

      //+ CLIMETHOD void Quiet
      //! Suppress the diagnostic output
      void Quiet() { verbose = false; }

      //+ CLIMETHOD void setPoffset double
      /**
	 Set the minimum log probability offset from the maximally
	 probably state for an outlier to be flagged in addition
	 Grubbs' test.  In other words, outliers must also have
	 \f$\log(P)<\max\{\log(P)\}+P_{offset}\f$.  This additional
	 requirement may be entirely disabled by setting
	 \f$P_{offset}=0\f$.
      */
      void setPoffset(double z) { poffset = z; }

      //+ CLIMETHOD int ConvergedIndex
      //! Return the number of pre-burn-in states to discard from Ensemble
      int ConvergedIndex() { 
	MPI_Bcast(&nconv, 1, MPI_INT, 0, MPI_COMM_WORLD);
	return nconv; 
      }

      //+ CLIMETHOD void setRhatMax double
      //! Set correlation coefficient threshold (default: 1.2)
      void setRhatMax(double r) { rtol = r; }

      //+ CLIMETHOD void setAlpha double
      //! Set outlier detection confidence: 1 - alpha (default: 0, means off)
      void setAlpha(double a) { alpha = a; }

      //+ CLIMETHOD void setMaxK int32
      //! Set the maximum number of states to retain if maxit=0 (default: 100000)
      static void setMaxK(int n) { maxkp = n; }

      //! Cumulates data (overridden in SampleDistribution)
      virtual bool AccumData(vector< vector<double> >& values, 
			     vector<State>& states);
      
      //! Return last state
      virtual bool GetLast(vector<double>& values, vector<State>& states)
      {
	if (prob.size()>0) {
	  values = prob.back();
	  states = data.back();
	  return true;
	} else {
	  return false;
	}
      }
      
      /** Access to ComputDistribution class.  Will be called by provided
	  helper class as needed */
      void ComputeDistribution() {
	_dist->ComputeDistribution();
      }
      
      //! Clone me
      GelmanRubinConverge* New() {
	GelmanRubinConverge *p = new GelmanRubinConverge(maxit, _dist, _id);
	return p;
      }
      
      #ifndef SWIG
    // AUTO GENERATED BY ../persistence/autopersist.py
    protected:
    GelmanRubinConverge() {}
    private:
    friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive & ar, const unsigned int version) {
        this->pre_serialize(ar, version);
         try {                                                         
          ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(Converge);            
          BIE_CATCH_BOOST_SERIALIZATION_EXCEPTION;                     
         }                                                             
         try {                                                         
          ar & BOOST_SERIALIZATION_NVP(next);                        
          BIE_CATCH_BOOST_SERIALIZATION_EXCEPTION;                     
         }                                                             
         try {                                                         
          ar & BOOST_SERIALIZATION_NVP(maxit);                        
          BIE_CATCH_BOOST_SERIALIZATION_EXCEPTION;                     
         }                                                             
         try {                                                         
          ar & BOOST_SERIALIZATION_NVP(nconv);                        
          BIE_CATCH_BOOST_SERIALIZATION_EXCEPTION;                     
         }                                                             
         try {                                                         
          ar & BOOST_SERIALIZATION_NVP(_id);                        
          BIE_CATCH_BOOST_SERIALIZATION_EXCEPTION;                     
         }                                                             
         try {                                                         
          ar & BOOST_SERIALIZATION_NVP(data);                        
          BIE_CATCH_BOOST_SERIALIZATION_EXCEPTION;                     
         }                                                             
         try {                                                         
          ar & BOOST_SERIALIZATION_NVP(prob);                        
          BIE_CATCH_BOOST_SERIALIZATION_EXCEPTION;                     
         }                                                             
         try {                                                         
          ar & BOOST_SERIALIZATION_NVP(chain_mask);                        
          BIE_CATCH_BOOST_SERIALIZATION_EXCEPTION;                     
         }                                                             
         try {                                                         
          ar & BOOST_SERIALIZATION_NVP(ninstance);                        
          BIE_CATCH_BOOST_SERIALIZATION_EXCEPTION;                     
         }                                                             
         try {                                                         
          ar & BOOST_SERIALIZATION_NVP(ngood);                        
          BIE_CATCH_BOOST_SERIALIZATION_EXCEPTION;                     
         }                                                             
         try {                                                         
          ar & BOOST_SERIALIZATION_NVP(nskip);                        
          BIE_CATCH_BOOST_SERIALIZATION_EXCEPTION;                     
         }                                                             
         try {                                                         
          ar & BOOST_SERIALIZATION_NVP(alpha);                        
          BIE_CATCH_BOOST_SERIALIZATION_EXCEPTION;                     
         }                                                             
         try {                                                         
          ar & BOOST_SERIALIZATION_NVP(maxoutlier);                        
          BIE_CATCH_BOOST_SERIALIZATION_EXCEPTION;                     
         }                                                             
         try {                                                         
          ar & BOOST_SERIALIZATION_NVP(noutlier);                        
          BIE_CATCH_BOOST_SERIALIZATION_EXCEPTION;                     
         }                                                             
         try {                                                         
          ar & BOOST_SERIALIZATION_NVP(rtol);                        
          BIE_CATCH_BOOST_SERIALIZATION_EXCEPTION;                     
         }                                                             
         try {                                                         
          ar & BOOST_SERIALIZATION_NVP(verbose);                        
          BIE_CATCH_BOOST_SERIALIZATION_EXCEPTION;                     
         }                                                             
         try {                                                         
          ar & BOOST_SERIALIZATION_NVP(poffset);                        
          BIE_CATCH_BOOST_SERIALIZATION_EXCEPTION;                     
         }                                                             
         try {                                                         
          ar & BOOST_SERIALIZATION_NVP(debug);                        
          BIE_CATCH_BOOST_SERIALIZATION_EXCEPTION;                     
         }                                                             
         try {                                                         
          ar & BOOST_SERIALIZATION_NVP(maxkp);                        
          BIE_CATCH_BOOST_SERIALIZATION_EXCEPTION;                     
         }                                                             
        this->post_serialize(ar, version);
    }
    #endif


    };

  /** @} */

}

#ifndef SWIG
BIE_CLASS_TYPE_INFO(BIE::GelmanRubinConverge)
BIE_CLASS_EXPORT_KEY(BIE::GelmanRubinConverge)
#endif
#endif
