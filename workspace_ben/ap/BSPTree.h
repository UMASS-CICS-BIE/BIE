// -*- C++ -*-

//=============================================================================
// BSP tree class
//=============================================================================

#ifndef _BSPTREE_h
#define _BSPTREE_h

#include <cmath>
#include <cassert>
#include <iostream>
#include <vector>
#include <map>
#include <set>

#include <boost/shared_ptr.hpp>

#include <RndInt.h>
#include <State.h>

class ACG;
class Normal;

#include "Serializable.h"


/**     
    \class LeafData
    \brief Store the data for the measure function

    \author Martin Weinberg
    \date July 2010

    An internal data structure for computing the measure function 
    for Lebesgue integration. */
class LeafData: public Serializable
{
 public:

  //! Scaling for Lebesgue volume (default: 0=logarithmic)
  static double exp_scale;


  //! Constructor
  LeafData() {}

  //! Probability values
  vector< pair<double, double> > pvals;

  //! Minimum probability value
  double minP;

  //! Median probability value
  double medP;

  //! Maximum probability value
  double maxP;

  //! Volume in this cell
  double vol;

  //! Compute fractional volume
  double volume(double p);

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
          ar & BOOST_SERIALIZATION_NVP(pvals);                        
          BIE_CATCH_BOOST_SERIALIZATION_EXCEPTION;                     
         }                                                             
         try {                                                         
          ar & BOOST_SERIALIZATION_NVP(minP);                        
          BIE_CATCH_BOOST_SERIALIZATION_EXCEPTION;                     
         }                                                             
         try {                                                         
          ar & BOOST_SERIALIZATION_NVP(medP);                        
          BIE_CATCH_BOOST_SERIALIZATION_EXCEPTION;                     
         }                                                             
         try {                                                         
          ar & BOOST_SERIALIZATION_NVP(maxP);                        
          BIE_CATCH_BOOST_SERIALIZATION_EXCEPTION;                     
         }                                                             
         try {                                                         
          ar & BOOST_SERIALIZATION_NVP(vol);                        
          BIE_CATCH_BOOST_SERIALIZATION_EXCEPTION;                     
         }                                                             
        this->post_serialize(ar, version);
    }
    #endif

};


/** 
    \class BSPTree
    \brief Abstract interface for BSP trees

    Interface for a BSP tree.  Currently implemented are KD-trees and
    ORB-trees (see KDTree and ORBTree, respectively).  Curently, only
    KD has density estimation capabilities.

    \see KDTree
    \see ORBTree
    \see TNTree
 */
class BSPTree: public Serializable
{

 public:

  /// A mnemonic for all indices
  typedef unsigned long indx;

  /// Lower error bound quantile (default: 0.10)
  static double lower_quant;

  /// Upper error bound quantile (default: 0.90)
  static double upper_quant;

  /** Use the volume determined by cell size rather than spanned points
      (this is the default) */
  static bool full_volume;

  /** Use the volume determined by the geometric mean of the cell size
      and the volume spanned by the points (default: false) */
  static bool geom_volume;

  /** Switch between full_volume and spanned volume depending on the
      filling factor */
  static bool fill_volume;

  /** Filling factor threshold for switching heuristic */
  static double fill_factor;

  //@{
  /** Default slice dimensions */

  //! First axis (Default: 0)
  static unsigned dim1;

  //! Second axis (Default: 1)
  static unsigned dim2;
  //@}

  /// Null constructor
  BSPTree();

  /// Destructor
  virtual ~BSPTree() {};

  //@{
  /// Functions for computing moments on the implied measure

  virtual void IntegralEval(std::vector< std::vector<double> >& data, 
			    std::vector<unsigned>& mult,
			    std::map<int, std::vector<double> >& result, 
			    std::map<int, std::vector<double> >& lower, 
			    std::map<int, std::vector<double> >& upper,
			    std::map<int, std::vector<double> >& vmean) = 0;

  virtual void LevelList(std::ostream& out, bool limits) = 0;

  virtual void IntegralList(std::ostream& out,
			    std::vector< std::vector<double> >& data,
			    std::vector<unsigned>& mult,
			    bool allcells=false) = 0;

  virtual void IntegralCellList(std::ostream& out, unsigned long ncut,
				std::vector< std::vector<double> >& data) = 0;

  //@}


  //@{
  /** Lebesgue integration

      Since this feature is in common between both KD and ORB, it is
      implemented as part of the BSPTree class.
  */

  /** Do the integral: returns the result of the lower Riemann, upper
      Riemann and the mean of the lower and upper Riemann sums.  Use
      LebesgueMeasure below to get an error analysis.
  */
  void LebesgueIntegral(std::map<int, double>& lower, 
			std::map<int, double>& mean, 
			std::map<int, double>& upper,
			unsigned knots, unsigned dindx, unsigned long ncut, 
			std::vector< std::vector<double> >& data,
			std::vector<unsigned>& mult,
			bool logM=false, double dlogM=1e20, double dlogP=1e20,
			bool linear=false);

  //@{
  /** For analysis: computes the measure function, contribution to the
      integral from each grid interval along with the absolute and
      relative errors

      @param P will return a vector of probability knots
      @param M will return a vector of measure function evals at those knots
      @param V will return a vector of measure integrand contributions
      @param R will return a vector of relative integrand errors
      @param A will return a vector of absolute integrand errors
      @param nval is the target numer of knots
      @param dindx is the index of the probability field in the data array
      @param ncut is the number of data points per volume element
      @param data is the input data array (each record has a number of fields
      		that may be selected using dindx)
      @param mult is the multiplicity for each unique point
      @param logM may be set to true for log scaling the prob value in the 
      		Lebesgue integral
      @param dlogM is the maximum offset from the peak probability value to 
      		the minimum probability value
      @param dlogP is the maximum interval between quadrature knots; nval 
      		will be increased to achieve this spacing
      @param linear uses linear rather than logarithmic mapping to assign 
      		cell volume to probability values
  */
  void LebesgueMeasure(std::vector<double>& P, std::vector<double>& M,
		       std::vector<double>& V, std::vector<double>& R, 
		       std::vector<double>& A,
		       unsigned nval, unsigned dindx, unsigned long ncut, 
		       std::vector< std::vector<double> >& data, 
		       std::vector<unsigned>& mult, int m, 
		       bool logM=false, double dlogM=1e20, double dlogP=1e20,
		       bool linear=false);

  void LebesgueMeasure(std::vector<double>& P, std::vector<double>& M,
		       std::vector<double>& V, std::vector<double>& R, 
		       std::vector<double>& A,
		       double Pmin, double Pmax,
		       unsigned nval, unsigned dindx, unsigned long ncut, 
		       std::vector< std::vector<double> >& data,
		       std::vector<unsigned>& mult, int m,
		       bool logM=false, double dlogM=1e20, double dlogP=1e20,
		       bool linear=false);
  //@}

  /// Get the volume of the cells (for debugging)
  double LebesgueVolume(int m);

  /// Get the median volume of the cells
  double MedianVolume(int m);

  /// Evaluate the measure function
  double LebesgueEval(int m, double y);

  //@}

  /// Set the volume trim
  virtual void SetVTrim(double fraction) 
  { std::cerr << "Volume trim not defined for this tree" << std::endl; }

  /// Get the volume fraction
  virtual double GetVTrim() { return vfrac; }

  /// Get prob--volume and count--volume historgrams for debugging
  void LebesgueHisto(std::vector<double>& P,
		     std::vector< std::pair<double, double> >& histoP,
		     std::map<int, double>& histoC,
		     unsigned dindx, unsigned long ncut,
		     std::vector< std::vector<double> >& data,
		     std::vector<unsigned>& mult, int m);

  /** Produce a gnuplot script illustrating the tessellation.
      
      Arguments are the same as in PDFList
  */
  virtual void gnuplotDemo(const std::vector< std::vector<double> >& data,
			   const std::vector< unsigned >& mult,
			   const std::vector< double >& center,
			   const unsigned dindx, const unsigned long ncut) = 0;
  
  /// Volume in the region enclosing the points
  virtual double EnclosedVolume();

  /// Set global range bounds
  static void SetRangeBounds(std::vector<double>& lo, std::vector<double>& up);

  //@{
  /// Get range bounds
  std::vector<double> GetLowerBound() { return lowerBound; }
  std::vector<double> GetUpperBound() { return upperBound; }
  //@}

protected:

  //@{
  /// Global range bounds
  static std::vector<double> lowerBound;
  static std::vector<double> upperBound;
  //@}

  /// Number in leaf buckets
  unsigned ncut;

  /// dimension of data
  unsigned long dims;

  /// scaling vector
  std::vector<double> scale;

  /// Find the bounding region
  virtual void makeBounds() = 0;

  /// dimension of result vector
  unsigned long sz;

  /// # of points 
  unsigned long num_points;

  //@{
  /// Lebesgue data and routines

  std::map<int, std::multimap<double, LeafData> > lebesgueLeaves;
  std::map<int, double> LminP, LmaxP;

  /** Create the volume structure for the measure function

      @param dindx is the index of the probability field in the data array
      @param ncut is the number of data points per volume element
      @param data is the input data array (each record has a number of 
      		fields that may be selected using dindx)
      @param mult is the multiplicity for each point in data
   */
  virtual void createLebesgueList(unsigned dindx, unsigned long ncut,
				  std::vector< std::vector<double> >& data,
				  std::vector<unsigned>& mult,
				  bool linear) = 0;

  //! See LebesgueMeasure
  void ProcessMeasure(std::vector<double>& Pval, 
		      std::vector<double>& Mval,
		      std::vector<double>& Yval, 
		      std::vector<double>& Erel, 
		      std::vector<double>& Eabs,
		      double Pmin, double Pmax,
		      unsigned nval, unsigned dindx, unsigned long ncut,
		      std::vector< std::vector<double> >& data,
		      std::vector<unsigned>& mult, int m,
		      bool logM, double dlogM, double dlogP);
  //@}

  typedef std::pair<double, double> dpair;
  
  //! Return random string file for debug file labeling
  static std::string getRandomFileTag() {
    std::ostringstream sout;
    sout << std::hex << std::setw(10) << std::setfill('0') << std::right 
	 << (*tags)();
    return sout.str();
  }

  /// Trimmed volume fraction
  double vfrac;

private:

  // This is a static shared pointer.  Since it will never go out of
  // scope, this may be overkill, but it should be clean at least.
  //
  static RndIntPtr tags;

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
          ar & BOOST_SERIALIZATION_NVP(lower_quant);                        
          BIE_CATCH_BOOST_SERIALIZATION_EXCEPTION;                     
         }                                                             
         try {                                                         
          ar & BOOST_SERIALIZATION_NVP(upper_quant);                        
          BIE_CATCH_BOOST_SERIALIZATION_EXCEPTION;                     
         }                                                             
         try {                                                         
          ar & BOOST_SERIALIZATION_NVP(full_volume);                        
          BIE_CATCH_BOOST_SERIALIZATION_EXCEPTION;                     
         }                                                             
         try {                                                         
          ar & BOOST_SERIALIZATION_NVP(geom_volume);                        
          BIE_CATCH_BOOST_SERIALIZATION_EXCEPTION;                     
         }                                                             
         try {                                                         
          ar & BOOST_SERIALIZATION_NVP(fill_volume);                        
          BIE_CATCH_BOOST_SERIALIZATION_EXCEPTION;                     
         }                                                             
         try {                                                         
          ar & BOOST_SERIALIZATION_NVP(fill_factor);                        
          BIE_CATCH_BOOST_SERIALIZATION_EXCEPTION;                     
         }                                                             
         try {                                                         
          ar & BOOST_SERIALIZATION_NVP(dim1);                        
          BIE_CATCH_BOOST_SERIALIZATION_EXCEPTION;                     
         }                                                             
         try {                                                         
          ar & BOOST_SERIALIZATION_NVP(dim2);                        
          BIE_CATCH_BOOST_SERIALIZATION_EXCEPTION;                     
         }                                                             
         try {                                                         
          ar & BOOST_SERIALIZATION_NVP(lowerBound);                        
          BIE_CATCH_BOOST_SERIALIZATION_EXCEPTION;                     
         }                                                             
         try {                                                         
          ar & BOOST_SERIALIZATION_NVP(upperBound);                        
          BIE_CATCH_BOOST_SERIALIZATION_EXCEPTION;                     
         }                                                             
         try {                                                         
          ar & BOOST_SERIALIZATION_NVP(ncut);                        
          BIE_CATCH_BOOST_SERIALIZATION_EXCEPTION;                     
         }                                                             
         try {                                                         
          ar & BOOST_SERIALIZATION_NVP(dims);                        
          BIE_CATCH_BOOST_SERIALIZATION_EXCEPTION;                     
         }                                                             
         try {                                                         
          ar & BOOST_SERIALIZATION_NVP(scale);                        
          BIE_CATCH_BOOST_SERIALIZATION_EXCEPTION;                     
         }                                                             
         try {                                                         
          ar & BOOST_SERIALIZATION_NVP(sz);                        
          BIE_CATCH_BOOST_SERIALIZATION_EXCEPTION;                     
         }                                                             
         try {                                                         
          ar & BOOST_SERIALIZATION_NVP(num_points);                        
          BIE_CATCH_BOOST_SERIALIZATION_EXCEPTION;                     
         }                                                             
         try {                                                         
          ar & BOOST_SERIALIZATION_NVP(lebesgueLeaves);                        
          BIE_CATCH_BOOST_SERIALIZATION_EXCEPTION;                     
         }                                                             
         try {                                                         
          ar & BOOST_SERIALIZATION_NVP(LminP);                        
          BIE_CATCH_BOOST_SERIALIZATION_EXCEPTION;                     
         }                                                             
         try {                                                         
          ar & BOOST_SERIALIZATION_NVP(LmaxP);                        
          BIE_CATCH_BOOST_SERIALIZATION_EXCEPTION;                     
         }                                                             
        this->post_serialize(ar, version);
    }
    #endif

};

#ifndef SWIG
BIE_CLASS_ABSTRACT(BSPTree)
BOOST_SERIALIZATION_SHARED_PTR(LeafData)
BOOST_SERIALIZATION_SHARED_PTR(BSPTree)

BIE_CLASS_TYPE_INFO(LeafData)
BIE_CLASS_EXPORT_KEY(LeafData)
BIE_CLASS_EXPORT_KEY(BSPTree)
#endif
#endif
