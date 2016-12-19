// -*- C++ -*-

//=============================================================================
// KD-tree base class to be used for density estimation
//=============================================================================

#ifndef _KDTREE_h
#define _KDTREE_h

#include <BSPTree.h>

class ACG;
class Uniform;
class Normal;

#include "Serializable.h"


/** An internal data structure for computing the measure function and
    density by KDTree.  This structure is a simple class to simplify
    serialization. */
class LeafElem: public Serializable
{
 public:
  
  //! Constructor
  LeafElem() {}

  /** This key is the fraction of the unit interval foliated into the
      KD structure */
  double key;

  //! This is the CDF value after constuction
  double cdf;
  
  //! The likelihood contribution of this cell
  double dL;
  
  //! The median likelihood value
  double L;

  //! The posterior value
  double P;

  //! Weight in cell
  double W;

  //! Volume in cell
  double V;

  //! Geometric volume in the cell
  double G;

  //! The cell index
  unsigned long cell;
  
  //! The index of the point with median likelihood value
  unsigned long point;

  //! The number of points in this cell
  unsigned long N;
  
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
          ar & BOOST_SERIALIZATION_NVP(key);                        
          BIE_CATCH_BOOST_SERIALIZATION_EXCEPTION;                     
         }                                                             
         try {                                                         
          ar & BOOST_SERIALIZATION_NVP(cdf);                        
          BIE_CATCH_BOOST_SERIALIZATION_EXCEPTION;                     
         }                                                             
         try {                                                         
          ar & BOOST_SERIALIZATION_NVP(dL);                        
          BIE_CATCH_BOOST_SERIALIZATION_EXCEPTION;                     
         }                                                             
         try {                                                         
          ar & BOOST_SERIALIZATION_NVP(L);                        
          BIE_CATCH_BOOST_SERIALIZATION_EXCEPTION;                     
         }                                                             
         try {                                                         
          ar & BOOST_SERIALIZATION_NVP(P);                        
          BIE_CATCH_BOOST_SERIALIZATION_EXCEPTION;                     
         }                                                             
         try {                                                         
          ar & BOOST_SERIALIZATION_NVP(W);                        
          BIE_CATCH_BOOST_SERIALIZATION_EXCEPTION;                     
         }                                                             
         try {                                                         
          ar & BOOST_SERIALIZATION_NVP(V);                        
          BIE_CATCH_BOOST_SERIALIZATION_EXCEPTION;                     
         }                                                             
         try {                                                         
          ar & BOOST_SERIALIZATION_NVP(G);                        
          BIE_CATCH_BOOST_SERIALIZATION_EXCEPTION;                     
         }                                                             
         try {                                                         
          ar & BOOST_SERIALIZATION_NVP(cell);                        
          BIE_CATCH_BOOST_SERIALIZATION_EXCEPTION;                     
         }                                                             
         try {                                                         
          ar & BOOST_SERIALIZATION_NVP(point);                        
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

/** 
    \class CellElem
    \brief A cell info type for producing diagnostics
    
    \author Martin Weinberg
    \date May 2011
*/
class CellElem: public Serializable
{
public:
  
  //! Constructor
  CellElem() {}
  
  //@{
  //! The parameter space intervals
  std::vector<double> lo, hi;
  //@}
  
  //! The data value
  double P;
  
  //! The index of the point with median likelihood value
  unsigned long point;
  
  //! The number of points in this cell
  unsigned long N;
  
  //! The total volume
  double V;

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
          ar & BOOST_SERIALIZATION_NVP(lo);                        
          BIE_CATCH_BOOST_SERIALIZATION_EXCEPTION;                     
         }                                                             
         try {                                                         
          ar & BOOST_SERIALIZATION_NVP(hi);                        
          BIE_CATCH_BOOST_SERIALIZATION_EXCEPTION;                     
         }                                                             
         try {                                                         
          ar & BOOST_SERIALIZATION_NVP(P);                        
          BIE_CATCH_BOOST_SERIALIZATION_EXCEPTION;                     
         }                                                             
         try {                                                         
          ar & BOOST_SERIALIZATION_NVP(point);                        
          BIE_CATCH_BOOST_SERIALIZATION_EXCEPTION;                     
         }                                                             
         try {                                                         
          ar & BOOST_SERIALIZATION_NVP(N);                        
          BIE_CATCH_BOOST_SERIALIZATION_EXCEPTION;                     
         }                                                             
         try {                                                         
          ar & BOOST_SERIALIZATION_NVP(V);                        
          BIE_CATCH_BOOST_SERIALIZATION_EXCEPTION;                     
         }                                                             
        this->post_serialize(ar, version);
    }
    #endif

};


/** Computes a KDTree structure from a vector of a vector of points
 */
class KDTree : public BSPTree
{
  friend class Kernel;
  friend class Delta;

 public:

  /// Verbose debugging flag (default: false)
  static bool verbose_debug;

  /// Indicates no further children
  constexpr static indx NO_CHILD = (unsigned long) -1;  
  
  /// Use range or disperion
  static bool use_range;


  /// Null constructor
  KDTree() {}

  /// Constructor
  KDTree( unsigned long d, unsigned long N, unsigned ncut,
	  const double* points_, const double* weights_,
	  const double* beta_);
  
  /// Constructor
  KDTree( unsigned long d, unsigned long N, unsigned ncut,
	  const std::vector< std::vector<double> > &points_, 
	  const std::vector<double> &weights_, 
	  const std::vector<double> &beta_);
  
  /// Destructor
  virtual ~KDTree() {};

  /// Factory
  KDTree *New();

  //@{
  /// Accessor Functions  
  KDTree::indx    root() const           { return 0;    }
  unsigned int    Ndim() const           { return dims; }
  unsigned long   Npts()                 const { return num_points; }
  unsigned long   Npts(KDTree::indx i)   const { return hi_leaf[i]-lo_leaf[i]+1; }
  
  const double* center(KDTree::indx i)   const {if (i<num_points) return &center0[i*dims]; else return &centers[i*dims]; }
  
  const double* range(KDTree::indx i)    const { return &ranges[i*dims]; }
  const double* Range(KDTree::indx i)    const { return &range0[i*dims]; }
  const double* Low  (KDTree::indx i)    const { return &Ledge[i*dims];  }
  const double* High (KDTree::indx i)    const { return &Hedge[i*dims];  }
  
  double  Vol       (KDTree::indx i);
  double  VOL       (KDTree::indx i);
  double  weight    (KDTree::indx i)     const { return *(&weights[0]+i); }
  bool    isLeaf    (KDTree::indx ind)   const { return ind >= num_points; }
  bool    validIndex(KDTree::indx ind)   const { return ((0<=ind) && (ind < 2*num_points)); }
  
  KDTree::indx left(KDTree::indx i)      const { return left_child[i]; }
  KDTree::indx right(KDTree::indx i)     const { return right_child[i]; }
  KDTree::indx leafFirst(KDTree::indx i) const { return lo_leaf[i]; }
  KDTree::indx leafLast(KDTree::indx i)  const { return hi_leaf[i]; }
  
  /// Convert a KDTree::indx to the numeric indx in the original data
  unsigned long getIndexOf(KDTree::indx i) const { return permutation[i]; }

  void movePoints(double*);
  void changeWeights(const double *);
  //@}

  //@{
  /// Test two sub-trees to see which is nearer another KDTree
  KDTree::indx closer(KDTree::indx, KDTree::indx, const KDTree&,KDTree::indx) const;  

  KDTree::indx closer(KDTree::indx i, KDTree::indx j, const KDTree& other_tree) const
  { return closer(i,j,other_tree,other_tree.root()); };
  //@}

  //@{
  /// Nearest neighbor functions
  KDTree::indx getClosestPoint(std::vector<double>& T, KDTree::indx root,
			       std::vector< std::vector<double> >& data);

  void kNearestNeighbors(indx *, double *, const double *, int, int) const;
  //@}

  //@{
  /// Diagnostic output
  void Diagnostics(std::ostream& out, unsigned long ncut);

  void Diagnostics(std::ostream& out, unsigned long ncut, int dim,
		   std::vector<double> (*func)(std::vector<double>&));
  //@}

  //@{
  /// Functions for computing moments on the implied measure

  void Expectation(unsigned long ncut,
		   std::vector< std::vector<double> >& data, 
		   std::vector<double>& result);

  void MeasureList(std::ostream& out, unsigned long ncut,
		   std::vector< std::vector<double> >& data);

  void createMeasureList(unsigned long ncut,
			 std::vector< std::vector<double> >& data);
  
  KDTree::indx find(const std::vector<double>& T, unsigned long ncut);
  KDTree::indx find(const std::vector<double>& T) { return find(T, ncutM); }
  
  double FTheta(const std::vector<double>& T);
  KDTree::indx FThetaInv(double u);
  double PosteriorProb(const std::vector<double>& T);
  
  void MedianEval(unsigned long ncut, 
		  std::vector< std::vector<double> >& data, 
		  std::vector<double>& result);

  void LevelList(std::ostream& out, bool limits);

  //@}

  /// Get the density in the cell containing @param T at teperature
  /// @param Temp (Temp=1 by default)
  double Density(const vector<double>& T, const double Temp=1.0);

  /// Return a sample from the density density field at teperature
  /// @param Temp (Temp=1 by default)
  std::vector<double> DensitySample(const double Temp=1.0);

  //@{
  /// For debugging

  /// Check to see that the particles are between their boundaries
  unsigned checkParticles(KDTree::indx low, KDTree::indx high,
			  const vector<double>& minC, 
			  const vector<double>& maxC);
  //@}
  
  
  /** Produce a gnuplot script illustrating the tessellation.
      
      Arguments are the same as in PDFList
  */
  void gnuplotDemo(const std::vector< std::vector<double> >& data,
		   const std::vector< unsigned >& mult,
		   const std::vector< double >& center,
		   const unsigned dindx, const unsigned long ncut);
  
  //! Trim the volume to include a fraction of the original points
  virtual void SetVTrim(double fraction)
  {
    vol_frac     = fraction;
    vol_trim     = true;
    vol_computed = false;
  }

protected:

  /// For jitter in range selection
  //@{
  /// Generator
  boost::shared_ptr<ACG>     gen;
  /// Variates
  boost::shared_ptr<Uniform> unit;
  boost::shared_ptr<Normal>  norm;
  //@}

  /// Split tally (for debugging)
  std::vector<unsigned> stally;

  /// construction recursion
  virtual void calcStats(KDTree::indx);


  //@{
  /// KD centers
  std::vector<double> centers, center0;
  //@}

  //@{
  /// bounding ranges for entire volume
  std::vector<double> minS, maxS;
  //@}

  //@{
  /// bounding box ranges, dims per KD, dist from center to one side
  std::vector<double> ranges, range0;
  //@}

  /// lower bounds on cells
  std::vector<double> Ledge;

  /// upper bounds on cells
  std::vector<double> Hedge;

  /// total weight in each KD 
  std::vector<double> weights;

  /// Beta value
  std::vector<double> beta;

  /// lower measure interval
  std::vector<double> lower;
  
  /// upper measure interval
  std::vector<double> upper;
  
  //@{
  /// left, right children; no parent indices
  std::vector<KDTree::indx> left_child,  right_child;
  //@}

  //@{
  /// lower & upper leaf indices for each KD
  std::vector<KDTree::indx> lo_leaf, hi_leaf;
  //@}

  /// point's position in the original data
  std::vector<KDTree::indx> permutation;

  /// internal var for placing the non-leaf nodes 
  KDTree::indx next;

  //@{
  /// Minimum and maximum values of log likelihood
  double minL, maxL;
  //@}

  //@{
  /// for computing the inverse distribution
  std::map<double, LeafElem> measureList;
  unsigned long ncutM;
  //@}

  /// for mapping the measure to Omega
  std::map<indx, LeafElem> measureMap;

  /// Density normalzation
  double dnorm;

  /// Density heating factor (power)
  double temp;

  /// Compute density norm
  void computeDensity(double T);

  /// for mapping cdf to cell index for density sampling
  std::map<double, unsigned long> densityMap;

  /// Build the KD tree
  void buildKD(KDTree::indx firstLeaf, KDTree::indx lastLeaf, 
	       KDTree::indx root,
	       const std::vector<double>& minc,
	       const std::vector<double>& maxc);

  /// Compute the data range
  void makeRange(KDTree::indx, KDTree::indx, KDTree::indx root);

  /** Find the dimension along which the leaves between low and high
      inclusive have the greatest variance */
  KDTree::indx max_dispersion(KDTree::indx, KDTree::indx) const;

  /** Find the dimension along which the leaves between low and high
      inclusive have the greatest range */
  KDTree::indx max_range(KDTree::indx) const;


  /// The protected method for build the tree beginning at index root
  KDTree::indx find(const std::vector<double>& T, unsigned long ncut, 
		    KDTree::indx root);

  /// leaf-swapping function
  virtual void swap(KDTree::indx, KDTree::indx);         

  /** Function to partition the data into two (equal-sized or near as
      possible) sets, one of which is uniformly greater than the other in
      the given dimension.  This is a partial binary sort. */
  void select(unsigned long dimension, unsigned long position, 
	      unsigned long low, unsigned long high);

  /// Returns distance squared (minimum)
  double minDist(indx, const double*) const;

  /// Returns distance squared (maximum)
  double maxDist(indx, const double*) const;

  /// build the non-leaf nodes from the leaves
  void buildTree();

  /// For computation against the tree
  void Volume(KDTree::indx root, unsigned long ncut, double& volume);

  /// For computation against the tree
  void Expect(KDTree::indx root, unsigned long ncut, 
	      std::vector< std::vector<double> >& data);

  /// For computation against the tree
  void Median(KDTree::indx root, unsigned long ncut, 
	      std::vector< std::vector<double> >& data);

  /// For computation against the tree
  void NodeList(KDTree::indx root, unsigned lev, unsigned long ncut,
		std::set< std::pair<indx, unsigned> >& leaves);
  
  /// Helper routine for gnuplot tessellation plot
  void cellList(KDTree::indx root, 
		const unsigned long ncut,
		const unsigned dindx, 
		const std::vector< std::vector<double> >& data,
		const std::vector<unsigned>& mult,
		std::vector<CellElem>& leaves);
  //@{
  /// Integration routines

  void PDFList(KDTree::indx root, unsigned long ncut,
	       std::vector< std::vector<double> >& tree,
	       std::multimap<double, LeafElem>& leaves);

  void Integral(KDTree::indx root, unsigned long ncut,
		std::vector<double>& result,
		std::vector<double> (*func)(std::vector<double>&));

  void IntegralEval(std::vector< std::vector<double> >& data, 
		    std::vector<unsigned>& mult,
		    std::map<int, std::vector<double> >& result,
		    std::map<int, std::vector<double> >& lower, 
		    std::map<int, std::vector<double> >& upper,
		    std::map<int, std::vector<double> >& vmean);
  
  void IntegralPC(KDTree::indx root, unsigned long ncut,
		  std::vector< std::vector<double> >& tree,
		  std::vector<unsigned>& mult,
		  std::map<int, std::vector< std::vector<double> > >& medn,
		  std::map<int, std::vector< std::vector<double> > >& lowr,
		  std::map<int, std::vector< std::vector<double> > >& uppr,
		  std::map<int, std::vector< std::vector<double> > >& mean);
  
  void IntegralList(std::ostream& out,
		    std::vector< std::vector<double> >& data,
		    std::vector<unsigned>& mult,
		    bool allcells);

  void IntegralCellList(std::ostream& out, unsigned long ncut,
			std::vector< std::vector<double> >& data);


  //@}		  
  
  /// Find the bounding region
  void makeBounds();
  
  
  //@{
  /// For tree trimming

  //! Return true if the cell is in the trimmed frontier
  bool TrimOK(indx root) {
    bool ok = true;
    if (vol_trim) {
      if (trim_list.find(root) == trim_list.end()) ok = false;
    }
    return ok;
  }

  //@{
  //! Compute the trimmed frontier
  double computeTrim(std::vector< std::vector<double> >& tree);
  void findTrimList(KDTree::indx root,
		    const unsigned long ncut,
		    std::vector< std::vector<double> >& tree,
		    std::multimap<double, KDTree::indx> &vol_frontier,
		    unsigned long &vol_count);
  std::set<KDTree::indx> trim_list;
  double vol_frac;
  bool vol_trim, vol_computed;
  //@}
  //@}

  //@{
  /** Create the volume structure for the measure function
      
      @param dindx is the index of the probability field in the data array
      @param ncut is the number of data points per volume element
      @param data is the input data array (each record has a number of 
      fields that may be selected using dindx)
  */
  void createLebesgueList(unsigned dindx, unsigned long ncut,
			  std::vector< std::vector<double> >& data,
			  std::vector<unsigned>& mult,
			  bool linear=false);
  
  void LebesgueList(BSPTree::indx root, unsigned dindx, unsigned long ncut, 
		    std::vector< std::vector<double> >& data, 
		    std::vector<unsigned>& mult,
		    std::map<int, std::multimap<double, LeafData> >& lebesgueLeaves,
		    bool linear);
  //@}
  
  #ifndef SWIG
    // AUTO GENERATED BY ../persistence/autopersist.py
    private:
    friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive & ar, const unsigned int version) {
        this->pre_serialize(ar, version);
         try {                                                         
          ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(BSPTree);            
          BIE_CATCH_BOOST_SERIALIZATION_EXCEPTION;                     
         }                                                             
         try {                                                         
          ar & BOOST_SERIALIZATION_NVP(verbose_debug);                        
          BIE_CATCH_BOOST_SERIALIZATION_EXCEPTION;                     
         }                                                             
         try {                                                         
          ar & BOOST_SERIALIZATION_NVP(use_range);                        
          BIE_CATCH_BOOST_SERIALIZATION_EXCEPTION;                     
         }                                                             
         try {                                                         
          ar & BOOST_SERIALIZATION_NVP(gen);                        
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
          ar & BOOST_SERIALIZATION_NVP(stally);                        
          BIE_CATCH_BOOST_SERIALIZATION_EXCEPTION;                     
         }                                                             
         try {                                                         
          ar & BOOST_SERIALIZATION_NVP(centers);                        
          BIE_CATCH_BOOST_SERIALIZATION_EXCEPTION;                     
         }                                                             
         try {                                                         
          ar & BOOST_SERIALIZATION_NVP(center0);                        
          BIE_CATCH_BOOST_SERIALIZATION_EXCEPTION;                     
         }                                                             
         try {                                                         
          ar & BOOST_SERIALIZATION_NVP(minS);                        
          BIE_CATCH_BOOST_SERIALIZATION_EXCEPTION;                     
         }                                                             
         try {                                                         
          ar & BOOST_SERIALIZATION_NVP(maxS);                        
          BIE_CATCH_BOOST_SERIALIZATION_EXCEPTION;                     
         }                                                             
         try {                                                         
          ar & BOOST_SERIALIZATION_NVP(ranges);                        
          BIE_CATCH_BOOST_SERIALIZATION_EXCEPTION;                     
         }                                                             
         try {                                                         
          ar & BOOST_SERIALIZATION_NVP(range0);                        
          BIE_CATCH_BOOST_SERIALIZATION_EXCEPTION;                     
         }                                                             
         try {                                                         
          ar & BOOST_SERIALIZATION_NVP(Ledge);                        
          BIE_CATCH_BOOST_SERIALIZATION_EXCEPTION;                     
         }                                                             
         try {                                                         
          ar & BOOST_SERIALIZATION_NVP(Hedge);                        
          BIE_CATCH_BOOST_SERIALIZATION_EXCEPTION;                     
         }                                                             
         try {                                                         
          ar & BOOST_SERIALIZATION_NVP(weights);                        
          BIE_CATCH_BOOST_SERIALIZATION_EXCEPTION;                     
         }                                                             
         try {                                                         
          ar & BOOST_SERIALIZATION_NVP(beta);                        
          BIE_CATCH_BOOST_SERIALIZATION_EXCEPTION;                     
         }                                                             
         try {                                                         
          ar & BOOST_SERIALIZATION_NVP(lower);                        
          BIE_CATCH_BOOST_SERIALIZATION_EXCEPTION;                     
         }                                                             
         try {                                                         
          ar & BOOST_SERIALIZATION_NVP(upper);                        
          BIE_CATCH_BOOST_SERIALIZATION_EXCEPTION;                     
         }                                                             
         try {                                                         
          ar & BOOST_SERIALIZATION_NVP(left_child);                        
          BIE_CATCH_BOOST_SERIALIZATION_EXCEPTION;                     
         }                                                             
         try {                                                         
          ar & BOOST_SERIALIZATION_NVP(right_child);                        
          BIE_CATCH_BOOST_SERIALIZATION_EXCEPTION;                     
         }                                                             
         try {                                                         
          ar & BOOST_SERIALIZATION_NVP(lo_leaf);                        
          BIE_CATCH_BOOST_SERIALIZATION_EXCEPTION;                     
         }                                                             
         try {                                                         
          ar & BOOST_SERIALIZATION_NVP(hi_leaf);                        
          BIE_CATCH_BOOST_SERIALIZATION_EXCEPTION;                     
         }                                                             
         try {                                                         
          ar & BOOST_SERIALIZATION_NVP(permutation);                        
          BIE_CATCH_BOOST_SERIALIZATION_EXCEPTION;                     
         }                                                             
         try {                                                         
          ar & BOOST_SERIALIZATION_NVP(next);                        
          BIE_CATCH_BOOST_SERIALIZATION_EXCEPTION;                     
         }                                                             
         try {                                                         
          ar & BOOST_SERIALIZATION_NVP(minL);                        
          BIE_CATCH_BOOST_SERIALIZATION_EXCEPTION;                     
         }                                                             
         try {                                                         
          ar & BOOST_SERIALIZATION_NVP(maxL);                        
          BIE_CATCH_BOOST_SERIALIZATION_EXCEPTION;                     
         }                                                             
         try {                                                         
          ar & BOOST_SERIALIZATION_NVP(measureList);                        
          BIE_CATCH_BOOST_SERIALIZATION_EXCEPTION;                     
         }                                                             
         try {                                                         
          ar & BOOST_SERIALIZATION_NVP(ncutM);                        
          BIE_CATCH_BOOST_SERIALIZATION_EXCEPTION;                     
         }                                                             
         try {                                                         
          ar & BOOST_SERIALIZATION_NVP(measureMap);                        
          BIE_CATCH_BOOST_SERIALIZATION_EXCEPTION;                     
         }                                                             
         try {                                                         
          ar & BOOST_SERIALIZATION_NVP(dnorm);                        
          BIE_CATCH_BOOST_SERIALIZATION_EXCEPTION;                     
         }                                                             
         try {                                                         
          ar & BOOST_SERIALIZATION_NVP(temp);                        
          BIE_CATCH_BOOST_SERIALIZATION_EXCEPTION;                     
         }                                                             
        this->post_serialize(ar, version);
    }
    #endif

};

/**
   Helper class for nearest neighbor sorting
*/
class Delta: public Serializable
{
private:

  double minimum, maximum;
  KDTree::indx Key;

  double minDist(KDTree::indx key, const double*, const KDTree* kd) const;
  double maxDist(KDTree::indx key, const double*, const KDTree* kd) const;

public:
  
  //! Null consturctor
  Delta() : minimum(0), maximum(0), Key(0) {}
  
  //! Main constructor
  Delta(KDTree::indx key, const double* p, const KDTree* q) : Key(key)
  {
    minimum = minDist(key, p, q);
    maximum = maxDist(key, p, q);
  }
  
  //! Copy constructor
  Delta(const Delta& p)
  {
    minimum = p.minimum;
    maximum = p.maximum;
    Key     = p.Key;
  }
  
  KDTree::indx key() const { return Key; }

  double Min() const { return minimum; }
  double Max() const { return maximum; }

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
          ar & BOOST_SERIALIZATION_NVP(minimum);                        
          BIE_CATCH_BOOST_SERIALIZATION_EXCEPTION;                     
         }                                                             
         try {                                                         
          ar & BOOST_SERIALIZATION_NVP(maximum);                        
          BIE_CATCH_BOOST_SERIALIZATION_EXCEPTION;                     
         }                                                             
         try {                                                         
          ar & BOOST_SERIALIZATION_NVP(Key);                        
          BIE_CATCH_BOOST_SERIALIZATION_EXCEPTION;                     
         }                                                             
        this->post_serialize(ar, version);
    }
    #endif

};


BOOST_SERIALIZATION_SHARED_PTR(LeafElem)
BOOST_SERIALIZATION_SHARED_PTR(LeafData)
BOOST_SERIALIZATION_SHARED_PTR(KDTree)
BOOST_SERIALIZATION_SHARED_PTR(Delta)

#ifndef SWIG
BIE_CLASS_TYPE_INFO(LeafElem)
BIE_CLASS_TYPE_INFO(CellElem)
BIE_CLASS_TYPE_INFO(KDTree)
BIE_CLASS_TYPE_INFO(Delta)
BIE_CLASS_EXPORT_KEY(LeafElem)
BIE_CLASS_EXPORT_KEY(CellElem)
BIE_CLASS_EXPORT_KEY(KDTree)
BIE_CLASS_EXPORT_KEY(Delta)
#endif
#endif