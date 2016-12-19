// This is really -*- C++ -*-

#ifndef Frontier_h
#define Frontier_h

#include <string>
#include <algorithm>

#include "Serializable.h"


#include <BIEmpi.h>
#include <Node.h>
#include <Tessellation.h>
#include <FrontierExpansionHeuristic.h>

namespace BIE {

  class FrontierExpansionHeuristic;

  //+ CLICLASS Frontier
  //! Used to represent a frontier in a tessellation tree.  which can be
  //! mutated based on user supplied functions or heuristics
  //! 
  //! There are two frontier modes: accumulated mode and normal mode.
  //! In accumulated mode, the active tree is the entire tree upto and
  //! including the frontier.  This choice is relevant for point sets
  //! with no clear spatial extent.  In normal mode, the active tree is
  //! the frontier itself.  Normally, the tessellation is exclusive so
  //! the tesselation denoted by the frontier is the entire space
  class Frontier: public Serializable 
  {
    
  public:
    
    //! Constructs a new Frontier object for the given tessellation.
    //! The frontier is initially set to be the set of root nodes in the
    //! tessellation.
    Frontier(Tessellation * tess);
    
    //! Copy constructor
    Frontier(const Frontier*);
    
    //+ CLIMETHOD Frontier* Copy
    //! Make a copy of myself
    Frontier* Copy();
    
    //! Destructor
    /* virtual */ ~Frontier() {};
    
    //! Returns the tessellation the frontier object applies to.
    //! Each frontier is associated with one (and only one) tessellation.
    inline Tessellation * getTessellation() 
    { return frontier_tessellation; }
    
    //! Returns a vector containing tileids of the tiles in the frontier.
    inline vector<int> ExportFrontier() { return frontier_tiles; }
    
    //+ CLIMETHOD void RetractToTopLevel
    //! Retract the frontier to the root nodes.
    void RetractToTopLevel();
    
    //+ CLIMETHOD void UpDownLevels int
    /** 
	Expand or contract frontier by the specified number of levels.
	\param n Number of levels to expand (positive N) or contract (negative N).
    */
    void UpDownLevels(int n);
    
    /**
     * Sets the nodes/tiles in the frontier explicitly.  
     * Nodes that are decended from another node in the frontier state are 
     * excluded (i.e. the most general tiles get preference).  
     * This doesn't just naively set the frontier since 
     * we risk invalid frontiers (e.g. a parent and child tile both in 
     * frontier), but recurses through the tree (so can be expensive).
     * \param FrontierState The new frontier specified as a vector of tile IDs.
     */
    void Set(vector<int> FrontierState);
    
    //+ CLIMETHOD bool IncreaseResolution FrontierExpansionHeuristic* int
    /**
       Increases the resolution of the frontier where recommended by a
       heuristic.  See the description of FrontierExpansionHeuristic and
       its subclasses for details of possible heuristics.
       
       \param heuristic An expansion heuristic.
       \param numlevels Number of levels to skip during expansion.
    */
    //+ CLIMETHOD bool IncreaseResolution FrontierExpansionHeuristic*
    /**
       Increases the resolution of the frontier where recommended by a heuristic.
       See the description of FrontierExpansionHeuristic and its subclasses
       or details of possible heuristics.
       
       param heuristic An expansion heuristic.
    */
    bool IncreaseResolution (FrontierExpansionHeuristic * heuristic, int numlevels = 1);
    
    //! Find the frontier node with the given coordinates,
    //! or return NULL if the coordinates are not covered by the frontier.
    Node * Find(double x, double y);
    
    //! Frontier iterator function: returns the tileid of the first tile/node
    //! in the frontier.
    inline int First(void) { return frontier_tiles[0]; }
    
    //! Frontier iterator function: returns the tileid of the last tile/node
    //! in the frontier.
    inline int Last(void) { return frontier_tiles[frontier_tiles.size()-1]; }
    
    //! Frontier iterator function: returns the tileid of the next tile/node
    //! in the frontier, and advances the position of the cursor to this tile.
    //! Returns the last tileid continually when the end is reached.
    int Next(void);
    
    //! Frontier iterator function: returns the tileid of the tile/node
    //! currently pointed to by the cursor.
    int CurrentItem(void);
    
    //! Frontier iterator function: Resets the cursor to the first tile/node
    //! in the frontier.
    inline void Reset(void) { frontier_index = 0; }
    
    //! Returns true once Next() has returned the final tile/node in the
    //! frontier.
    inline bool IsDone(void) { return (frontier_index == frontier_tiles.size()); }
    
    //! Frontier size
    unsigned int Size() { return frontier_tiles.size(); }
    
    //+ CLIMETHOD void printSize
    //! Print the frontier size to the console
    void printSize() { 
      if (myid==0) cout << "Frontier size=" << frontier_tiles.size() << endl;
    }
    
    //! Is frontier in accumulate mode?
    inline bool AccumulateMode() { return accumulate_mode; }
    
  private:
    //! Frontier state expressed as a list of tileids.
    vector<int> frontier_tiles;
    
    //! Frontier state expressed as a vector of bools, where the field
    //! index by a tileid is true if the tile is in the frontier.
    vector<bool> frontier_infrontier;
    
    //! Frontier index     
    unsigned int frontier_index;
    
    //! Reference to the tessellation this frontier applies to.
    Tessellation * frontier_tessellation;
    
    //! Assume an accumulating frontier: that is, everything at lower
    //! levels remains on the frontier.  Currently, only PointTessellations
    //! use this mode.
    bool accumulate_mode;
    
    //! Increase resolution one level.
    void ExpandOneLevel();
    
    //! Decrease resolution one level.
    void ContractOneLevel();
    
    //! Returns a node object given a tile id.
    Node * GetNode(int tileid);
    
    //! Internal function for erasing the frontier.
    void EraseFrontier();
    
    //! Internal function for adding a node to the frontier.
    //! It assumes caller knows what it is doing.
    void AddToFrontier(int tileid);
    
    //! Look for new frontier recursively in accumulation mode
    bool WalkAndAdd(vector<bool>& frontiercopy, int& parentid);
    
    
    //! Find the true frontier when in accumulation mode
    //@{
    //! Top level member
    vector<int> TrueFrontier();
    //! Recursion member
    void TrueFrontier(vector<int>& real_frontier, int& parentid);
    //@}
    
    //! Recursively select frontier by walking tessllation tree
    void RecursivelySelect(Node *node, vector<bool> & nodesetbitvector);
    
    //! Finds the frontier node containing the coordinates descended 
    //! from the node passed in, or returns NULL if no such node exists.
    Node * FindInTree(Node* node, double x, double y);
    
    #ifndef SWIG
    // AUTO GENERATED BY ../persistence/autopersist.py
    protected:
    Frontier() {}
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
          ar & BOOST_SERIALIZATION_NVP(frontier_tiles);                        
          BIE_CATCH_BOOST_SERIALIZATION_EXCEPTION;                     
         }                                                             
         try {                                                         
          ar & BOOST_SERIALIZATION_NVP(frontier_infrontier);                        
          BIE_CATCH_BOOST_SERIALIZATION_EXCEPTION;                     
         }                                                             
         try {                                                         
          ar & BOOST_SERIALIZATION_NVP(frontier_index);                        
          BIE_CATCH_BOOST_SERIALIZATION_EXCEPTION;                     
         }                                                             
         try {                                                         
          ar & BOOST_SERIALIZATION_NVP(frontier_tessellation);                        
          BIE_CATCH_BOOST_SERIALIZATION_EXCEPTION;                     
         }                                                             
         try {                                                         
          ar & BOOST_SERIALIZATION_NVP(accumulate_mode);                        
          BIE_CATCH_BOOST_SERIALIZATION_EXCEPTION;                     
         }                                                             
        this->post_serialize(ar, version);
    }
    #endif

      };
  
} // namespace BIE

#ifndef SWIG
BIE_CLASS_TYPE_INFO(BIE::Frontier)
BIE_CLASS_EXPORT_KEY(BIE::Frontier)
#endif
#endif
