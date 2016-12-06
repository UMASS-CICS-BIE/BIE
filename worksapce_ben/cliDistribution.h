#ifndef CLIDISTRIBUTION_H
#define CLIDISTRIBUTION_H

#include <vector>
#include "Distribution.h"

using namespace std;

namespace BIE {
  
  //+ CLICLASS vector<BinnedDistribution*> [1]
  //+ CLICONSTR int

  /// Wrapper class for vector<BinnedDistribution*> for cli instantiation
  //+ CLICLASS cliDistribution SUPER vector<BinnedDistribution*>
  class cliDistribution : public vector<BinnedDistribution*>
  {
  public:
    //+ CLICONSTR
    //! Constructor
    cliDistribution() : vector<BinnedDistribution*>() {
    }

    //+ CLICONSTR void*
    //! Constructor
    cliDistribution(void *) {};

    //+ CLICONSTR int
    //! an array of pointers to BinnedDistributions
    cliDistribution(int i) : vector<BinnedDistribution*>(i) {
    }

    //+ CLIMETHOD BinnedDistribution* getval int
    //! get a value from the distribution list
    BinnedDistribution* getval(int pos) { 
      return  (*this)[pos];
    }

    //+ CLIMETHOD void setval int BinnedDistribution* 
    //! Assign a particular distribution
    void setval(int i, BinnedDistribution* val) {
      (*this)[i] = val;
    }    
    
  };  
}
  
#endif
