#ifndef LikelihoodComputation_h
#define LikelihoodComputation_h

namespace BIE {

  /**
     Abstract base class for Likelihood computation.
     This class also supplies a default serial implementation.

  */
  class LikelihoodComputation {
  public:
    
    ///
    LikelihoodComputation();


    /// create virtual destructor for subclasses
    virtual ~LikelihoodComputation() {}

    /** Determine liklihood for given state.
	usually overriden in derived class.

	@param state is the parameter vector and meta info
	@param indx  is the "temperature" level
    */
    virtual double Likelihood(State &state, int indx);

  };

}
#endif // LikelihoodComputation_h
