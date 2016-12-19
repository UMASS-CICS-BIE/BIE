// This is really -*- C++ -*-

#ifndef CLIVECTOR_H
#define CLIVECTOR_H

#include <iostream>
#include <vector>
#include <string>
#include <typeinfo>

#include <Serializable.h>

namespace BIE {

  class Tessellation;
  class Distribution;
  class Simulation;
  class RJMapping;

  /** 
      A template class to wrap the STL vector and handle CLI member
      functions for accessing the elements.  The template contains the
      serizliation code and the CLI interface member functions.  The
      () operator returns the managed vector<T> instance for the
      programmer.  In other words, all the STL functions are available
      for instance <code>clivectord v</code> by using <code>v()</code>.

      An earlier implementation used a typedef for the STL vector<T>
      directly.  There were two disadvantages:
      <ol>
      <li> We had to use a manipulator "helper" class to defined elements 
      of the vector<T> from within the CLI
      <li> The persistence system had to treat each of these as a special 
      case in the CLI symbol table.
      </ol>
      The wrapper allows the standard upcasting restoration method to
      work for STL types and simplifies the CLI manipulation at the
      expense of a bit more work for the programmer.  Seems the like
      the better tradeoff, to me.
  */
  template<class T>
  class clivector : public Serializable
  {
  private:

    std::vector<T> v;

  public:

    //@{
    //! Constructors
    clivector<T>()           {                      }
    clivector<T>(int n)      { v = std::vector<T>(n   ); }
    clivector<T>(int n, T z) { v = std::vector<T>(n, z); }
    //! Copy a STL vector
    clivector<T>(const std::vector<T>& z) { v = z; }
    //@}

    //! Prepend an STL vector
    void prepend(clivector<T>* a)
    { v.insert(v.begin(), a->v.begin(), a->v.end()); }

    //! Append an STL vector
    void append(clivector<T>* a)
    { v.insert(v.end(), a->v.begin(), a->v.end()); }

    //! Prepend a value
    void prepend(T a)
    { v.insert(v.begin(), a); }

    //! Append a value
    void append(T a)        
    { v.insert(v.end(), a);   }

    //! Set vector elements
    void setval(int n, T z) { 
      if (n < static_cast<int>(v.size())) v[n] = z;
      else std::cout << "Index out of bounds"
		     << " (must be < " << v.size() << ")" << std::endl;
    }
  
    //! Get the vector
    std::vector<T>& operator()() { return v; }

    //! Print vector elements
    void showval(int n) 
    {
      if (n < static_cast<int>(v.size())) std::cout << v[n] << std::endl;
      else std::cout << "Index out of bounds"
		     << " (must be < " << v.size() << ")" << std::endl;
    }
    
    //! Print all vector elements
    void showall() 
    {
      if (v.size()) {
	for (unsigned i=0; i<v.size(); i++)
	  std::cout << "v[" << i << "] = " << v[i] << std::endl;
      } else {
	std::cout << "empty" << std::endl;
      }
    }
    

  private:

    friend class boost::serialization::access;

    template<class Archive>
    void serialize(Archive &ar, const unsigned int file_version) {
      try {
	ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(Serializable);
	BIE_CATCH_BOOST_SERIALIZATION_EXCEPTION;
      }
      try {
	ar & BOOST_SERIALIZATION_NVP(v);
	BIE_CATCH_BOOST_SERIALIZATION_EXCEPTION;
      }
    }
  };

  /*
    These classes are simply typedefs to expanded templates.  The
    stanzas below define the typedef-named classes to the reflection
    system. Too bad the template can't expand the stylized comments as
    well!
  */

  //! A wrapper for vector<double>
  typedef clivector<double> clivectord;

  //+ CLICLASS clivectord
  //+ CLICONSTR
  //+ CLICONSTR int
  //+ CLICONSTR int double
  //+ CLIMETHOD void append  double
  //+ CLIMETHOD void append  clivectord*
  //+ CLIMETHOD void prepend double
  //+ CLIMETHOD void prepend clivectord*
  //+ CLIMETHOD void setval  int double
  //+ CLIMETHOD void showval int
  //+ CLIMETHOD void showall
  
  //! A wrapper for vector<int>
  typedef clivector<int> clivectori;

  //+ CLICLASS clivectori
  //+ CLICONSTR
  //+ CLICONSTR int
  //+ CLICONSTR int int
  //+ CLIMETHOD void append  int
  //+ CLIMETHOD void append  clivectori*
  //+ CLIMETHOD void prepend int
  //+ CLIMETHOD void prepend clivectori*
  //+ CLIMETHOD void setval  int int
  //+ CLIMETHOD void showval int
  //+ CLIMETHOD void showall
  
  //! A wrapper for vector<string>
  typedef clivector<std::string> clivectors;

  //+ CLICLASS clivectors
  //+ CLICONSTR
  //+ CLICONSTR int
  //+ CLICONSTR int string
  //+ CLIMETHOD void append  string
  //+ CLIMETHOD void append  clivectors*
  //+ CLIMETHOD void prepend string
  //+ CLIMETHOD void prepend clivectors*
  //+ CLIMETHOD void setval  int string
  //+ CLIMETHOD void showval int
  //+ CLIMETHOD void showall
  
  //! A wrapper for vector<Tessellation*>
  typedef clivector<Tessellation*> clivectortess;

  //+ CLICLASS clivectortess
  //+ CLICONSTR
  //+ CLICONSTR int
  //+ CLICONSTR int Tessellation*
  //+ CLIMETHOD void append  Tessellation*
  //+ CLIMETHOD void append  clivectortess*
  //+ CLIMETHOD void prepend Tessellation*
  //+ CLIMETHOD void prepend clivectortess*
  //+ CLIMETHOD void setval  int Tessellation*
  //+ CLIMETHOD void showval int
  //+ CLIMETHOD void showall
  
  //! A wrapper for vector<Distribution*>
  typedef clivector<Distribution*> clivectordist;

  //+ CLICLASS clivectordist
  //+ CLICONSTR
  //+ CLICONSTR int
  //+ CLICONSTR int Distribution*
  //+ CLIMETHOD void append  Distribution*
  //+ CLIMETHOD void append  clivectordist*
  //+ CLIMETHOD void prepend Distribution*
  //+ CLIMETHOD void prepend clivectordist*
  //+ CLIMETHOD void setval  int Distribution*
  //+ CLIMETHOD void showval int
  //+ CLIMETHOD void showall

  //! A wrapper for vector<Simulation*>
  typedef clivector<Simulation*> clivectorsim;

  //+ CLICLASS clivectorsim
  //+ CLICONSTR
  //+ CLICONSTR int
  //+ CLICONSTR int Simulation*
  //+ CLIMETHOD void append  Simulation*
  //+ CLIMETHOD void append  clivectorsim*
  //+ CLIMETHOD void prepend Simulation*
  //+ CLIMETHOD void prepend clivectorsim*
  //+ CLIMETHOD void setval  int Simulation*
  //+ CLIMETHOD void showval int
  //+ CLIMETHOD void showall

  //! A wrapper for vector<RJMapping*>
  typedef clivector<RJMapping*> clivectorRJmap;

  //+ CLICLASS clivectorRJmap
  //+ CLICONSTR
  //+ CLICONSTR int
  //+ CLICONSTR int RJMapping*
  //+ CLIMETHOD void append  RJMapping*
  //+ CLIMETHOD void append  clivectorRJmap*
  //+ CLIMETHOD void prepend RJMapping*
  //+ CLIMETHOD void prepend clivectorRJmap*
  //+ CLIMETHOD void setval  int RJMapping*
  //+ CLIMETHOD void showval int
  //+ CLIMETHOD void showall

}
  
#ifndef SWIG
  BIE_CLASS_TYPE_INFO (BIE::clivectord)
  BIE_CLASS_EXPORT_KEY(BIE::clivectord)

  BIE_CLASS_TYPE_INFO (BIE::clivectori)
  BIE_CLASS_EXPORT_KEY(BIE::clivectori)

  BIE_CLASS_TYPE_INFO (BIE::clivectors)
  BIE_CLASS_EXPORT_KEY(BIE::clivectors)

  BIE_CLASS_TYPE_INFO (BIE::clivectortess)
  BIE_CLASS_EXPORT_KEY(BIE::clivectortess)

  BIE_CLASS_TYPE_INFO (BIE::clivectordist)
  BIE_CLASS_EXPORT_KEY(BIE::clivectordist)

  BIE_CLASS_TYPE_INFO (BIE::clivectorsim)
  BIE_CLASS_EXPORT_KEY(BIE::clivectorsim)

  BIE_CLASS_TYPE_INFO (BIE::clivectorRJmap)
  BIE_CLASS_EXPORT_KEY(BIE::clivectorRJmap)
#endif


#endif
