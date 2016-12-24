// This is really -*- C++  -*-

#ifndef Serializable_h
#define Serializable_h

#include <iostream>
#include <typeinfo>

#include <boost/version.hpp>

//
// Archive headers
//
#include <boost/archive/archive_exception.hpp>

#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <boost/archive/xml_oarchive.hpp>
#include <boost/archive/xml_iarchive.hpp>
#include <boost/archive/impl/archive_serializer_map.ipp>

//
// Serialization
//
#include <boost/serialization/export.hpp>
#include <boost/serialization/serialization.hpp>
#include <boost/serialization/split_free.hpp>
#include <boost/serialization/shared_ptr.hpp>
#include <boost/serialization/assume_abstract.hpp>

//
// STL serialization headers
//
#include <boost/serialization/deque.hpp>
#include <boost/serialization/list.hpp>
#include <boost/serialization/set.hpp>
#include <boost/serialization/vector.hpp>
#include <boost/serialization/string.hpp>
#include <boost/serialization/map.hpp>

#include "PersistenceException.h"

namespace BIE {

  //! The is the base class for a serializable class using boost
  class Serializable 
  { 
  public:
    //! Default constructor
    Serializable();

    //! Destructor
    virtual ~Serializable() {}

  protected:
    friend class boost::serialization::access;

    //! This must be overridden by derived classes . . .
    template<class Archive>
      void serialize(Archive &ar, const unsigned int file_version) {
    }
    

    /** These functions are to be overridded by classes that need special
	treatment for serialization, saving or loading.  They do nothing
	by default */
    //@{

    /** These may be used when the save/load function is NOT split, in
	other words, one may use the same calls on saving and loading */
    //@{
    //! Called BEFORE serializing a class 
    template<class Archive>
      void pre_serialize(Archive &ar, const unsigned int file_version) {}
    
    //! Called AFTER serializing a class 
    template<class Archive>
      void post_serialize(Archive &ar, const unsigned int file_version) {}
    //@}

    /**  These may be used when the save/load function is IS split, in
	 other words, one must use different calls for saving and loading */
    //@{
    //! Called BEFORE saving a class 
    template<class Archive>
      void pre_save(Archive &ar, const unsigned int file_version) const {}

    //! Called AFTER saving a class 
    template<class Archive>
      void post_save(Archive &ar, const unsigned int file_version) const {}

    //! Called BEFORE loading a class 
    template<class Archive>
      void pre_load(Archive &ar, const unsigned int file_version) {}

    //! Called AFTER loading a class 
    template<class Archive>
      void post_load(Archive &ar, const unsigned int file_version) {}
    //@}
    //@}

  };

} // namespace BIE


//
// For convenience
//
#ifndef BIE_CLASS_TYPE_INFO
#define BIE_CLASS_TYPE_INFO(type) \
  BOOST_CLASS_TYPE_INFO(::type, extended_type_info_typeid< ::type>)
#endif

//
// Deal with BOOST version interface changes
//
#ifndef SWIG

BOOST_SERIALIZATION_ASSUME_ABSTRACT(BIE::Serializable);

#define BIE_CLASS_EXPORT_KEY(type) \
BOOST_CLASS_EXPORT_KEY(type)

#define BIE_CLASS_EXPORT_IMPLEMENT(type) \
BOOST_CLASS_EXPORT_IMPLEMENT(type)

#define BIE_CLASS_ABSTRACT(type) \
BOOST_SERIALIZATION_ASSUME_ABSTRACT(type)

#define BIE_CLASS_EXPORT_KEY(type) \
BOOST_CLASS_EXPORT_KEY(type)
#define BIE_CLASS_EXPORT_IMPLEMENT(type) \
BOOST_CLASS_EXPORT_IMPLEMENT(type)

//
// Register this class . . .
//
BIE_CLASS_TYPE_INFO(BIE::Serializable)
BOOST_CLASS_EXPORT_KEY(BIE::Serializable)

//
// Nota Bene:
//
// By default, object tracking is performed if and only if an object
// of the the class is anywhere serialized through a pointer.  This
// almost certanly occurs in BIE, but if it doesn't, the shared base
// object wouldn't be tracked.  This would lead to multiple save/load
// operation of the data in this shared base class, and it would be a
// big waste of time.  So the tracking behavior trait of the base
// class is set to always track serialized objects for this base
// class.  This permits the system to detect and elminate redundent
// save/load operations.
//
BOOST_CLASS_TRACKING(BIE::Serializable, boost::serialization::track_always)

//
// Series of catch clauses needed by the persistence stack trace (for
// debugging)
//
#ifndef BIE_CATCH_BOOST_SERIALIZATION_EXCEPTION
#define BIE_CATCH_BOOST_SERIALIZATION_EXCEPTION				\
  } catch (::boost::archive::archive_exception & e) {			\
    throw new BoostSerializationException(string(e.what()), 		\
					  __FILE__,__LINE__);		\
  } catch (BoostSerializationException * e) {				\
    throw new BoostSerializationException(e,__FILE__,__LINE__); 
#endif

#endif

#endif