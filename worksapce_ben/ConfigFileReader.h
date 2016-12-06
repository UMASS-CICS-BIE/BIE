// This is really -*- C++ -*-

#ifndef ConfigFileReader_h
#define ConfigFileReader_h

#include <string>
using namespace std;
#include <cc++/misc.h>

using namespace ost;

/**
 * This class is a wrapper to the Common C++ class used for INI file
 * processing and provides a few additional conveniences for config
 * file management.  The class provides methods that give typed access
 * to the variables.  Also, undefined variables are handled a bit more
 * gracefully than in the provided interface.  If a variable is not
 * found, a suitably typed default return value can be provided rather
 * than NULL.  The class will also take care of creating a .bierc file
 * for a user that doesn't have one already.
 *
 * Typical usage:
 * ConfigFileReader * kd = new ConfigFileReader("cli");
 * bool mybool = kd->getValue("mybool", true);
 * string mystring = kd->getValue("mystring", "Alistair");
 */
class ConfigFileReader : public Keydata
{
public:
  /// Constructs an object for reading values of variables from INI files.
  /// The section of the INI file must be specified (e.g. "cli" is a section
  /// in the .bierc config file).  Will also create a ~/.bierc file
  /// for users who are missing this file.
  ///
  /// section [] specifies the file section we want to refer to 
  /// (e.g. "cli" refers to ~username/.bierc section [cli] followed by
  /// the [cli] stanza in a .bierc defined in the working dir, if it exists)
  ConfigFileReader(const string& keyword);
  
  /// Get the last set value for a configuration variable as a string.
  /// The value will probably come from a config file.  When the variable
  /// does not exist, the default value is returned.
  const char* getValue(const char *symbol, const char * defaultvalue);
  
  /// Get the last set value for a configuration variable as a boolean.
  /// The value will probably come from a config file.  When the variable
  /// does not exist, the default value is returned.
  bool getValue(const char *symbol, bool defaultvalue);
  
  /// Get the last set value for a configuration variable as an int.
  /// The value will probably come from a config file.  When the variable
  /// does not exist, the default value is returned.
  int getValue(const char *symbol, int defaultvalue);
  
  /// Get the last set value for a configuration variable as a double.
  /// The value will probably come from a config file.  When the variable
  /// does not exist, the default value is returned.
  double getValue(const char *symbol, double defaultvalue);
  
  /// Try to find BIE source
  std::string getBIEDIR(const string&);

private:
  /// Creates a configuration file at .bierc for users who do not have this
  /// already.
  void createConfigFile(const string&);
};

#endif
