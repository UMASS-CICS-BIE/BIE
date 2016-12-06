// This is really -*- C++ -*-

#ifndef DataAnyDim_h
#define DataAnyDim_h

#include <vector>

namespace BIE {

  /**
     Defines a single data point, now with multiple attributes
  */
  //+ CLICLASS Data
  class Data {

  private:
    double X, Y, wght;
    int len;
    vector<double> attrib;

  public:
    //+ CLICONSTR int
    //! Constructor
    Data(int length) {
      len = length;
      wght = 1.0;
      attrib.resize(len);
      for (int i=0; i<len; i++) 
	attrib[i] = 0.0;
    }
        
    //! Overloaded [] to get attributes of data point
    //@{
    //! Assignable
    double &operator[] (int i) { return attrib[i]; }
    //! Constant
    double operator[] (int i) const { return attrib[i]; }
    //@}


    //! Get X value of data
    //@{
    //! Assignable
    double &x() { return X; }
    //! Constant
    double x() const { return X; }
    //@}

    //! Get X value of data
    //@{
    //! Assignable
    double &y() { return Y; }
    //! Constant
    double y() const { return Y; }
    //@}

    //! Get weight of data point
    double &weight() { return wght; }

    //! Resize data
    void resize(int size) { attrib.resize(size); }  

    //! Get attribute data
    vector<double> &attribute() { return attrib;}

    //! Set the x-value
    void set_x(double x_in) {X = x_in;}
    //! Set the y-value
    void set_y(double y_in) {Y = y_in;}
    //! Set the weight value
    void set_w(double w_in) {wght = w_in;}
 
    //! Set the attribute at at @param position to @param value
    void set_attrib(int position, double value) { attrib[position] = value; }

  };
}

#endif


