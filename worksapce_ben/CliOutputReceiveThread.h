#ifndef CliOutputReceiveThread_h
#define CliOutputReceiveThread_h

#include <string>
#include <cc++/thread.h>
#include <mpi.h>
#include <sstream>

#ifdef  CCXX_NAMESPACES
 using namespace std;
 using namespace ost;
#endif

namespace BIE { 

  //! Receive data through pipe from running simulation
  class CliOutputReceiveThread: public Thread {
  public:
    //! Constructor
    CliOutputReceiveThread(string hostname, int port, Semaphore *consumer);
    //! Constructor
    CliOutputReceiveThread(int socketid, Semaphore *consumer);
    //! Destructor
    ~CliOutputReceiveThread();

    //! Start the thread
    void run();
    //! Get the stream with the CLI output
    istream *getReply();
    //! Return my controller semaphore
    Semaphore *getSemaphore();
    //! Set my controller semaphore
    Semaphore *setConsumerSemaphore(Semaphore *consumer) {return (_consumer = consumer);};
    //! I have data
    bool isReady(){return _isReady;};
    //! Exit as soon as possible
    void setStop() {_stop=true;}

  private:
    int readReply();
    void initialize(Semaphore *consumer);
    
    char *_hname;
    int _port;
    int _socketid;
    char *_receiveBuffer;
    int _receiveBufferSize;
    Semaphore *_qsem;
    Semaphore *_consumer;
    bool _stop, _isReady;
    int _size;
  };
}
#endif
