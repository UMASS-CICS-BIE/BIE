// This is really -*- C++ -*-

#ifndef EvaluationRequest_h
#define EvaluationRequest_h

#include <BIEmpi.h>
#include <cc++/thread.h>
#include <stdio.h>  // for sprintf
#include <TSLog.h>
#include "ElapsedTime.h"
#include <vector>
using namespace std;
using namespace ost;

namespace BIE {

  //! Request model evaluation
  class EvaluationRequest {
  public:
    /// Constructor
    EvaluationRequest() {
       sent = false;
       received = false;
    }
    /// Destructor
    virtual ~EvaluationRequest(){}

    /// Send to node via a list of communicators
    virtual void Send(int destination, vector<MPI_Comm>* comms)=0;

    /// Receive a message
    virtual void Recv(MPI_Status *status)=0;

    /// Wait for MPI completion
    virtual void WaitForCompletion()
      { finished->wait();}


    /// The "quick and dirty" log instance
    TSLog log;

    /// Write a message to the log
    void LogIt(const char* buf) { log.LogIt(buf);}
    /// Write a message to the log with an integer value
    void LogInt(const char* buf, int i) { log.LogInt(buf, i);}
    /// Write a message to the log with an double value
    void LogDouble(const char* buf, double d) { log.LogDouble(buf, d);}
    /// Print the log
    void PrintLog() { log.PrintIt();}

    /// even quicker and dirtier timing
    //@{
    //! Time in mutex waits
    ElapsedTime sendMutexWaitTime;
    //! Time in sending commands
    ElapsedTime sendCmdTime;
    //! Time evaluating something useful
    ElapsedTime sendEvalTime;
    //@}

    //! True if value sent
    bool sent;
    //! True if value received
    bool received;
    //! Caller passes in a semphore to be signaled when request is complete
    Semaphore *finished;

  protected:

  };
}
#endif // EvaluationRequest_h
