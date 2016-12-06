#ifndef GETFILE_h
#define GETFILE_h

#include <stdlib.h>
#include <stdio.h>
#include <sys/types.h>
#include <sys/socket.h>
#include <netinet/in.h>
#include <string.h>
#include <sys/socket.h>

namespace BIE {

#define MAXBUFSIZE 512
#define MAX_REQUEST_LEN 1024

  //@{
  //! Functions for tesstool
  void graphrequest_proc(int new_id);
  int sendmsg(int socket, char *msg, int size);
  //@}

}
#endif

