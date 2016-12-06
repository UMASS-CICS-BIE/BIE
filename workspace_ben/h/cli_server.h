#ifndef CLI_SERVER_h
#define CLI_SERVER_h

#include <sys/types.h>
#include <sys/socket.h>
#include <netinet/in.h>
#include <arpa/inet.h>
#include <errno.h>
#include <signal.h>
#include <string.h>
#include <sys/wait.h>
#include <sys/stat.h>

#include <gfunction.h>
#include <cc++/thread.h>

#include <cstdlib>  // for exit
#include <iostream>
#include <iomanip>
#include <string>


using namespace std;
using namespace ost;


#define CLI_INPUT   "cli_input_fifo"
#define CLI_OUTPUT  "cli_output_fifo"
#define CLI_CONSOLE "cli_console_fifo"
#define MY_PORT_ID  8000

static int wait_for_cli;
static int input_ctr=0;
static int output_ctr=0;

/** data from GUI to cli via server
*/
class input_redirector : public Thread 
{
 private:
  int  sock_id;
  char pipe_name[255];
  int  pipe_id;
  char char_buf[255];
  int  num_read;
  int  i;

  ofstream debug;

public:
  /// initialize the socket and named pipe with sockid and pipeid
  void Setup(int in_sock_id, int in_pipe_id) {
    sock_id = in_sock_id;
    pipe_id = in_pipe_id;
  }
  /// initialize the socket and named pipe with sockid and pipename
  void Setup(int in_sock_id, const char * new_pipe_name) {
    sock_id = in_sock_id;
    strcpy(pipe_name, new_pipe_name);
    ostringstream sout;
    sout << "debug_input_" << input_ctr++ << ".err";
    debug.open(sout.str().c_str());
    debug << "Pipe name=" << new_pipe_name << endl;
  }

  /// read from the socket and write to the pipe until terminated
  void run() {
    pipe_id = open(pipe_name, O_WRONLY);
    debug << "Run input::pipe_id=" << pipe_id << endl << flush;
    if (pipe_id < 0) {
      perror("problem opening pipe in input_redirector");
      system("/usr/bin/killall cli");
      std::exit(-1);  // std exit to kill app not just thread with thread::exit
    }

    debug << "Run input::pipe_name="    << pipe_name    << endl << flush;
    debug << "Run input::sock_id="      << sock_id      << endl << flush;
    debug << "Run input::wait_for_cli=" << wait_for_cli << endl << flush;

    while (wait_for_cli) {
      debug << " sock_id=" << sock_id << endl << flush;
      debug << " input_redirector:: reading from pipe" << endl << flush;
      num_read = read(sock_id, char_buf, sizeof(char_buf));
      /* */
      debug << "copying from socket to pipe start num_read=" << num_read << endl << flush;
      for (int i = 0; i < num_read; i++)
        debug << char_buf[i];
      debug << "copying from socket to pipe end" << endl << flush;
      /* */
      if (write(pipe_id, char_buf, num_read) < 0) {
        if (true) {
          perror("problem writing data to pipe in input_redirect");
          system("/usr/bin/killall cli");
          std::exit(-1);
        }
      }
    }
  }
  /// close named pipe (why not sockid?)
  void Final() {
    ::close(pipe_id);  // avoid conflict with Thread::close
  }
};

/** data out from cli to GUI via server
*/
class output_redirector : public Thread 
{
public:
  //! Socket id
  int sock_id;

  //! Pipe name
  char pipe_name[255];

  //! Pipe id
  int pipe_id;

  //! Character buffer
  char char_buf[255];

  //! Number read
  int num_read;

  //! Counter
  int i;

  //! Debug output stream
  ofstream debug;

public:
  /// setup the named pipe and socket with pipeid and sockid
  void Setup(int in_pipe_id, int in_sock_id) {
    sock_id = in_sock_id;
    pipe_id = in_pipe_id;
  }

  /// setup the named pipe and socket with pipename and sockid
  void Setup(const char * new_pipe_name, int in_sock_id) {
    sock_id = in_sock_id;
    strcpy(pipe_name, new_pipe_name);
    ostringstream sout;
    sout << "debug_output_" << output_ctr++ << ".err";
    debug.open(sout.str().c_str());
    debug << "Pipe name=" << new_pipe_name << endl;
  }

  /// read from the pipe and write to the socket until terminated
  void run() {
    pipe_id = open(pipe_name, O_RDONLY);
    debug << "output::pipe_id=" << pipe_id << endl << flush;
    if (pipe_id < 0) {
      perror("problem opening pipe in output_redirector");
      system("/usr/bin/killall cli");
      std::exit(-1);
    }

    debug << "output::pipe_name="    << pipe_name    << endl << flush;
    debug << "output::sock_id="      << sock_id      << endl << flush;
    debug << "output::wait_for_cli=" << wait_for_cli << endl << flush;

    while (wait_for_cli) {
      debug << " output_redirector:: reading from pipe" << endl;
      num_read = read(pipe_id, char_buf, sizeof(char_buf));
       
      debug << "copying from pipe to socket start num_read=" << num_read << endl;
      for (int i = 0; i < num_read; i++)
        debug << char_buf[i];
      debug << "copying from pipe to socket start" << endl;
      
      if (send(sock_id, char_buf, num_read, 0) < 0) {
        if (true) {
          perror("problem sending data to socket in output_redirect");
          system("/usr/bin/killall cli");
          std::exit(-1);
        }
      }
    }
  }
  /// close named pipe
  void Final() {
    ::close(pipe_id);  // avoid conflict with Thread::close
  }
};


#endif
