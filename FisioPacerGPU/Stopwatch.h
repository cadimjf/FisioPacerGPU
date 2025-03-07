//            QMcBeaver
//
//         Constructed by 
//
//     Michael Todd Feldmann 
//              and 
//   David Randall "Chip" Kent IV
//
// Copyright 2000.  All rights reserved.
//
// drkent@users.sourceforge.net mtfeldmann@users.sourceforge.net

#ifndef Stopwatch_H
#define Stopwatch_H

#ifdef _WIN32

#include <sys/timeb.h>
#include <sys/types.h>
#include <Windows.h>
#include <iostream>

struct timezone {
    int tz_minuteswest;
    int tz_dsttime;
};


#ifndef _WINSOCKAPI_
struct timeval2 {
        long    tv_sec;         /* seconds */
        long    tv_usec;        /* and microseconds */
};
#endif

static void gettimeofday(struct timeval* t,void* timezone){
	struct _timeb timebuffer;
    _ftime( &timebuffer );
	t->tv_sec=timebuffer.time;
    t->tv_usec=1000*timebuffer.millitm;
}
#else
#include <sys/time.h>
#endif


#include <string>
//#include <strstream>
#include <iostream>

#ifdef PARALLEL
#include <mpi.h>
#endif

using namespace std;

/**
  An accurate software stopwatch.
  */

class Stopwatch
{
 private:
  long stime1, stime2, result_ms, total_ms;
  int milli1, milli2;
  bool running;
  struct timeval tp;
  struct timezone tz;
  
 public:
  /**
    Creates an instance of the stopwatch that is zeroed and not running.
    */

  Stopwatch();


  /**
    Resets and stops the stopwatch.
    */

  void reset();


  /**
    Starts the stopwatch.
    */

  void start();


  /**
    Stops the stopwatch.
    */

  void stop();


  /**
    Gets the time in milliseconds.
    */

  long timeMS();


  /**
    Returns true if the stopwatch is running and false otherwise.
    */

  bool isRunning();


  /**
     Sets two objects equal.
    */

  void operator = ( const Stopwatch &rhs);

  /**
    Returns a stopwatch which contains the total time from two
    stopwatch objects.
    */
  Stopwatch operator+(Stopwatch & rhs);


#ifdef PARALLEL

private:
  /**
    A flag which tells if MPI_TYPE has been generated.
    */

  static bool mpiTypeCreated;

  
  /**
    Build MPI_TYPE.
    */

  static void buildMpiType();

  /**
    Build MPI_REDUCE.
    */

  static void buildMpiReduce();


  /**
    An MPI function which allows MPI_Reduce to be used in adding Stopwatches.
    */

  static void Reduce_Function(Stopwatch *in, Stopwatch *inout, 
                              int *len, MPI_Datatype *dptr);

public:

  /**
    The MPI data type for a Stopwatch.
    */
  static MPI_Datatype MPI_TYPE;


  /** 
    The MPI operation for performing MPI_Reduce on Stopwatch objects.
    */
  
  static MPI_Op MPI_REDUCE;

#endif

};

#endif
