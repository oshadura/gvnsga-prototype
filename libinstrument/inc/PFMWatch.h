////////////////////////////////////////
//
// Sandro Wenzel - sandro.wenzel@cern.ch
//
////////////////////////////////////////
//#ifdef ENABLE_PERF

#ifndef PFMWATCH_H
#define PFMWATCH_H

// A "timer" measuring hardware performance counter between start and stop
// this is based directly on the libpfm4.4 (perf) library
// this code is adapted from the self-basic.c example of libpfm

#include <cstdlib>

#include <err.h>
#include <errno.h>
#include <inttypes.h>
#include <iomanip>
#include <locale.h>
#include <perfmon/pfmlib_perf_event.h>
#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/ioctl.h>
#include <sys/types.h>
#include <unistd.h>
#include <vector>

#define GROUPSIZE 8 // 8 is the maximum number of counted events

class PFMWatch {

public:
  PFMWatch();
  ~PFMWatch();
  void Start();
  void Stop();
  void HeatUp();
  double getDeltaSecs();
  double getOverhead(int N);
  void printSummary();
  uint64_t getCounter(unsigned int event) { return count[event]; }
  std::string getEventName(unsigned int event) { return EVENTSTRING[event]; }
  uint64_t getNumberOfEvents() { return GROUPSIZE; }
  uint64_t GetNICS() const { return this->getCounter(0)/1000000; }
  uint64_t GetNCS() const { return this->getCounter(1); }
  uint64_t GetNC() const { return this->getCounter(2)/1000000; }
  uint64_t GetNI() const { return this->getCounter(3)/1000000; }
  uint64_t GetNBM() const { return this->getCounter(4)/1000000; }
  uint64_t GetNDC() const { return this->getCounter(5)/1000000; }
  uint64_t GetNIC() const { return this->getCounter(6)/1000000; }
  uint64_t GetNB() const { return this->getCounter(7)/1000000; }

private:
  struct perf_event_attr attr[GROUPSIZE];
  int fd[GROUPSIZE];
  int ret;
  std::vector<std::string> EVENTSTRING;
  uint64_t count[GROUPSIZE];
  uint64_t startvalues[3 * GROUPSIZE];
  uint64_t stopvalues[3 * GROUPSIZE];
  double countoverhead[GROUPSIZE];
};

#endif

//#endif
