////////////////////////////////////////
//
// Sandro Wenzel - sandro.wenzel@cern.ch
//
////////////////////////////////////////
#ifdef ENABLE_PERFMON
#ifndef PFMWATCH_H
#define PFMWATCH_H

// A "timer" measuring hardware performance counter between start and stop
// this is based directly on the libpfm4.4 (perf) library
// this code is adapted from the self-basic.c example of libpfm

#include <cstdlib>

#include <sys/types.h>
#include <inttypes.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <errno.h>
#include <unistd.h>
#include <string.h>
#include <locale.h>
#include <sys/ioctl.h>
#include <err.h>
#include <vector>
#include <perfmon/pfmlib_perf_event.h>
#include <iomanip>

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
  uint64_t GetNInstructions() const { return this->getCounter(3); }
  uint64_t GetBranchMisses() const { return this->getCounter(4); }

private:
  struct perf_event_attr attr[GROUPSIZE];
  int fd[GROUPSIZE];
  int ret;
  std::vector<std::string> EVENTSTRING;
  uint64_t count[GROUPSIZE];
  uint64_t startvalues[3 * GROUPSIZE];
  uint64_t stopvalues[3 * GROUPSIZE];
  double countoverhead[GROUPSIZE];

  //.ClassDef(PFMWatch,1)
};

#endif
#endif
