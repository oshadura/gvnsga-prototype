#ifdef ENABLE_PERFMON

#include "PFMWatch.h"
#include <iostream>

PFMWatch::PFMWatch() {
  EVENTSTRING.resize(GROUPSIZE);

  //    EVENTSTRING[0]="cs";
  //    EVENTSTRING[1]="migrations";
  EVENTSTRING[2] = "cycles";
  EVENTSTRING[3] = "instructions";
  EVENTSTRING[4] = "branch-misses";
  EVENTSTRING[5] = "L1-dcache-load-misses";
  EVENTSTRING[6] = "L1-icache-load-misses";
  EVENTSTRING[7] = "branches";
  // EVENTSTRING[2]="cache-misses";

  // EV/ENTSTRING[4]="cycles";
  // EVENTSTRING[5]="instructions";
  EVENTSTRING[0] = "idle-cycles-frontend";
  EVENTSTRING[1] = "cs"; // idle-cycles-backend";

  ret = pfm_initialize();
  if (ret != PFM_SUCCESS)
    errx(1, "cannot initialize library: %s", pfm_strerror(ret));

  for (unsigned int i = 0; i < GROUPSIZE; ++i) {
    memset(&attr[i], 0, sizeof(perf_event_attr)); // initialize attr;
    count[i] = 0;

    /*
     * 1st argument: event string
     * 2nd argument: default privilege level (used if not specified in the
     * event string)
     * 3rd argument: the perf_event_attr to initialize
     */
    ret = pfm_get_perf_event_encoding(
        EVENTSTRING[i].c_str(), PFM_PLM0 | PFM_PLM3, &attr[i], NULL, NULL);
    if (ret != PFM_SUCCESS)
      errx(1, "cannot find encoding: %s", pfm_strerror(ret));

    /*
     * request timing information because event may be multiplexed
     * and thus it may not count all the time. The scaling information
     * will be used to scale the raw count as if the event had run all
     * along
     */
    attr[i].read_format =
        PERF_FORMAT_TOTAL_TIME_ENABLED | PERF_FORMAT_TOTAL_TIME_RUNNING;
    attr[i].disabled = 1;

    /*
     * create the event and attach to self
     * Note that it attaches only to the main thread, there is no inheritance
     * to threads that may be created subsequently.
     *
     * if mulithreaded, then getpid() must be replaced by gettid()
     */
    if (i == 0) {
      fd[i] = perf_event_open(&attr[i], 0, -1, -1, 0);
      if (fd[i] < 0)
        err(1, "cannot create main event");
    } else {
      // create other events in same event group
      fd[i] = perf_event_open(&attr[i], 0, -1, fd[0],
                              0); // this is now an event group with fd
      if (fd[i] < 0)
        err(1, "cannot create child event");
    }
  }

  for (unsigned int i = 0; i < GROUPSIZE; ++i) {
    // start the counter for each event
    ret = ioctl(fd[i], PERF_EVENT_IOC_ENABLE, 0);
    if (ret)
      err(1, "ioctl(enable) failed");
  }
}

PFMWatch::~PFMWatch() {
  /*
   * stop counting
   */
  for (unsigned int i = 0; i < GROUPSIZE; ++i) {
    ret = ioctl(fd[i], PERF_EVENT_IOC_DISABLE, 0);
    if (ret)
      err(1, "ioctl(disable) failed");

    close(fd[i]);
  }

  /* free libpfm resources cleanly */
  pfm_terminate();
}

void PFMWatch::Start() {
  /*
   * start counting now
   */
  for (unsigned int i = 0; i < GROUPSIZE; ++i) {
    ret = read(fd[i], startvalues + i * 3, sizeof(uint64_t) * 3);
    if (ret != sizeof(uint64_t) * 3)
      err(1, "cannot read results: %s", strerror(errno));
  }
}

void PFMWatch::Stop() {
  /*
   * read out and calculate counters
   */
  for (unsigned int i = 0; i < GROUPSIZE; ++i) {
    ret = read(fd[i], stopvalues + i * 3, sizeof(uint64_t) * 3);
    if (ret != sizeof(uint64_t) * 3)
      err(1, "cannot read results: %s", strerror(errno));

    if (stopvalues[3 * i + 2])
      count[i] =
          (uint64_t)((double)(stopvalues[3 * i + 0] - startvalues[3 * i + 0]) *
                     (stopvalues[3 * i + 1] - startvalues[3 * i + 1]) /
                     (1. * (stopvalues[3 * i + 2] - startvalues[3 * i + 2])));
  }
}

void PFMWatch::HeatUp() {
  for (unsigned int i = 0; i < 5; ++i) {
    Start();
    Stop();
  }
}

double PFMWatch::getDeltaSecs() { return (count[1]); }
//  double getDeltaSecs() { return
//  (count[0]-countoverhead[0])/(count[1]-countoverhead[1]); }

// find out about the background overhead
double PFMWatch::getOverhead(int N) {
  HeatUp();
  unsigned long long Taccum = 0L;
  for (unsigned int i = 0; i < GROUPSIZE; ++i) {
    countoverhead[i] = 0;
  }

  for (unsigned int i = 0; i < sizeof(N); i++) {
    Start();
    Stop();
    for (unsigned int i = 0; i < GROUPSIZE; ++i) {
      fprintf(stderr, "#count of event %d is %ld\n", i, count[i]);
      countoverhead[i] += count[i];
    }
  }
  for (unsigned int i = 0; i < GROUPSIZE; ++i) {
    countoverhead[i] /= (1. * N);
    fprintf(stderr, "#OVERHEAD of %s = %lf\n", EVENTSTRING[i].c_str(),
            countoverhead[i]);
  }

  return countoverhead[1];
  //    return Taccum/(1.*N);
}

void PFMWatch::printSummary() {
  std::cout << "# PFM Summary (for measured code section) " << std::endl;
  std::cout << "# --------------------------------------- " << std::endl;
  for (unsigned int i = 0; i < GROUPSIZE; i++) {
    std::cout << std::left << "# " << std::setw(25) << EVENTSTRING[i]
              << "\t: " << count[i] << std::endl;
  }
}
#endif
