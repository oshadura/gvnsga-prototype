#ifndef __INSTRUMENTATION__
#define __INSTRUMENTATION__

#include "stdlib.h"
#include "stdio.h"
#include "string.h"
#include "sys/times.h"
#include "sys/vtimes.h"

class CPUManager {

private:
  clock_t lastCPU, lastSysCPU, lastUserCPU;
  int numProcessors;
  unsigned long long lastTotalUser, lastTotalUserLow, lastTotalSys,
      lastTotalIdle;

public:
  CPUManager() {}
  ~CPUManager() {}

  void InitCPU() {
    FILE *file;
    struct tms timeSample;
    char line[128];

    lastCPU = times(&timeSample);
    lastSysCPU = timeSample.tms_stime;
    lastUserCPU = timeSample.tms_utime;

    file = fopen("/proc/cpuinfo", "r");
    numProcessors = 0;
    while (fgets(line, 128, file) != NULL) {
      if (strncmp(line, "processor", 9) == 0)
        numProcessors++;
    }
    fclose(file);
  }

  void InitAllCPU() {
    FILE *file = fopen("/proc/stat", "r");
    fscanf(file, "cpu %llu %llu %llu %llu", &lastTotalUser, &lastTotalUserLow,
           &lastTotalSys, &lastTotalIdle);
    fclose(file);
  }

  double GetCurrentValueAllCPU() {
    double percent;
    FILE *file;
    unsigned long long totalUser, totalUserLow, totalSys, totalIdle, total;

    file = fopen("/proc/stat", "r");
    fscanf(file, "cpu %llu %llu %llu %llu", &totalUser, &totalUserLow,
           &totalSys, &totalIdle);
    fclose(file);

    if (totalUser < lastTotalUser || totalUserLow < lastTotalUserLow ||
        totalSys < lastTotalSys || totalIdle < lastTotalIdle) {
      // Overflow detection. Just skip this value.
      percent = -1.0;
    } else {
      total = (totalUser - lastTotalUser) + (totalUserLow - lastTotalUserLow) +
              (totalSys - lastTotalSys);
      percent = total;
      total += (totalIdle - lastTotalIdle);
      percent /= total;
      percent *= 100;
    }

    lastTotalUser = totalUser;
    lastTotalUserLow = totalUserLow;
    lastTotalSys = totalSys;
    lastTotalIdle = totalIdle;

    return percent;
  }

  double GetCurrentValueCPU() {
    struct tms timeSample;
    clock_t now;
    double percent;
    now = times(&timeSample);
    if (now <= lastCPU || timeSample.tms_stime < lastSysCPU ||
        timeSample.tms_utime < lastUserCPU) {
      percent = -1.0;
    } else {
      percent = (timeSample.tms_stime - lastSysCPU) +
                (timeSample.tms_utime - lastUserCPU);
      percent /= (now - lastCPU);
      percent /= numProcessors;
      percent *= 100;
    }
    lastCPU = now;
    lastSysCPU = timeSample.tms_stime;
    lastUserCPU = timeSample.tms_utime;

    return percent;
  }
};

#endif