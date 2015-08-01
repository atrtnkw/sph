#include <sys/stat.h>
#include <sys/time.h>

PS::F64 getWallclockTime() {
    struct timeval tv;
    gettimeofday(& tv, NULL);
    return ((double)(tv.tv_sec) + (double)(tv.tv_usec) * 1e-6);
}
