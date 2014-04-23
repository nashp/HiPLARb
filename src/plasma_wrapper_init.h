#ifndef PLASMA_WRAPPER_INIT
#define PLASMA_WRAPPER_INIT

int plasma_wrapper_init();
int plasma_wrapper_deinit();

extern int R_PLASMA_MAJOR;
extern int R_PLASMA_MINOR;
extern int R_PLASMA_MICRO;

extern int R_PLASMA_NUM_THREADS;
extern int R_PLASMA_SCHED;

#endif
