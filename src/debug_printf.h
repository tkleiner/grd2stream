#ifndef DEBUG_PRINTF_H
#define DEBUG_PRINTF_H

/**
 *
 */
enum {
  DEBUG_INFO,
  DEBUG_WARNING,
  DEBUG_ERROR
};

/**
 *
 */
void debug_printf(int dp, char *format, ...);

#endif
