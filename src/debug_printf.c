
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <stdio.h>
#include <string.h>
#include <stdarg.h>

#include "debug_printf.h"


void debug_printf(int dp,char *format, ...)
{
#ifdef DEBUG
  va_list arglist;
  va_start(arglist,format);
  if(dp == DEBUG_INFO)
    {
      printf("#INFO:\t");
    }
  else if(dp == DEBUG_WARNING)
    {
      printf("#WARNING:\t");
    }
  else
    {
      printf("#ERROR:\t");
    }
  vprintf(format, arglist);
  if(format[strlen(format)-1] != '\n')
    {
      printf("\n");
    }
  va_end(arglist);
#endif
}
