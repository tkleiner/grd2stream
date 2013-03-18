/*
 * Copyright (c) 2002 by Louis Zechtzer
 *
 * Permission to use, copy and distribute this software is hereby granted
 * under the terms of version 2 or any later version of the GNU General Public
 * License, as published by the Free Software Foundation.
 *
 * THIS SOFTWARE IS PROVIDED IN ITS "AS IS" CONDITION, WITH NO WARRANTY
 * WHATSOEVER. NO LIABILITY OF ANY KIND FOR ANY DAMAGES WHATSOEVER RESULTING
 * FROM THE USE OF THIS SOFTWARE WILL BE ACCEPTED.
 */
/* 
 * Authors: Louis Zechtzer (lou@clarity.net)
 */

#ifndef __LOG_H
#define __LOG_H

#include <stdio.h>

#define LOG_FILE "/var/log/gmttools"
#define LOG_BUF_SIZE 160 /* allows 2-line log message */

/*
 * Define logging severity levels.  These values are based on 
 * those in syslog.h.
 * 
 * The low 4 bits are for finer grained debugging messages.
 */

#define OM_LOG_CRITICAL		0	/* critical conditions */
#define OM_LOG_ALERT		(1<<8)	/* immediate action required */
#define OM_LOG_NOTICE		(2<<8)  /* normal/important */
#define OM_LOG_INFO		(3<<8)  /* information */
#define OM_LOG_DEBUG		(4<<8)  /* debugging information */

#define OM_LOG_PERROR		(1<<7)	/* write strerr(errno) to log also */

#define DEBUG_TRACE_NONE	0x0	/* no debug tracing displayed */
#define DEBUG_TRACE_CALLS	0x1	/* leaving or entering functions */
#define DEBUG_TRACE_SEND        0x2	/* sending messages/datagrams */
#define DEBUG_TRACE_RECV	0x4	/* receiving messages/datagrams */
#define DEBUG_TRACE_PRINT	0x8	/* nonspecifc messages */
#define DEBUG_TRACE_NET		0x10	/* network module messages */


#define DEBUG_TRACE_ALL		0xff	/* constant to display all debugging */ 

/* return error codes */
#define LOG_SUCCESS		0x0
#define LOG_ERR_BAD_FILE 	0x1
#define LOG_ERR_TYPE_UNKNOWN	0x2
#define LOG_FAILURE		0x3

#define OM_LOG_TO_TXT(log_type) log_types_txt[log_type >> 8]

/* following constants passed to log_initialize function */
#define LOG_TO_SYSLOG 		0x0
#define LOG_TO_FILE		0x1
#define LOG_TO_STDERR		0x2

/* 
 * This can be generalized to have multiple logging mechanisms in the future.
 * (e.g. log one set of information to a file, and other information to 
 * syslog.
 */
typedef struct log_ctx {
	FILE *logfile;
	int debug_mask;
	int save_mask;
	int log_dest;
} log_ctx_t;

extern void log_set_debug(int);
extern void log_toggle_debug();
/*extern void log(int type, const char *, ...);*/
extern void log_printf(int type, const char *, ...);
extern int log_initialize(int);
extern int log_finalize();

#ifdef TESTING
extern void log_test();
#endif

extern const char *program_name; /* from main.c */

#endif /* __LOG_H */
