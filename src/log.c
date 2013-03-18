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

#include <stdio.h>
#include <stdarg.h>
#include <syslog.h>
#include <string.h>
#include <errno.h>

#include "log.h"



static int log_initialized;
static log_ctx_t ctx;

static const char *log_types_txt[] = 
	{ "CRITICAL", "ALERT", "NOTICE", "INFO", "DEBUG" };

static int log_config_syslog();
static int log_config_file();
static int log_config_stderr();
static int om_log_to_syslog(int);

/*
 * log_initialize: initialize the logging module.  
 * returns:
 *   LOG_ERR_BAD_FILE: output log file cannot be opened
 *   LOG_SUCCESS: upon success.
 *
 * XXX - add check to ensure that sys_initialize was called before here?
 */
int
log_initialize(int style)
{
	int rc = 0;

	/* debug messages are turned off by default */
	ctx.debug_mask = DEBUG_TRACE_NONE;
	ctx.save_mask = -1;

	switch(style) {
	case LOG_TO_SYSLOG:
		rc = log_config_syslog();
		break;
	case LOG_TO_FILE:
		rc = log_config_file();	
		break;
	case LOG_TO_STDERR:
		rc = log_config_stderr();
		break;
	}

	ctx.log_dest = style;
	
	log_initialized = (rc == LOG_SUCCESS) ? 1 : 0;

	return rc;
}

/*
 * log_finalize: clean up by closing selected logging mechanism.
 */
int
log_finalize()
{
	int rc;

	if (!log_initialized) {
		return LOG_FAILURE;
	}

	rc = LOG_SUCCESS;

	switch (ctx.log_dest) {
	case LOG_TO_SYSLOG:
		closelog();
		break;
	case LOG_TO_FILE:
		if (fclose(ctx.logfile) == EOF) {
			rc = LOG_FAILURE;
		}
		break;
	}

	if (rc == LOG_SUCCESS) {
		log_initialized = 0; 
	}

	return rc;
}

/*
 * log: send log messages to their configured location
 *
 * XXX - prepend a time/info to messages
 */
void
log_printf(int type, const char *fmt, ...)
{
	va_list vl;
	char buf[LOG_BUF_SIZE + 1];

	/* log file will be null if not opened properly */
	if (ctx.log_dest != LOG_TO_SYSLOG && ctx.logfile == NULL) {
		return;
	}

	/* If debug message, check debug level and ignore if necessary */	
	if ((type & OM_LOG_DEBUG) && ((type & ctx.debug_mask) == 0)) {
		return;
	}

	va_start(vl, fmt);
	if (vsnprintf(buf, LOG_BUF_SIZE, fmt, vl) == -1) {
		buf[LOG_BUF_SIZE] = '\0';
	}
	va_end(vl);

	switch(ctx.log_dest) {
	case LOG_TO_SYSLOG:
		{
		if (type & OM_LOG_PERROR) {
			syslog(om_log_to_syslog(type), "%s: %s\n", 
				buf, strerror(errno));
		} else {
			syslog(om_log_to_syslog(type), "%s\n", buf);
		}
		break;
		}	
	case LOG_TO_FILE:
	case LOG_TO_STDERR:
		{
		if (type & OM_LOG_PERROR) {
			fprintf(ctx.logfile, "%s: %s: %s\n", 
				OM_LOG_TO_TXT(type), buf, strerror(errno));
		} else {
			fprintf(ctx.logfile, "%s: %s\n", OM_LOG_TO_TXT(type), 
				buf);
		}
		fflush(ctx.logfile);
		break;
		}
	}
}


/* 
 * log_set_debug: set mask that controls which types of debug messages 
 * are logged.
 */
void 
log_set_debug(int level)
{
	ctx.debug_mask = level;
}

/*
 * log_toggle_debug: toggle between full debugging and pre-set level. 
 * can be attached to a SIGUSR1 so logging can be changed at run-time.
 */
void
log_toggle_debug()
{
	if (ctx.save_mask == -1) {
		ctx.save_mask = ctx.debug_mask;
		log_set_debug(DEBUG_TRACE_ALL);
	} else {
		log_set_debug(ctx.save_mask);
		ctx.save_mask = -1;
	}
}

static int
log_config_syslog()
{
	openlog(program_name, LOG_CONS|LOG_PID, LOG_DAEMON);
	return LOG_SUCCESS;
}

static int
log_config_file()
{
	int rc;

	/* XXX - rotate out old file? */
	ctx.logfile = fopen(LOG_FILE, "w");
	if (ctx.logfile == NULL) {
		rc = LOG_ERR_BAD_FILE;
	} else {
		rc = LOG_SUCCESS;
	}

	return rc;
}

static int
log_config_stderr()
{
	ctx.logfile = stderr;
	return LOG_SUCCESS;
}

/* 
 * om_log_to_syslog: convert this module's log message type abstraction to
 * syslog's.
 */
static int
om_log_to_syslog(int logtype)
{
	switch (logtype) {
	case OM_LOG_CRITICAL:
		return LOG_CRIT;
	case OM_LOG_ALERT:
		return LOG_ALERT;
	case OM_LOG_NOTICE:
		return LOG_NOTICE;
	case OM_LOG_INFO:
		return LOG_INFO;
	case OM_LOG_DEBUG:
		return LOG_DEBUG;
	default:
		return LOG_ERR_TYPE_UNKNOWN;
	}
}

#ifdef TESTING
#include <assert.h>

void 
log_test()
{
	int i = 12345;
	
	log_initialize(LOG_TO_SYSLOG);
	log(OM_LOG_INFO, "Logging number %d to syslog", i);
	log_finalize();

	log_initialize(LOG_TO_FILE);
	log(OM_LOG_INFO, "Logging number %d to %s", i, LOG_FILE);
	log_finalize();

	log_initialize(LOG_TO_STDERR);
	log_set_debug(DEBUG_TRACE_RECV);
	log(OM_LOG_DEBUG|DEBUG_TRACE_RECV, "Displaying DEBUG_TRACE_RECV msg");
	log_set_debug(DEBUG_TRACE_NONE);
	log(OM_LOG_DEBUG|DEBUG_TRACE_RECV, "If displayed, there is a bug");

	printf("Expect to see an error message:\n");
	errno = EBADF;
	log(OM_LOG_ALERT|OM_LOG_PERROR, "Displaying EBADF perror");
	log_finalize();
}
#endif
