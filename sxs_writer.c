/*
 * =====================================================================================
 *
 *       Filename:  sxs_writer.c
 *
 *    Description:  write sxs records
 *
 *        Version:  1.0
 *        Created:  21/02/2019 07:54:35
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Dengfeng Guan (D. Guan), dfguan@hit.edu.cn
 *   Organization:  Center for Bioinformatics, Harbin Institute of Technology
 *
 * =====================================================================================
 */

#include <stdio.h>
#include <stdarg.h>
#include "color.h"
#include "sxslib.h"

int sxs_write(FILE *stream, char *format, ...)
{
	va_list args;
	va_start(args, format);
	char fmt = format[0];
	if (fmt == 'A') {
		fprintf(stream, COLOR_RED); //print two integer
		fprintf(stream, "A\t");
		vfprintf(stream, "%u\t%u\n", args);	
	} else if (fmt == 'I') {
		fprintf(stream, COLOR_GREEN);
		fprintf(stream, "I\t"); //print two integer
		vfprintf(stream, "%u\t%u\t%u\t%u\n", args);	
	
	} else if (fmt == 'D') {
		fprintf(stream, COLOR_MAGENTA);
		fprintf(stream, "D\t"); //print two integer
		vfprintf(stream, "%u\n", args);	
	} else if (fmt == 'C') {
		fprintf(stream, COLOR_YELLOW);
		fprintf(stream, "C\t"); //print two integer
		vfprintf(stream, "%s\n", args);	
	} else if (fmt == 'M') {
		fprintf(stream, COLOR_BLUE);
		fprintf(stream, "M\t"); //print two integer
		vfprintf(stream, "%u\n", args);	
	} else if (fmt == 'Q') {
		fprintf(stream, COLOR_CYAN);
		fprintf(stream, "Q\t"); //print two integer
		vfprintf(stream, "%u\n", args);	
	} else if (fmt == 'T' || fmt == 'U' || fmt == 'V') {
		fprintf(stream, "T");
		int nt = va_arg(args, int);
		int *ntps = va_arg(args, int *);
		int i;
		for ( i = 0; i < nt; ++i) 
			fprintf(stream, "\t%d", ntps[i]);	
		fprintf(stream, "\n");	
	} else {
		va_end(args);
		fprintf(stream, COLOR_RESET);
		return -1;
	}	
	fprintf(stream, COLOR_RESET);
	va_end(args);
	return 0;
}




