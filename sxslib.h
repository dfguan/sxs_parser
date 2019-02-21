/*
 * =====================================================================================
 *
 *       Filename:  sxslib.h
 *
 *    Description:  library header for sxs file 
 *
 *        Version:  1.0
 *        Created:  21/02/2019 10:39:53
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Dengfeng Guan (D. Guan), dfguan@hit.edu.cn
 *   Organization:  Center for Bioinformatics, Harbin Institute of Technology
 *
 * =====================================================================================
 */

#ifndef SXS_LIB
#define SXS_LIB

#include <stdio.h>
#include <stdint.h>
#include <sys/types.h>

#ifndef KSTRING_T
#define KSTRING_T kstring_t 
typedef struct __kstring_t {
	size_t l, m;
	char *s;
}kstring_t;
#endif

typedef struct {
	int n_fns;
	kstring_t fns[2];
	int n_rds[2];
	int types[2];
}sxs_hdr_t;

typedef struct {
	int aid, bid;
	int as, ae, bs, be;
	kstring_t cigar;
	int match, diff;
	int maq;
	int n_tp, m_tp;
	int *tps;
}sxs_unit_t;

typedef struct {
	void *fp;
	kstring_t buf;
	char last_char;
}sxs_file_t;


#ifdef __cplusplus
extern "C" {
#endif 

sxs_file_t *sxs_open(const char *fn); //open sxs file and return file handle
int sxs_close(sxs_file_t *fp); //close sxs file

int sxs_read_unit(sxs_file_t *fp, sxs_unit_t *su); //read a single sxs record
int sxs_read_blk(sxs_file_t *fp, sxs_unit_t *su, int n); //read bulk sxs records 

sxs_unit_t *sxs_unit_init(int n); //initiate sxs unit data structure
int sxs_unit_destroy(); //release memory for sxs units

int sxs_write(FILE *stream, char *format, ...);
#ifdef __cplusplus
}
#endif


#endif

