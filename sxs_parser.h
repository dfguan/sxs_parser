/*
 * =====================================================================================
 *
 *       Filename:  sxs_parser.h
 *
 *    Description:	hdr parser for sxs 
 *
 *        Version:  1.0
 *        Created:  20/02/19 11:13:26
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Dengfeng Guan (D. Guan), dfguan@hit.edu.cn
 *   Organization:  Center for Bioinformatics, Harbin Institute of Technology
 *
 * =====================================================================================
 */

#ifndef SXS_PARSER
#define SXS_PARSER

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
	int n_fns, m_fns;
	kstring_t *fns;
}sxs_hdr_t;

typedef struct {
	int aid, bid;
	int as, ae;	
	int bs, be;
	kstring_t cigar;
	int match;
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

sxs_file_t *sxs_open(const char *fn);
int sxs_close(sxs_file_t *fp);

int sxs_read_unit(sxs_file_t *fp, sxs_unit_t *su);
int sxs_read_blk(sxs_file_t *fp, sxs_unit_t *su, int n);

#ifdef __cplusplus
}
#endif


#endif
