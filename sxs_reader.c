/*
 * =====================================================================================
 *
 *       Filename:  sxs_parser.c
 *
 *    Description:  parser for sxs 
 *
 *        Version:  1.0
 *        Created:  20/02/19 11:30:57
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Dengfeng Guan (D. Guan), dfguan@hit.edu.cn
 *   Organization:  Center for Bioinformatics, Harbin Institute of Technology
 *
 * =====================================================================================
 */

#include <zlib.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>


#include "sxslib.h"
#include "kseq.h"

KSTREAM_INIT(gzFile, gzread, gzseek, 0x10000)

sxs_file_t *sxs_open(const char *fn)
{
	kstream_t *ks;
	gzFile fp;
	fp = fn && strcmp(fn, "-") ? gzopen(fn, "r") : gzdopen(fileno(stdin), "r");
	if (fp == 0) return 0;
	ks = ks_init(fp);
	sxs_file_t *sf = (sxs_file_t *)calloc(1, sizeof(sxs_file_t));
	sf->fp = ks;	
	return sf;
}
int sxs_close(sxs_file_t *fp) {
	if (!fp) return 0;
	if (fp->buf.s) free(fp->buf.s);
	kstream_t *ks = (kstream_t *)fp->fp;
	gzclose(ks->f);
	ks_destroy(ks);
	free(fp);
	return 0;
}

sxs_unit_t *sxs_unit_init(int n)
{
	sxs_unit_t *sus = (sxs_unit_t *)calloc(n, sizeof(sxs_unit_t));	
	return sus;
}

int sxs_unit_destroy(sxs_unit_t *sus,int n) 
{
	int i;
	if (sus) {
		for ( i = 0; i < n; ++i) {
			if (sus[i].cigar.s) free(sus[i].cigar.s);
			if (sus[i].tps) free(sus[i].tps);	
		}
		free(sus);	
	} 
	return 0;
}

int parse_a(char *s, int l, sxs_unit_t *su)
{
	char *q;
	int i, t;
	/*fprintf(stderr, "%s\n", s);*/
	i = 0; while (i < l && s[i] == '\t') ++i; 
	for (t = 0, q = s + i; i <= l; ++i) {
		if (i < l && s[i] != '\t') continue;
		s[i] = 0;
		if (t == 0) su->aid = atoi(q);
		else if (t == 1) su->bid = atoi(q);
		++t, q = i < l? &s[i+1] : 0;
	}
	if (t < 2) return -1;
	return 0;
}

int parse_m(char *s, int l, sxs_unit_t *su)
{
	char *q;
	int i, t;
	/*fprintf(stderr, "%s\n", s);*/
	i = 0; while (i < l && s[i] == '\t') ++i; //skip all tabs
	for (t = 0, q = s + i; i <= l; ++i) {
		if (i < l && s[i] != '\t') continue;
		s[i] = 0;
		if (t == 0) su->match = atoi(q);
		++t, q = i < l? &s[i+1] : 0;
	}
	if (t < 1) return -1;
	return 0;
}

int parse_i(char *s, int l, sxs_unit_t *su)
{
	char *q;
	int i, t;
	/*fprintf(stderr, "%s\n", s);*/
	i = 0; while (i < l && s[i] == '\t') ++i; 
	for (t = 0, q = s + i; i <= l; ++i) {
		if (i < l && s[i] != '\t') continue;
		s[i] = 0;
		if (t == 0) su->as = atoi(q);
		else if (t == 1) su->ae = atoi(q);
		else if (t == 2) su->bs = atoi(q);
		else if (t == 3) su->be = atoi(q);
		++t, q = i < l? &s[i+1] : 0;
	}
	if (t < 4) return -1;
	return 0;
}


int parse_t(char *s, int l, sxs_unit_t *su)
{
	char *q;
	int i, t;
	su->n_tp = 0;
	i = 0; while (i < l && s[i] == '\t') ++i; 
	for (t = 0, q = s + i; i <= l; ++i) {
		if (i < l && s[i] != '\t') continue;
		s[i] = 0;
		if (t == 0) su->n_tp = atoi(q);
		++t, q = i < l? &s[i+1] : 0;
	}

	if (t < 1 || !su->n_tp) return -1; //no tracing points
	//allocate memory space if necessary 
	if (su->n_tp >= su->m_tp) {
		su->m_tp = su->n_tp + 16;
		su->tps = (int *)(sizeof(int) * su->m_tp);
	}	
	for ( t = 0, q = s + i; i <= l; ++i) {
		if (i < l && s[i] != '\t') continue;
		s[i] = 0;
		if (t < su->n_tp) // more than expected 
			su->tps[t++] = atoi(q);
		else {
			++t;
			break;	
		}
		q = i < l? &s[i+1] : 0;
	}	
	if (t != su->n_tp) return -1;
	return 0;

}

int parse_q(char *s, int l, sxs_unit_t *su)
{
	char *q;
	int i, t;
	/*fprintf(stderr, "%s\n", s);*/
	i = 0; while (i < l && s[i] == '\t') ++i; 
	for (t = 0, q = s + i; i <= l; ++i) {
		if (i < l && s[i] != '\t') continue;
		s[i] = 0;
		if (t == 0) su->maq = atoi(q);
		++t, q = i < l? &s[i+1] : 0;
	}
	if (t < 1) return -1;
	return 0;
}	

int parse_c(char *s, int l, sxs_unit_t *su)
{
	char *q;
	int i, t;
	su->cigar.l = 0;	
	i = 0; while (i < l && s[i] == '\t') ++i; 
	for (t = 0, q = s + i; i <= l; ++i) {
		if (i < l && s[i] != '\t') continue;
		s[i] = 0;
		if (t == 0) {
			size_t cigar_len = &s[i] - q;
			if (cigar_len >= su->cigar.m) {
				su->cigar.m = cigar_len + 16;
				su->cigar.s = (char *)realloc(su->cigar.s, sizeof(char) * su->cigar.m);			
			}
			strncpy(su->cigar.s, q, cigar_len+1);
			su->cigar.l = cigar_len;
		}
		++t, q = i < l? &s[i+1] : 0;
	}
	if (t < 1) return -1;
	return 0;
}	

int parse_hdr(char *s, int l, sxs_hdr_t *sh)
{
	char *q;
	int i, t;
	if (sh->n_fns >= 2) return -1;
	/*fprintf(stderr, "%s\n", s);*/
	i = 0; while (i < l && s[i] == '\t') ++i; 
	kstring_t *fn = &sh->fns[sh->n_fns];
	int *nrds = &sh->n_rds[sh->n_fns];
	int *type = &sh->types[sh->n_fns];
	for (t = 0, q = s + i; i <= l; ++i) {
		if (i < l && s[i] != '\t') continue;
		s[i] = 0;
		if (t == 0) {
			size_t name_len = &s[i] - q;
			if (name_len >= sh->fns[sh->n_fns].m) {
				fn->m = name_len + 16;
				fn->s = (char *)realloc(fn->s, sizeof(char) * fn->m);			
			}
			strncpy(fn->s, q, name_len+1);
			fn->l = name_len;
		}
		else if (t == 1)  *type = q[0];
		else if (t == 2)  *nrds = atoi(q);
		++t, q = i < l? &s[i+1] : 0;
	}
	if (t < 2) return -1;
	else 
		++sh->n_fns;
	return 0;


}

int sxs_read_hdr(sxs_file_t *fp, sxs_hdr_t *sh)
{
	//should be called at the first time 
	int c;
	kstream_t *ks = (kstream_t *)fp->fp;
	int ret = 0;	
	while ((c = ks_getc(ks)) != -1) {
		if (c == '<') {
			if (ks_getuntil(ks, KS_SEP_LINE, &fp->buf, 0) < 0) {
				ret = -1; break;	
			} else {
				parse_hdr(fp->buf.s, fp->buf.l, sh);	
			} 
		} else if (c == 'A') {
			fp->last_char = c;
			break;
		} else 
			break;
	} 
	if (c == 'A' && sh->n_fns == 2 && ret == 0) return 0;
	else 
		return -1;	
}

int sxs_read_unit(sxs_file_t *fp, sxs_unit_t *su)
{	
	int req = 0, ret, dret; 	
	kstream_t *ks = (kstream_t *)fp->fp;
	int c = fp->last_char;
	if (fp->last_char == 0) {
		while ((c = ks_getc(ks)) != -1 && c != 'A'); // skip all contents before a line not start with A;  
		if (c == -1) return -1; // the end of the file
		fp->last_char = c;	
	}
	su->cigar.l = su->n_tp = 0;	
	while (1) {
		ret = ks_getuntil(ks, KS_SEP_LINE, &fp->buf, &dret); //read the A line 
		if (ret < 0) return ret; // reach file end and A line missing  
		//don't know which one is required yet	
		switch (c) {
			case 'A':
				ret = parse_a(fp->buf.s, fp->buf.l, su); 
				++req;
				break;
			case 'I':
				ret = parse_i(fp->buf.s, fp->buf.l, su);
				++req;
				break;
			case 'M':
				ret = parse_m(fp->buf.s, fp->buf.l, su);
				break;
			case 'C':
				ret = parse_c(fp->buf.s, fp->buf.l, su);
				break;
			case 'T':
				ret = parse_t(fp->buf.s, fp->buf.l, su);
				break;
			case 'Q':
				ret = parse_q(fp->buf.s, fp->buf.l, su);
				++req;
				break;
			default:
				ret = -1; //undefined line
		}
		if (ret < 0) // syntax error line 
			break;
		c = ks_getc(ks);
		if (c == -1) { //check we got all required? 
			//this is the end of file 
			fp->last_char = 0;
			break;
		} else if (c == 'A') {
			fp->last_char = c;
			break;	
		}				
	}	
	if (ret >= 0 && req == 3) return ret;	
	else 
		return -1;
}

int sxs_read_blk(sxs_file_t *fp, sxs_unit_t *su, int n)
{
	int i;
	for ( i = 0; i < n; ++i) {
		if (sxs_read_unit(fp, &su[i]) < 0) 
			break;
	}
	return i;	
}

int sxs_print_hdr(sxs_hdr_t *s)
{
	if (s) {
		int i;
		for ( i = 0; i < s->n_fns; ++i) 
			fprintf(stdout, "< %s", s->fns[i].s);
	}
	return 0;
}

int sxs_print_unit(sxs_unit_t *s, int n)
{
	int i;
	for ( i = 0; i < n; ++i ) {
		sxs_write(stdout, "A", s[i].aid, s[i].bid);
		sxs_write(stdout, "I", s[i].as, s[i].ae, s[i].bs, s[i].be);
		if (s[i].cigar.l) sxs_write(stdout, "C", s[i].cigar.s);
		if (s[i].n_tp) sxs_write(stdout, "T", s[i].n_tp, s[i].tps);
		sxs_write(stdout, "M", s[i].match);
		sxs_write(stdout, "Q", s[i].maq);
	}	
	return 0;
}

#ifdef MAIN

int main(int argc, char *argv[])
{
	sxs_file_t *sf = sxs_open(argv[1]);
	sxs_unit_t *su = sxs_unit_init(1);
	
	while (sxs_read_unit(sf, su) >= 0) 
		sxs_print_unit(su, 1);
	
	sxs_unit_destroy(su, 1);
	sxs_close(sf);
}
#endif
