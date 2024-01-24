/* miranda.h									  */
/*--------------------------------------------------------------------------------*/
/* miRanda- An miRNA target scanner, aims to predict mRNA targets for microRNAs,  */
/* using dynamic-programming alignment and thermodynamics                         */
/*                                                                                */
/* Copyright (C) (2003) Memorial Sloan-Kettering Cancer Center, New York          */
/*                                                                                */
/* Distributed under the GNU Public License (GPL)                                 */
/* See the files 'COPYING' and 'LICENSE' for details                              */
/*                                                                                */
/* Authors: Anton Enright, Bino John, Chris Sander and Debora Marks               */
/* Email: mirnatargets@cbio.mskcc.org - reaches all authors                       */
/*                                                                                */
/* Written By: Anton Enright (enrighta@mskcc.org)                                 */
/*                                                                                */
/* Please send bug reports to: miranda@cbio.mskcc.org                             */
/*                                                                                */
/* If you use miRanda in your research please cite:                               */
/* Enright AJ, John B, Gaul U, Tuschl T, Sander C and Marks DS;                   
*/
/* (2003) Genome Biology; 5(1):R1.                                                */
/*                                                                                */
/* This software will be further developed under the open source model,           */
/* coordinated by Anton Enright and Chris Sander:                                 */
/* miranda@cbio.mskcc.org (reaches both).                                         */
/*--------------------------------------------------------------------------------*/
/*
 * Copyright (C) (2003) Memorial Sloan-Kettering Cancer Center
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either
 * version 2 of the License, or (at your option) any later
 * version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
 */
#ifndef __MIRANDA_H
#define __MIRANDA_H

/* Setup Traceback Flags */
#define DIAG 0
#define UP 1
#define LEFT -1
extern int CURR;
/* Setup Strands */
#define FORWARD 0
#define REVERSE 1

/* Adjustable Algorithm Parameters */
extern double scale;
extern int nomodel;
extern double gap_open;
extern double gap_extend;
extern double average;
extern double stdev;
extern double z_threshold;
extern double score_threshold;
extern double energy_threshold;
extern int length_5p_for_weighting;
extern int length_3p_for_weighting;
extern int do_shuffle;
extern int no_energy;
extern int shuffle_window;
extern int total_shuffles;
extern int verbosity;
extern unsigned int uniform;
extern int outfile;
extern int truncated;
extern int total_hits;
extern int seqlen1_global;
extern int overlap_cutoff;
extern FILE *fpout;

/* Generic structure to store individual hit information */
typedef struct hit_struct
{

  double score;
  int query_start;
  int query_end;
  int ref_start;
  int ref_end;
  char *alignment[3];
  char *rest[6];

}
hit_struct;

/* Generic structure to store consecutive UTR hit information */
typedef struct final_score
{
  double scan_score;
  int no_hits;
  double max_hit;
  double max_score;
  double total_score;
  char positional[1024];
}
final_score;

/* Score structure for non-optimal path detection */
typedef struct score_struct
{

  int score;
  int path;
  int i;
  int j;

}
score_struct;

/* Declare Functions */
double max (double, double, double);
int max_finder_fourstates (int, int, int);
int max_finder_threestates (int, int);
int max_finder_and_track_threestates (int, int, int);

int score5p (char, char);
int score (char, char);

void print2dmatrix (int **, int ,int);
void print1dmatrix (int *, int);

void traceback (int **, int ***, char *, char *, int, int, int,
		hit_struct *, int);
int testfor_overlap (int *,int *,int *,int,int);

long readinseq (long, FILE *, char *, char *, char *);
void shuffle (char *, int, int);
void build_matrix (int **, int ***, int **, int **,int **,int **, char *, char *, int, int,
		   score_struct *);

void get_nt_nt_seq_scores (int **,char *, char *, int, int);
double build_matrix_quick (int **, int **, char *, char *, int, int);



void revstring (char s[]);
void irand (int);
int nrand (int);
void printhit (char *, char *, hit_struct *, char *, char *, int, double, double, FILE *);
void print_license();
void nullify_j_overlaps (int,score_struct *);
void print_small_license();
double get_energy (hit_struct *);
double vfold (char *);
void clear_hit (hit_struct *, int, int);
void clear_matrix (double **, int, int, int, int);
void initialize_bases ();
int cmpscores (const void *, const void *);
int find_targets (FILE *, FILE *, FILE *, char *);
double do_alignment (int **, int ***, int **, int **, int **, int **, char *, char *,
		     score_struct *, hit_struct *, int, int, int,
		     final_score *, int, char *, char *, FILE *);
void print_banner();
void print_usage();
void print_options();

/* output.c */
int parse_command_line (int argc, char *argv[], char *filename1, char *filename2,
        char *fileout);

#endif