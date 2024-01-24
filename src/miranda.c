/* miranda.c */
/*--------------------------------------------------------------------------------*/
/* miRanda- An miRNA target scanner, aims to predict mRNA targets for microRNAs,  */
/* using dynamic-programming alignment and thermodynamics                         */
/* 										  */
/* Copyright (C) (2003) Memorial Sloan-Kettering Cancer Center, New York          */
/*                                                                                */
/* Distributed under the GNU Public License (GPL)                                 */
/* See the files 'COPYING' and 'LICENSE' for details				  */
/*                                                                                */
/* Authors: Anton Enright, Bino John, Chris Sander and Debora Marks               */
/* Email: mirnatargets@cbio.mskcc.org - reaches all authors                       */
/*                                                                                */
/* Written By: Anton Enright (enrighta@mskcc.org)                                 */
/*                                                                                */
/* Please send bug reports to: miranda@cbio.mskcc.org                             */
/*                                                                                */
/* If you use miRanda in your research please cite:                               */
/* Enright AJ, John B, Gaul U, Tuschl T, Sander C and Marks DS;                   */
/* (2003) Genome Biology; 5(1):R1.                                                */
/*									          */
/* This software will be further developed under the open source model,           */
/* coordinated by Anton Enright and Chris Sander:				  */
/* miranda@cbio.mskcc.org (reaches both).					  */
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

#include <config.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "utils.h"
#include "miranda.h"

/* Set Default Parameter Values */

int length_5p_for_weighting = 8;   /* The 5' sequence length to be weighed  except for the last residue*/
int length_3p_for_weighting = 0;
int overlap_cutoff = 0;
int CURR = 0;
double scale = 4.0;     /* The 5' miRNA scaling parameter             */
int nomodel = 0;      /* Strict alignment model on/off              */
double gap_open = -8;   /* Gap-open Penalty                           */
double gap_extend = -2;   /* Gap-extend Penalty                         */
double score_threshold = 50;    /* SW Score Threshold for reporting hits      */
double energy_threshold = -20;  /* Energy Threshold (DG) for reporting hits   */
int verbosity = 1;    /* Verbose mode on/off                        */
int outfile = 0;      /* Dump to file on/off                        */
int truncated = 0;    /* Truncate sequences on/off                  */
int do_shuffle = 0;   /* Generate statistics using seq shuffling    */
int no_energy = 0;    /* Turn off Vienna Energy Calcs - FASTER      */
double average = 0;     /* Some statistics for shuffled searches      */
double stdev = 0;
double z_threshold = 5.0;   /* Z-Score threshold >=           */
int shuffle_window = 10;    /* Size of shuffling window       */
int total_shuffles = 100;   /* Total number of shuffles       */
unsigned int uniform = 0;     /* Uniform Shuffling mode on/off  */
int total_hits = 0;   /* Generic counter for alignments */
FILE *fpout = NULL;

int
main (int argc, char *argv[])
{

  char filename1[200];
  char filename2[200];
  char fileout[200];
  FILE *fp1 = 0;
  FILE *fp2 = 0;
  fpout = stdout; 

  /* Command-line parsing begins here */
  parse_command_line (argc, argv, filename1, filename2, fileout);

  /* Now check our input and output files can be accessed / created */

  if ((fp1 = fopen (filename1, "r")) == NULL)
    {
      fprintf (stderr, "Error: Cannot open file %s\n", filename1);
      exit (1);
    }

  if ((fp2 = fopen (filename2, "r")) == NULL)
    {
      fprintf (stderr, "Error: Cannot open file %s\n", filename2);
      exit (1);
    }


  if ((outfile) && ((fpout = fopen (fileout, "w")) == NULL))
    {
      fprintf (stderr, "Error: Cannot create output file %s\n", fileout);
      exit (1);
    }

  if (verbosity)
    {
      print_parameters (filename1, filename2, fpout);
    }

  /* Everything looks good.... Start the Scan! */

  find_targets (fp1, fp2, fpout, filename2);
  exit (0);
}
