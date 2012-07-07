/* suftest.c
   Copyright N. Jesper Larsson 1999.
   Program to test suffix sorting function. Reads a sequence of bytes from a
   file and calls suffixsort. This is the program used in the experiments in
   "Faster Suffix Sorting" by N. Jesper Larsson (jesper@cs.lth.se) and Kunihiko
   Sadakane (sada@is.s.u-tokyo.ac.jp) to time the suffixsort function in the
   file qsufsort.c.

   This software may be used freely for any purpose. However, when distributed,
   the original source must be clearly stated, and, when the source code is
   distributed, the copyright notice must be retained and any alterations in
   the code must be clearly marked. No warranty is given regarding the quality
   of this software.

   CHANGES:
   1999-06-14: Fixed preprocessor conditions so that it's possible to have
   CHECK==0 and PRINT==2 simultaneously. Added to comment about PRINT==2.*/

#include <stdlib.h>
#include <stdio.h>
#include <limits.h>

#ifndef CHECK
/* Nonzero CHECK means check that sort is correct. Very slow for repetitive
   files.*/
#define CHECK 0
#endif

#ifndef PRINT
/* 0 for no printing. 1 to print suffix numbers in sorted order. 2 to print
   first MAXPRINT characters of each suffix in sorted order (makes sense only
   if the input is text and COMPACT is 0).*/
#define PRINT 0
#endif
#ifndef MAXPRINT
#define MAXPRINT 10
#endif

#ifndef COMPACT
/* 0 to use alphabet 0...UCHAR_MAX without checking what appears. 1 to limit
   the alphabet to the range l...k-1 that actually appears in the input. 2 to
   transform the alphabet into 1...k-1 with no gaps and minimum k, preserving
   the order.*/
#define COMPACT 2
#endif

#define MIN(a, b) ((a)<=(b) ? (a) : (b))

void suffixsort(int *x, int *p, int n, int k, int l);

int scmp3(unsigned char *p, unsigned char *q, int *l, int maxl)
{
   int i;
   i = 0;
   while (maxl>0 && *p==*q) {
      p++; q++; i++;
      maxl--;
   }
   *l = i;
   if (maxl>0) return *p-*q;
   return q-p;
}

int main(int argc, char *argv[])
{
   int current_char, i, j, *original_text, *suffix_array, *pi;
   int file_size, max_char, min_char;
#if COMPACT==2
   unsigned char compressed_alphabet[UCHAR_MAX+1];
#endif
#if CHECK || PRINT==2
   unsigned char *s;
#endif
   char *filename;
   FILE *file;

   if (argc != 2) {
      fprintf(stderr, "syntax: suftest file\n");
      return 1;
   }

   filename = argv[1];
   if (!(file = fopen(filename, "rb"))) {
      perror(filename);
      return 1;
   }

   if (fseek(file, 0L, SEEK_END)) {
      perror(filename);
      return 1;
   }

   file_size = ftell(file);
   if (file_size == 0) {
      fprintf(stderr, "%s: file empty\n", filename);
      return 0;
   }

   suffix_array = malloc((file_size + 1) * sizeof *suffix_array);
   original_text = malloc((file_size + 1) * sizeof *original_text);
   if (!suffix_array || !original_text) {
      fprintf(stderr, "malloc failed\n");
      return 1;
   }

// TODO: Remove if s that are not used

#if COMPACT==1
   min_char = UCHAR_MAX;
   max_char = 1;
   // TODO: optimize
   for (rewind(file), pi = original_text; pi < original_text + file_size; ++pi) {
      *pi = current_char = getc(file);
      if (current_char < min_char)
         min_char = current_char;
      if (current_char >= max_char)
         max_char = current_char + 1;
   }
#else
   for (rewind(file), pi = original_text; pi < original_text + file_size; ++pi)
      *pi = getc(file);
#if COMPACT==0
   min_char = 0;
   max_char = UCHAR_MAX + 1;
#elif COMPACT==2
   for (i = 0; i <= UCHAR_MAX; ++i)
      compressed_alphabet[i] = 0;
   for (pi = original_text; pi < original_text + file_size; ++pi)
      compressed_alphabet[*pi] = 1;
   for (i = max_char = 0; i <= UCHAR_MAX; ++i)
      if (compressed_alphabet[i])
         compressed_alphabet[i] = max_char++;
   for (pi = original_text; pi < original_text + file_size; ++pi)
      *pi = compressed_alphabet[*pi] + 1;
   min_char = 1;
   ++max_char;
#endif
#endif
   if (ferror(file)) {
      perror(filename);
      return 1;
   }
#if CHECK || PRINT==2
   s = malloc(file_size * sizeof *s);
   if (! s) {
      fprintf(stderr, "malloc failed\n");
      return 1;
   }
   for (i=0; i< file_size; ++i)
      s[i] = (unsigned char) (original_text[i] - min_char);
#endif

   suffixsort(original_text, suffix_array, file_size, max_char, min_char);

#if CHECK

   fprintf(stderr, "checking...\n");
   for (i=0; i< file_size; ++i) {
      if (scmp3(s + suffix_array[i], s + suffix_array[i + 1], & j, MIN(file_size - suffix_array[i], file_size - suffix_array[i+1])) >= 0)
         fprintf(stderr, "i %d suffix_array[i] %d suffix_array[i+1] %d\n", i, suffix_array[i], suffix_array[i+1]);
   }
   fprintf(stderr, "done.\n");

#endif

#if PRINT

   for (i = 0; i<= file_size; ++i) {
#if PRINT==1
      printf("%d\n", suffix_array[i]);
#else
      printf("%3d \"", suffix_array[i]);
      for (j = suffix_array[i]; j < file_size && j - suffix_array[i] < MAXPRINT; ++j)
         switch(current_char = s[j]) {
         case '\n':
            printf("\\n");
            break;
         case '\t':
            printf("\\t");
            break;
         default:
            putchar(current_char);
         }
      printf("\"\n");
#endif
   }

#endif

   return 0;
}
