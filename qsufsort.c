/* qsufsort.c
   Copyright 1999, N. Jesper Larsson, all rights reserved.

   This file contains an implementation of the algorithm presented in "Faster
   Suffix Sorting" by N. Jesper Larsson (jesper@cs.lth.se) and Kunihiko
   Sadakane (sada@is.s.u-tokyo.ac.jp).

   This software may be used freely for any purpose. However, when distributed,
   the original source must be clearly stated, and, when the source code is
   distributed, the copyright notice must be retained and any alterations in
   the code must be clearly marked. No warranty is given regarding the quality
   of this software.*/

#include <limits.h>

static int *array_buffer,  // end: store the suffix array
   *text_buffer,  // begin: original text; end: inverse of array_buffer
   n_aggregated_chars, // number of symbols aggregated by transform
   n_sorted_prefixes;  // length of already-sorted prefixes.

#define KEY(p)          (text_buffer[*(p)+(n_sorted_prefixes)])
#define SWAP(p, q)      (tmp=*(p), *(p)=*(q), *(q)=tmp)
#define MED3(a, b, c)   (KEY(a)<KEY(b) ?                        \
        (KEY(b)<KEY(c) ? (b) : KEY(a)<KEY(c) ? (c) : (a))       \
        : (KEY(b)>KEY(c) ? (b) : KEY(a)>KEY(c) ? (c) : (a)))

/* Subroutine for select_sort_split and sort_split. Sets group numbers for a
   group whose lowest position in array_buffer is pl and highest position is pm.*/

static void update_group(int *pl, int *pm)
{
   int g;

   g = pm - array_buffer;                      /* group number.*/
   text_buffer[*pl]=g;                    /* update group number of first position.*/
   if (pl==pm)
      *pl=-1;                   /* one element, sorted group.*/
   else
      do                        /* more than one element, unsorted group.*/
         text_buffer[*++pl]=g;            /* update group numbers.*/
      while (pl<pm);
}

/* Quadratic sorting method to use for small subarrays. To be able to update
   group numbers consistently, a variant of selection sorting is used.*/

static void select_sort_split(int *p, int n) {
   int *pa, *pb, *pi, *pn;
   int f, v, tmp;

   pa=p;                        /* pa is start of group being picked out.*/
   pn=p+n-1;                    /* pn is last position of subarray.*/
   while (pa<pn) {
      for (pi=pb=pa+1, f=KEY(pa); pi<=pn; ++pi)
         if ((v=KEY(pi))<f) {
            f=v;                /* f is smallest key found.*/
            SWAP(pi, pa);       /* place smallest element at beginning.*/
            pb=pa+1;            /* pb is position for elements equal to f.*/
         } else if (v==f) {     /* if equal to smallest key.*/
            SWAP(pi, pb);       /* place next to other smallest elements.*/
            ++pb;
         }
      update_group(pa, pb-1);   /* update group values for new group.*/
      pa=pb;                    /* continue sorting rest of the subarray.*/
   }
   if (pa==pn) {                /* check if last part is single element.*/
      text_buffer[*pa] = pa - array_buffer;
      *pa=-1;                   /* sorted group.*/
   }
}

/* Subroutine for sort_split, algorithm by Bentley & McIlroy.*/

static int choose_pivot(int *p, int n) {
   int *pl, *pm, *pn;
   int s;

   pm=p+(n>>1);                 /* small arrays, middle element.*/
   if (n>7) {
      pl=p;
      pn=p+n-1;
      if (n>40) {               /* big arrays, pseudomedian of 9.*/
         s=n>>3;
         pl=MED3(pl, pl+s, pl+s+s);
         pm=MED3(pm-s, pm, pm+s);
         pn=MED3(pn-s-s, pn-s, pn);
      }
      pm=MED3(pl, pm, pn);      /* midsize arrays, median of 3.*/
   }
   return KEY(pm);
}

/* Sorting routine called for each unsorted group. Sorts the array of integers
   (suffix numbers) of length n starting at p. The algorithm is a ternary-split
   quicksort taken from Bentley & McIlroy, "Engineering a Sort Function",
   Software -- Practice and Experience 23(11), 1249-1265 (November 1993). This
   function is based on Program 7.*/

static void sort_split(int *p, int n)
{
   int *pa, *pb, *pc, *pd, *pl, *pm, *pn;
   int f, v, s, t, tmp;

   if (n<7) {                   /* multi-selection sort smallest arrays.*/
      select_sort_split(p, n);
      return;
   }

   v=choose_pivot(p, n);
   pa=pb=p;
   pc=pd=p+n-1;
   while (1) {                  /* split-end partition.*/
      while (pb<=pc && (f=KEY(pb))<=v) {
         if (f==v) {
            SWAP(pa, pb);
            ++pa;
         }
         ++pb;
      }
      while (pc>=pb && (f=KEY(pc))>=v) {
         if (f==v) {
            SWAP(pc, pd);
            --pd;
         }
         --pc;
      }
      if (pb>pc)
         break;
      SWAP(pb, pc);
      ++pb;
      --pc;
   }
   pn=p+n;
   if ((s=pa-p)>(t=pb-pa))
      s=t;
   for (pl=p, pm=pb-s; s; --s, ++pl, ++pm)
      SWAP(pl, pm);
   if ((s=pd-pc)>(t=pn-pd-1))
      s=t;
   for (pl=pb, pm=pn-s; s; --s, ++pl, ++pm)
      SWAP(pl, pm);

   s=pb-pa;
   t=pd-pc;
   if (s>0)
      sort_split(p, s);
   update_group(p+s, p+n-t-1);
   if (t>0)
      sort_split(p+n-t, t);
}

/* Bucketsort for first iteration.

   Input: x[0...n-1] holds integers in the range 1...k-1, all of which appear
   at least once. x[n] is 0. (This is the corresponding output of transform.) k
   must be at most n+1. p is array of size n+1 whose contents are disregarded.

   Output: x is text_buffer and p is array_buffer after the initial sorting stage of the refined
   suffix sorting algorithm.*/

static void bucketsort(int *x, int *p, int n, int k)
{
   int *pi, i, c, d, g;

   for (pi=p; pi<p+k; ++pi)
      *pi=-1;                   /* mark linked lists empty.*/
   for (i=0; i<=n; ++i) {
      x[i]=p[c=x[i]];           /* insert in linked list.*/
      p[c]=i;
   }
   for (pi=p+k-1, i=n; pi>=p; --pi) {
      d=x[c=*pi];               /* c is position, d is next in list.*/
      x[c]=g=i;                 /* last position equals group number.*/
      if (d>=0) {               /* if more than one element in group.*/
         p[i--]=c;              /* p is permutation for the sorted x.*/
         do {
            d=x[c=d];           /* next in linked list.*/
            x[c]=g;             /* group number in x.*/
            p[i--]=c;           /* permutation in p.*/
         } while (d>=0);
      } else
         p[i--]=-1;             /* one element, sorted group.*/
   }
}

/* Transforms the alphabet of original_text by attempting to aggregate several symbols into
   one, while preserving the suffix order of original_text. The alphabet may also be
   compacted, so that original_text on output comprises all integers of the new alphabet
   with no skipped numbers.

   Input: original_text is an array of size original_text_size+1 whose first original_text_size elements are positive
   integers in the range min_char...max_char-1. p is array of size original_text_size+1, used for temporary
   storage. max_char_limit controls aggregation and compaction by defining the maximum value
   for any symbol during transformation: max_char_limit must be at least max_char-l; if max_char_limit <= original_text_size,
   compaction is guaranteed; if max_char-min_char>original_text_size, compaction is never done; if max_char_limit is
   INT_MAX, the maximum number of symbols are aggregated into one.

   Output: Returns an integer size_new_alphabet in the range 1...max_char_limit representing the size of the
   new alphabet. If size_new_alphabet <= original_text_size+1, the alphabet is compacted. The global variable r is
   set to the number of old symbols grouped into one. Only original_text[original_text_size] is 0.*/

//transform(text_buffer, array_buffer, original_text_size, max_char, min_char, original_text_size)
static int transform(int *original_text, int *p, int original_text_size, int max_char, int min_char, int max_char_limit)
{
   int b, c, max_char_new_alphabet, overflow, i, size_new_alphabet, m, n_bits_old_alphabet;
   int *pi, *pj;

   // find highest bit set, based on the range of chars
   for (n_bits_old_alphabet = 0, i = max_char - min_char; i; i >>= 1)
      ++ n_bits_old_alphabet;

   overflow = INT_MAX >> n_bits_old_alphabet; // limits the new alphabet size

   // walk over the input array
   b = 0;
   max_char_new_alphabet = 0;
   for (i = 0; i < original_text_size; ++i) {
      if (max_char_new_alphabet > overflow)
         break;
      c = max_char_new_alphabet << n_bits_old_alphabet | (max_char - min_char);
      if (c > max_char_limit)
         break;
      b = b << n_bits_old_alphabet | (original_text[i] - min_char + 1);        /* b is start of x in chunk alphabet.*/
      max_char_new_alphabet = c;
   }
   n_aggregated_chars = i;

   m = (1 << (n_aggregated_chars - 1) * n_bits_old_alphabet) - 1;            /* m masks off top old symbol from chunk.*/
   original_text[original_text_size] = min_char - 1;  // emulate zero terminator
   if (max_char_new_alphabet <= original_text_size) {                  /* if bucketing possible, compact alphabet.*/
      for (pi = p; pi <= p + max_char_new_alphabet; ++pi)
         *pi = 0;                 /* zero transformation table.*/
      for (pi = original_text + n_aggregated_chars, c = b; pi <= original_text + original_text_size; ++pi) {
         p[c] = 1;                /* mark used chunk symbol.*/
         c = (c & m) << n_bits_old_alphabet | (*pi - min_char + 1);  /* shift in next old symbol in chunk.*/
      }
      for (i = 1; i < n_aggregated_chars; ++i) {     /* handle last n_aggregated_chars-1 positions.*/
         p[c] = 1;                /* mark used chunk symbol.*/
         c = (c & m) << n_bits_old_alphabet;            /* shift in next old symbol in chunk.*/
      }
      for (pi = p, size_new_alphabet = 1; pi <= p + max_char_new_alphabet; ++pi)
         if (*pi)
            *pi = size_new_alphabet ++;
      for (pi = original_text, pj= original_text + n_aggregated_chars, c=b; pj <= original_text + original_text_size; ++pi, ++pj) {
         *pi = p[c];              /* transform to new alphabet.*/
         c = (c&m) << n_bits_old_alphabet | (*pj - min_char + 1);  /* shift in next old symbol in chunk.*/
      }
      while (pi < original_text + original_text_size) {          /* handle last r-1 positions.*/
         *pi++ = p[c];            /* transform to new alphabet.*/
         c = (c & m)<< n_bits_old_alphabet;            /* shift right-end zero in chunk.*/
      }
   } else {                     /* bucketing not possible, don't compact.*/
      for (pi = original_text, pj = original_text + n_aggregated_chars, c=b; pj <= original_text + original_text_size; ++pi, ++pj) {
         *pi=c;                 /* transform to new alphabet.*/
         c = (c & m) << n_bits_old_alphabet | (*pj - min_char + 1);  /* shift in next old symbol in chunk.*/
      }
      while (pi < original_text + original_text_size) {          /* handle last r-1 positions.*/
         *pi ++= c;               /* transform to new alphabet.*/
         c = (c & m) << n_bits_old_alphabet;            /* shift right-end zero in chunk.*/
      }
      size_new_alphabet = max_char_new_alphabet + 1;
   }
   original_text[original_text_size] = 0;                      /* end-of-string symbol is zero.*/
   return size_new_alphabet;
}

/* Makes suffix array p of x. x becomes inverse of p. p and x are both of size
   n+1. Contents of x[0...n-1] are integers in the range l...k-1. Original
   contents of x[n] is disregarded, the n-th symbol being regarded as
   end-of-string smaller than all other symbols.*/

void suffixsort(int *original_text, int *suffix_array, int original_text_size, int max_char, int min_char)
{
   int *pi, *pk;
   int i, alphabet_size, s, sl;

   text_buffer = original_text;
   array_buffer = suffix_array;

   if (original_text_size >= max_char - min_char) {  // if bucket sort makes sense, apply it
      alphabet_size = transform(text_buffer, array_buffer, original_text_size, max_char, min_char, original_text_size);
      bucketsort(text_buffer, array_buffer, original_text_size, alphabet_size);
   }
   else { // otherwise, apply quicksort
      transform(text_buffer, array_buffer, original_text_size, max_char, min_char, INT_MAX);
      for (i=0; i<=original_text_size; ++i)
         array_buffer[i] = i;
      n_sorted_prefixes = 0;
      sort_split(array_buffer, original_text_size + 1);  // quicksort on first n_aggregated_chars positions
   }
   n_sorted_prefixes = n_aggregated_chars;                         /* number of symbols aggregated by transform.*/

   while (*array_buffer >= -original_text_size) {
      pi = array_buffer;                     /* pi is first position of group.*/
      sl=0;                     /* sl is negated length of sorted groups.*/
      do {
         if ((s=*pi)<0) {
            pi -= s;              /* skip over sorted group.*/
            sl += s;              /* add negated length to sl.*/
         } else {
            if (sl) {
               *(pi+sl)=sl;     /* combine sorted groups before pi.*/
               sl=0;
            }
            pk = array_buffer + text_buffer[s]+1;        /* pk-1 is last position of unsorted group.*/
            sort_split(pi, pk-pi);
            pi=pk;              /* next group.*/
         }
      } while (pi <= array_buffer + original_text_size);
      if (sl)                   /* if the array ends with a sorted group.*/
         *(pi+sl)=sl;           /* combine sorted groups at end of array_buffer.*/
      n_sorted_prefixes = 2 * n_sorted_prefixes;
   }

   for (i=0; i<=original_text_size; ++i)         /* reconstruct suffix array from inverse.*/
      array_buffer[text_buffer[i]] = i;
}
