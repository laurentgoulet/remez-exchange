/**************************************************************************
 * Parks-McClellan algorithm for FIR filter design (C version)
 *-------------------------------------------------
 *  Copyright (c) 1995,1998  Jake Janovetz <janovetz@uiuc.edu>
 *
 *  This library is free software; you can redistribute it and/or
 *  modify it under the terms of the GNU Library General Public
 *  License as published by the Free Software Foundation; either
 *  version 2 of the License, or (at your option) any later version.
 *
 *  This library is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 *  Library General Public License for more details.
 *
 *  You should have received a copy of the GNU Library General Public
 *  License along with this library; if not, write to the Free
 *  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 *
 *
 *  Sep 1999 - Paul Kienzle (pkienzle@cs.indiana.edu)
 *      Modified for use in octave as a replacement for the matlab function
 *      remez.mex.  In particular, magnitude responses are required for all
 *      band edges rather than one per band, griddensity is a parameter,
 *      and errors are returned rather than printed directly.
 *  Mar 2000 - Kai Habel (kahacjde@linux.zrz.tu-berlin.de)
 *      Change: ColumnVector x=arg(i).vector_value();
 *      to: ColumnVector x(arg(i).vector_value());
 *  There appear to be some problems with the routine Search. See comments
 *  therein [search for PAK:].  I haven't looked closely at the rest
 *  of the code---it may also have some problems.
 *************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <R.h>

#define CONST const
#define BANDPASS       1
#define DIFFERENTIATOR 2
#define HILBERT        3

#define NEGATIVE       0
#define POSITIVE       1

#define Pi             3.1415926535897932
#define Pi2            6.2831853071795865

#define GRIDDENSITY    16
#define MAXITERATIONS  40

/*******************
 * CreateDenseGrid
 *=================
 * Creates the dense grid of frequencies from the specified bands.
 * Also creates the Desired Frequency Response function (D[]) and
 * the Weight function (W[]) on that dense grid
 *
 *
 * INPUT:
 * ------
 * int      r        - 1/2 the number of filter coefficients
 * int      numtaps  - Number of taps in the resulting filter
 * int      numband  - Number of bands in user specification
 * double   bands[]  - User-specified band edges [2*numband]
 * double   des[]    - Desired response per band [2*numband]
 * double   weight[] - Weight per band [numband]
 * int      symmetry - Symmetry of filter - used for grid check
 * int      griddensity
 *
 * OUTPUT:
 * -------
 * int    gridsize   - Number of elements in the dense frequency grid
 * double Grid[]     - Frequencies (0 to 0.5) on the dense grid [gridsize]
 * double D[]        - Desired response on the dense grid [gridsize]
 * double W[]        - Weight function on the dense grid [gridsize]
 *******************/

void CreateDenseGrid(int r, int numtaps, int numband, const double bands[],
                     const double des[], const double weight[], int* gridsize,
                     double Grid[], double D[], double W[],
                     int symmetry, int griddensity)
{
   int i, j, k, band;
   double delf, lowf, highf, grid0, df, interp;

   /* nominal frequency step */
   delf = 0.5/(griddensity*r);

/*
 * For differentiator, hilbert,
 *   symmetry is odd and Grid[0] = max(0.5*delf, bands[0])
 */
   grid0 = (symmetry == NEGATIVE) && (0.5*delf > bands[0]) ? 0.5*delf : bands[0];

   j=0;
   for (band=0; band < numband; band++)
   {
      lowf = (band==0 ? grid0 : bands[2*band]);
      highf = bands[2*band + 1];
      k = (int)(2*r*griddensity*(bands[2*band+1] - bands[2*band]) + 0.5); /* use same formula as for allocation */
      if (band == 0 && symmetry == NEGATIVE)
        k--;
      /* actual frequency step */
      df = (highf-lowf)/(k-1);
      interp = (des[2*band+1]-des[2*band])/(highf-bands[2*band]);
      for (i=0; i<k; i++)
      {
         D[j] = des[2*band] + (lowf-bands[2*band])*interp;
         W[j] = weight[band];
         Grid[j] = lowf;
         lowf += df;
	 /* Don't exceed allocated gridsize */
         if (j < *gridsize)
         {
           j++;
         }
      }
      Grid[j-1] = highf;
   }

   /* In case of mismatch, propagate initialized gridsize */
   *gridsize = j;
/*
 * Similar to above, if odd symmetry, last grid point can't be .5
 *  - but, if there are even taps, leave the last grid point at .5
 */
   if ((symmetry == NEGATIVE) &&
       (Grid[*gridsize-1] > (0.5 - delf)) &&
       (numtaps % 2))
   {
      Grid[*gridsize-1] = 0.5-delf;
   }
}


/********************
 * InitialGuess
 *==============
 * Places Extremal Frequencies evenly throughout the dense grid.
 *
 *
 * INPUT: 
 * ------
 * int r        - 1/2 the number of filter coefficients
 * int gridsize - Number of elements in the dense frequency grid
 *
 * OUTPUT:
 * -------
 * int Ext[]    - Extremal indexes to dense frequency grid [r+1]
 ********************/

void InitialGuess(int r, int Ext[], int gridsize)
{
   int i;

   for (i=0; i<=r; i++)
      Ext[i] = i * (gridsize-1) / r;
}


/***********************
 * CalcParms
 *===========
 *
 *
 * INPUT:
 * ------
 * int    r      - 1/2 the number of filter coefficients
 * int    Ext[]  - Extremal indexes to dense frequency grid [r+1]
 * double Grid[] - Frequencies (0 to 0.5) on the dense grid [gridsize]
 * double D[]    - Desired response on the dense grid [gridsize]
 * double W[]    - Weight function on the dense grid [gridsize]
 *
 * OUTPUT:
 * -------
 * double ad[]   - 'b' in Oppenheim & Schafer [r+1]
 * double x[]    - [r+1]
 * double y[]    - 'C' in Oppenheim & Schafer [r+1]
 ***********************/

void CalcParms(int r, int Ext[], double Grid[], double D[], double W[],
                double ad[], double x[], double y[])
{
   int i, j, k, ld;
   double sign, xi, delta, denom, numer;

/*
 * Find x[]
 */
   for (i=0; i<=r; i++)
      x[i] = cos(Pi2 * Grid[Ext[i]]);

/*
 * Calculate ad[]  - Oppenheim & Schafer eq 7.132
 */
   ld = (r-1)/15 + 1;         /* Skips around to avoid round errors */
   for (i=0; i<=r; i++)
   {
       denom = 1.0;
       xi = x[i];
       for (j=0; j<ld; j++)
       {
          for (k=j; k<=r; k+=ld)
             if (k != i)
                denom *= 2.0*(xi - x[k]);
       }
       if (fabs(denom)<1e-10)
          denom = 1e-10 * (denom < 0 ? -1 : 1);
       ad[i] = 1.0/denom;
   }

/*
 * Calculate delta  - Oppenheim & Schafer eq 7.131
 */
   numer = denom = 0;
   sign = 1;
   for (i=0; i<=r; i++)
   {
      numer += ad[i] * D[Ext[i]];
      denom += sign * ad[i]/W[Ext[i]];
      sign = -sign;
   }
   delta = numer/denom;
   sign = 1;

/*
 * Calculate y[]  - Oppenheim & Schafer eq 7.133b
 */
   for (i=0; i<=r; i++)
   {
      y[i] = D[Ext[i]] - sign * delta/W[Ext[i]];
      sign = -sign;
   }
}


/*********************
 * ComputeA
 *==========
 * Using values calculated in CalcParms, ComputeA calculates the
 * actual filter response at a given frequency (freq).  Uses
 * eq 7.133a from Oppenheim & Schafer.
 *
 *
 * INPUT:
 * ------
 * double freq - Frequency (0 to 0.5) at which to calculate A
 * int    r    - 1/2 the number of filter coefficients
 * double ad[] - 'b' in Oppenheim & Schafer [r+1]
 * double x[]  - [r+1]
 * double y[]  - 'C' in Oppenheim & Schafer [r+1]
 *
 * OUTPUT:
 * -------
 * Returns double value of A[freq]
 *********************/

double ComputeA(double freq, int r, double ad[], double x[], double y[])
{
   int i;
   double xc, c, denom, numer;

   denom = numer = 0;
   xc = cos(Pi2 * freq);
   for (i=0; i<=r; i++)
   {
      c = xc - x[i];
      if (fabs(c) < 1.0e-7)
      {
         numer = y[i];
         denom = 1;
         break;
      }
      c = ad[i]/c;
      denom += c;
      numer += c*y[i];
   }
   return numer/denom;
}


/************************
 * CalcError
 *===========
 * Calculates the Error function from the desired frequency response
 * on the dense grid (D[]), the weight function on the dense grid (W[]),
 * and the present response calculation (A[])
 *
 *
 * INPUT:
 * ------
 * int    r      - 1/2 the number of filter coefficients
 * double ad[]   - [r+1]
 * double x[]    - [r+1]
 * double y[]    - [r+1]
 * int gridsize  - Number of elements in the dense frequency grid
 * double Grid[] - Frequencies on the dense grid [gridsize]
 * double D[]    - Desired response on the dense grid [gridsize]
 * double W[]    - Weight function on the desnse grid [gridsize]
 *
 * OUTPUT:
 * -------
 * double E[]    - Error function on dense grid [gridsize]
 ************************/

void CalcError(int r, double ad[], double x[], double y[],
               int gridsize, double Grid[],
               double D[], double W[], double E[])
{
   int i;
   double A;

   for (i=0; i<gridsize; i++)
   {
      A = ComputeA(Grid[i], r, ad, x, y);
      E[i] = W[i] * (D[i] - A);
   }
}

/************************
 * Search
 *========
 * Searches for the maxima/minima of the error curve.  If more than
 * r+1 extrema are found, it uses the following heuristic (thanks
 * Chris Hanson):
 * 1) Adjacent non-alternating extrema deleted first.
 * 2) If there are more than one excess extrema, delete the
 *    one with the smallest error.  This will create a non-alternation
 *    condition that is fixed by 1).
 * 3) If there is exactly one excess extremum, delete the smaller
 *    of the first/last extremum
 *
 *
 * INPUT:
 * ------
 * int    r        - 1/2 the number of filter coefficients
 * int    Ext[]    - Indexes to Grid[] of extremal frequencies [r+1]
 * int    gridsize - Number of elements in the dense frequency grid
 * double E[]      - Array of error values.  [gridsize]
 * OUTPUT:
 * -------
 * int    Ext[]    - New indexes to extremal frequencies [r+1]
 ************************/
int Search(int r, int Ext[],
            int gridsize, double E[])
{
   int i, j, k, l, extra;     /* Counters */
   int up, alt;
   int *foundExt;             /* Array of found extremals */

/*
 * Allocate enough space for found extremals.
 */
   foundExt = (int *)malloc((2*r) * sizeof(int));
   k = 0;

/*
 * Check for extremum at 0.
 */
   if (((E[0]>0.0) && (E[0]>E[1])) ||
       ((E[0]<0.0) && (E[0]<E[1])))
      foundExt[k++] = 0;

/*
 * Check for extrema inside dense grid
 */
   for (i=1; i<gridsize-1; i++)
   {
      if (((E[i]>=E[i-1]) && (E[i]>E[i+1]) && (E[i]>0.0)) ||
          ((E[i]<=E[i-1]) && (E[i]<E[i+1]) && (E[i]<0.0))) {
    // PAK: we sometimes get too many extremal frequencies
    if (k >= 2*r) return -3;
    foundExt[k++] = i;
      }
   }

/*
 * Check for extremum at 0.5
 */
   j = gridsize-1;
   if (((E[j]>0.0) && (E[j]>E[j-1])) ||
       ((E[j]<0.0) && (E[j]<E[j-1]))) {
     if (k >= 2*r) return -3;
     foundExt[k++] = j;
   }

   // PAK: we sometimes get not enough extremal frequencies
   if (k < r+1) return -2;


/*
 * Remove extra extremals
 */
   extra = k - (r+1);
   //   assert(extra >= 0);

   while (extra > 0)
   {
      if (E[foundExt[0]] > 0.0)
         up = 1;                /* first one is a maxima */
      else
         up = 0;                /* first one is a minima */

      l=0;
      alt = 1;
      for (j=1; j<k; j++)
      {
         if (fabs(E[foundExt[j]]) < fabs(E[foundExt[l]]))
            l = j;               /* new smallest error. */
         if ((up) && (E[foundExt[j]] < 0.0))
            up = 0;             /* switch to a minima */
         else if ((!up) && (E[foundExt[j]] > 0.0))
            up = 1;             /* switch to a maxima */
         else
     { 
            alt = 0;
        // PAK: break now and you will delete the smallest overall
        // extremal.  If you want to delete the smallest of the
        // pair of non-alternating extremals, then you must do:
            //
        // if (fabs(E[foundExt[j]]) < fabs(E[foundExt[j-1]])) l=j;
        // else l=j-1;
            break;              /* Ooops, found two non-alternating */
         }                      /* extrema.  Delete smallest of them */
      }  /* if the loop finishes, all extrema are alternating */

/*
 * If there's only one extremal and all are alternating,
 * delete the smallest of the first/last extremals.
 */
      if ((alt) && (extra == 1))
      {
         if (fabs(E[foundExt[k-1]]) < fabs(E[foundExt[0]]))
       /* Delete last extremal */
       l = k-1;
       // PAK: changed from l = foundExt[k-1]; 
         else
       /* Delete first extremal */
       l = 0;
       // PAK: changed from l = foundExt[0];     
      }

      for (j=l; j<k-1; j++)        /* Loop that does the deletion */
      {
         foundExt[j] = foundExt[j+1];
         //  assert(foundExt[j]<gridsize);
      }
      k--;
      extra--;
   }

   for (i=0; i<=r; i++)
   {
     //      assert(foundExt[i]<gridsize);
      Ext[i] = foundExt[i];       /* Copy found extremals to Ext[] */
   }

   free(foundExt);
   return 0;
}


/*********************
 * FreqSample
 *============
 * Simple frequency sampling algorithm to determine the impulse
 * response h[] from A's found in ComputeA
 *
 *
 * INPUT:
 * ------
 * int      N        - Number of filter coefficients
 * double   A[]      - Sample points of desired response [N/2]
 * int      symmetry - Symmetry of desired filter
 *
 * OUTPUT:
 * -------
 * double h[] - Impulse Response of final filter [N]
 *********************/
void FreqSample(int N, double A[], double h[], int symm)
{
   int n, k;
   double x, val, M;

   M = (N-1.0)/2.0;
   if (symm == POSITIVE)
   {
      if (N%2)
      {
         for (n=0; n<N; n++)
         {
            val = A[0];
            x = Pi2 * (n - M)/N;
            for (k=1; k<=M; k++)
               val += 2.0 * A[k] * cos(x*k);
            h[n] = val/N;
         }
      }
      else
      {
         for (n=0; n<N; n++)
         {
            val = A[0];
            x = Pi2 * (n - M)/N;
            for (k=1; k<=(N/2-1); k++)
               val += 2.0 * A[k] * cos(x*k);
            h[n] = val/N;
         }
      }
   }
   else
   {
      if (N%2)
      {
         for (n=0; n<N; n++)
         {
            val = 0;
            x = Pi2 * (n - M)/N;
            for (k=1; k<=M; k++)
               val += 2.0 * A[k] * sin(x*k);
            h[n] = -val/N;
         }
      }
      else
      {
          for (n=0; n<N; n++)
          {
             val = A[N/2] * sin(Pi * (n - M));
             x = Pi2 * (n - M)/N;
             for (k=1; k<=(N/2-1); k++)
                val += 2.0 * A[k] * sin(x*k);
             h[n] = -val/N;
          }
      }
   }
}

/*******************
 * isDone
 *========
 * Checks to see if the error function is small enough to consider
 * the result to have converged.
 *
 * INPUT:
 * ------
 * int    r     - 1/2 the number of filter coeffiecients
 * int    Ext[] - Indexes to extremal frequencies [r+1]
 * double E[]   - Error function on the dense grid [gridsize]
 *
 * OUTPUT:
 * -------
 * Returns 1 if the result converged
 * Returns 0 if the result has not converged
 ********************/

int isDone(int r, int Ext[], double E[])
{
   int i;
   double min, max, current;

   min = max = fabs(E[Ext[0]]);
   for (i=1; i<=r; i++)
   {
      current = fabs(E[Ext[i]]);
      if (current < min)
         min = current;
      if (current > max)
         max = current;
   }
   return (((max-min)/max) < 0.0001);
}

/********************
 * remez
 *=======
 * Calculates the optimal (in the Chebyshev/minimax sense)
 * FIR filter impulse response given a set of band edges,
 * the desired reponse on those bands, and the weight given to
 * the error in those bands.
 *
 * INPUT:
 * ------
 * int     *numtaps     - Number of filter coefficients
 * int     *numband     - Number of bands in filter specification
 * double  bands[]      - User-specified band edges [2 * numband]
 * double  des[]        - User-specified band responses [2 * numband]
 * double  weight[]     - User-specified error weights [numband]
 * int     *type        - Type of filter
 * int     *griddensity - ??
 *
 * OUTPUT:
 * -------
 * double h[]      - Impulse response of final filter [numtaps]
 ********************/

void remez(double h[], int *numtaps,
      int *numband, const double bands[], 
      const double des[], const double weight[],
      int *type, int *griddensity)
{
   double *Grid, *W, *D, *E;
   int    i, iter, gridsize, r, *Ext;
   double *taps, c;
   double *x, *y, *ad;
   int    symmetry;

   if (*type == BANDPASS)
      symmetry = POSITIVE;
   else
      symmetry = NEGATIVE;

   r = *numtaps / 2;                  /* number of extrema */
   if ((*numtaps % 2) && (symmetry == POSITIVE))
      r++;
   h[0] = 32;
/*
 * Predict dense grid size in advance for memory allocation
 *   .5 is so we round up, not truncate
 */
   gridsize = 0;
   for (i=0; i < *numband; i++)
   {
      gridsize += (int)(2 * r * (*griddensity) * (bands[2*i+1] - bands[2*i]) + .5);
   }
   if (symmetry == NEGATIVE)
   {
      gridsize--;
   }

/*
 * Dynamically allocate memory for arrays with proper sizes
 */
   Grid = (double *)malloc(gridsize * sizeof(double));
   D = (double *)malloc(gridsize * sizeof(double));
   W = (double *)malloc(gridsize * sizeof(double));
   E = (double *)malloc(gridsize * sizeof(double));
   Ext = (int *)malloc((r+1) * sizeof(int));
   taps = (double *)malloc((r+1) * sizeof(double));
   x = (double *)malloc((r+1) * sizeof(double));
   y = (double *)malloc((r+1) * sizeof(double));
   ad = (double *)malloc((r+1) * sizeof(double));

/*
 * Create dense frequency grid
 */
   CreateDenseGrid(r, *numtaps, *numband, bands, des, weight,
                   &gridsize, Grid, D, W, symmetry, *griddensity);
   InitialGuess(r, Ext, gridsize);

/*
 * For Differentiator: (fix grid)
 */
   if (*type == DIFFERENTIATOR)
   {
      for (i=0; i<gridsize; i++)
      {
/* D[i] = D[i]*Grid[i]; */
         if (D[i] > 0.0001)
            W[i] = W[i]/Grid[i];
      }
   }

/*
 * For odd or Negative symmetry filters, alter the
 * D[] and W[] according to Parks McClellan
 */
   if (symmetry == POSITIVE)
   {
      if (*numtaps % 2 == 0)
      {
         for (i=0; i<gridsize; i++)
         {
            c = cos(Pi * Grid[i]);
            D[i] /= c;
            W[i] *= c; 
         }
      }
   }
   else
   {
      if (*numtaps % 2)
      {
         for (i=0; i<gridsize; i++)
         {
            c = sin(Pi2 * Grid[i]);
            D[i] /= c;
            W[i] *= c;
         }
      }
      else
      {
         for (i=0; i<gridsize; i++)
         {
            c = sin(Pi * Grid[i]);
            D[i] /= c;
            W[i] *= c;
         }
      }
   }

/*
 * Perform the Remez Exchange algorithm
 */
   for (iter=0; iter<MAXITERATIONS; iter++)
   {
      CalcParms(r, Ext, Grid, D, W, ad, x, y);
      CalcError(r, ad, x, y, gridsize, Grid, D, W, E);
      int err = Search(r, Ext, gridsize, E);
      if (err) error("error, %i, %i", err, gridsize);
      //      for(i=0; i <= r; i++) assert(Ext[i]<gridsize);
      if (isDone(r, Ext, E))
         break;
   }

   CalcParms(r, Ext, Grid, D, W, ad, x, y);

/*
 * Find the 'taps' of the filter for use with Frequency
 * Sampling.  If odd or Negative symmetry, fix the taps
 * according to Parks McClellan
 */
   for (i=0; i <= *numtaps / 2; i++)
   {
      if (symmetry == POSITIVE)
      {
         if (*numtaps % 2)
            c = 1;
         else
            c = cos(Pi * (double)i / *numtaps);
      }
      else
      {
         if (*numtaps % 2)
            c = sin(Pi2 * (double)i / *numtaps);
         else
            c = sin(Pi * (double)i / *numtaps);
      }
      taps[i] = ComputeA((double)i / *numtaps, r, ad, x, y) * c;
   }

/*
 * Frequency sampling design with calculated taps
 */
   FreqSample(*numtaps, taps, h, symmetry);

/*
 * Delete allocated memory
 */
   free(Grid);
   free(W);
   free(D);
   free(E);
   free(Ext);
   free(x);
   free(y);
   free(ad);
}


/* == Octave interface starts here ====================================== */
/*******************
DEFUN_DLD (remez, args, ,
  "b = remez(n, f, a [, w] [, ftype] [, griddensity])\n\
Parks-McClellan optimal FIR filter design.\n\
n gives the number of taps in the returned filter\n\
f gives frequency at the band edges [ b1 e1 b2 e2 b3 e3 ...]\n\
a gives amplitude at the band edges [ a(b1) a(e1) a(b2) a(e2) ...]\n\
w gives weighting applied to each band\n\
ftype is 'bandpass', 'hilbert' or 'differentiator'\n\
griddensity determines how accurately the filter will be\n\
    constructed. The minimum value is 16, but higher numbers are\n\
    slower to compute.\n\
\n\
Frequency is in the range (0, 1), with 1 being the nyquist frequency")
{
  octave_value_list retval;
  int i;

  int nargin = args.length();
  if (nargin < 3 || nargin > 6) {
    print_usage("remez");
    return retval;
  }

  int numtaps = NINT (args(0).double_value()) + 1; // #coeff = filter order+1
  if (numtaps < 4) {
    error("remez: number of taps must be an integer greater than 3");
    return retval;
  }

  ColumnVector o_bands(args(1).vector_value());
  int numbands = o_bands.length()/2;
  OCTAVE_LOCAL_BUFFER(double, bands, numbands*2);
  if (numbands < 1 || o_bands.length()%2 == 1) {
    error("remez: must have an even number of band edges");
    return retval;
  }
  for (i=1; i < o_bands.length(); i++) {
    if (o_bands(i)<o_bands(i-1)) {
      error("band edges must be nondecreasing");
      return retval;
    }
  }
  if (o_bands(0) < 0 || o_bands(1) > 1) {
    error("band edges must be in the range [0,1]");
    return retval;
  }
  for(i=0; i < 2*numbands; i++) bands[i] = o_bands(i)/2.0;

  ColumnVector o_response(args(2).vector_value());
  OCTAVE_LOCAL_BUFFER (double, response, numbands*2);
  if (o_response.length() != o_bands.length()) {
    error("remez: must have one response magnitude for each band edge");
    return retval;
  }
  for(i=0; i < 2*numbands; i++) response[i] = o_response(i);

  std::string stype = std::string("bandpass");
  int density = 16;
  OCTAVE_LOCAL_BUFFER (double, weight, numbands);
  for (i=0; i < numbands; i++) weight[i] = 1.0;
  if (nargin > 3) {
    if (args(3).is_real_matrix()) {
      ColumnVector o_weight(args(3).vector_value());
      if (o_weight.length() != numbands) {
    error("remez: need one weight for each band [=length(band)/2]");
    return retval;
      }
      for (i=0; i < numbands; i++) weight[i] = o_weight(i);
    }
    else if (args(3).is_string())
      stype = args(3).string_value();
    else if (args(3).is_real_scalar())
      density = NINT(args(3).double_value());
    else {
      error("remez: incorrect argument list");
      return retval;
    }
  }
  if (nargin > 4) {
    if (args(4).is_string() && !args(3).is_string())
      stype = args(4).string_value();
    else if (args(4).is_real_scalar() && !args(3).is_real_scalar())
      density = NINT(args(4).double_value());
    else {
      error("remez: incorrect argument list");
      return retval;
    }
  }
  if (nargin > 5) {
    if (args(5).is_real_scalar() 
    && !args(4).is_real_scalar() 
    && !args(3).is_real_scalar())
      density = NINT(args(4).double_value());
    else {
      error("remez: incorrect argument list");
      return retval;
    }
  }

  int itype;
  if (stype == "bandpass") 
    itype = BANDPASS;
  else if (stype == "differentiator") 
    itype = DIFFERENTIATOR;
  else if (stype == "hilbert") 
    itype = HILBERT;
  else {
    error("remez: unknown ftype '%s'", stype.data());
    return retval;
  }

  if (density < 16) {
    error("remez: griddensity is too low; must be greater than 16");
    return retval;
  }

  OCTAVE_LOCAL_BUFFER (double, coeff, numtaps+5);
  int err = remez(coeff,numtaps,numbands,bands,response,weight,itype,density);

  if (err == -1)
    warning("remez: -- failed to converge -- returned filter may be bad.");
  else if (err == -2) {
    error("remez: insufficient extremals--cannot continue");
    return retval;
  }
  else if (err == -3) {
    error("remez: too many extremals--cannot continue");
    return retval;
  }

  ColumnVector h(numtaps);
  while(numtaps--) h(numtaps) = coeff[numtaps];

  return octave_value(h);
}
****************/
/*
%!test
%! b = [
%!    0.0416104133896269
%!    0.0581257393237330
%!   -0.0281289583629051
%!   -0.0535690839042500
%!   -0.0617659121718172
%!    0.0507522359162575
%!    0.2078720591778160
%!    0.3326716366583849
%!    0.3326716366583849
%!    0.2078720591778160
%!    0.0507522359162575
%!   -0.0617659121718172
%!   -0.0535690839042500
%!   -0.0281289583629051
%!    0.0581257393237330
%!    0.0416104133896269];
%! assert(remez(15, [0,0.3,0.4,1], [1,1,0,0]), b, 1e-14);

%!test
%! % This test is based on
%! %   "Digital Signal Processing - Principles, Algorithms, and Applications",
%! %   by John G. Proakis and Dimitri G. Manolakis, Third Edition.
%! %   Example 8.2.3 (pp 648-649):
%! %    "Design a lowpass filter of length M = 61 with a passband edge frequency
%! %    f_p = 0.1 and a stopband edge frequency f_s = 0.15."
%! % Note however that this implementation comes within ~1e-5 of the expected
%! % values provided in the reference, but here we match a sample execution
%! % of the implementation to catch subtle variations of the output of
%! % more than 1e-14 precision.
%! b = [
%!   -0.0012125398190740
%!   -0.0006741952600982
%!    0.0000970651050650
%!    0.0013531775338392
%!    0.0022977033858626
%!    0.0019988980782303
%!    0.0001001314465370
%!   -0.0026450655232223
%!   -0.0045143625006237
%!   -0.0037734736912287
%!    0.0000090484550353
%!    0.0051767880724393
%!    0.0084897796253818
%!    0.0069575168320132
%!    0.0000755501040643
%!   -0.0090367960665303
%!   -0.0147234691921247
%!   -0.0119644057670757
%!   -0.0000361740508339
%!    0.0157089621934589
%!    0.0256579017640421
%!    0.0210635925949618
%!    0.0000762130350280
%!   -0.0288966964064803
%!   -0.0491186067913813
%!   -0.0427197958722350
%!   -0.0000578335706833
%!    0.0735682089586026
%!    0.1578205244268486
%!    0.2246608424873270
%!    0.2500781371537121
%!    0.2246608424873270
%!    0.1578205244268486
%!    0.0735682089586026
%!   -0.0000578335706833
%!   -0.0427197958722350
%!   -0.0491186067913813
%!   -0.0288966964064803
%!    0.0000762130350280
%!    0.0210635925949618
%!    0.0256579017640421
%!    0.0157089621934589
%!   -0.0000361740508339
%!   -0.0119644057670757
%!   -0.0147234691921247
%!   -0.0090367960665303
%!    0.0000755501040643
%!    0.0069575168320132
%!    0.0084897796253818
%!    0.0051767880724393
%!    0.0000090484550353
%!   -0.0037734736912287
%!   -0.0045143625006237
%!   -0.0026450655232223
%!    0.0001001314465370
%!    0.0019988980782303
%!    0.0022977033858626
%!    0.0013531775338392
%!    0.0000970651050650
%!   -0.0006741952600982
%!   -0.0012125398190740];
%! assert(remez(60, [0,0.2,0.3,1], [1,1,0,0]), b, 1e-14);

%!test
%! % This test is based on
%! %   "Digital Signal Processing - Principles, Algorithms, and Applications", 
%! %   by John G. Proakis and Dimitri G. Manolakis, Third Edition.
%! %   Example 8.2.4 (pp 650-651)
%! %     "Design a bandpass filter of length M = 32 with a passband edge
%! %      frequencies f_{p1} = 0.2 and a f_{p2} = 0.35 and stopband edge
%! %      frequencies f_{s1} = 0.1 and f_{s2} = 0.425."
%! %   The given solution then also mention the use of a [10,1,10] weights
%! %   on the bands.
%! % Note however that this implementation comes within ~2e-5 of the expected
%! % values provided in the reference, but here we match a sample execution
%! % of the implementation to catch subtle variations of the output of
%! % more than 1e-14 precision.
%! b = [
%!   -0.0057449593160819
%!    0.0009886964787895
%!    0.0075675123109956
%!   -0.0065107266931741
%!    0.0139589163577450
%!    0.0022986511991473
%!   -0.0200018667208152
%!    0.0071313101796869
%!   -0.0396505837710808
%!    0.0112654327508972
%!    0.0662458749269309
%!   -0.0105061313665467
%!    0.0851251628680909
%!   -0.1202492905983940
%!   -0.2968008471295933
%!    0.3041282776771977
%!    0.3041282776771977
%!   -0.2968008471295933
%!   -0.1202492905983940
%!    0.0851251628680909
%!   -0.0105061313665467
%!    0.0662458749269309
%!    0.0112654327508972
%!   -0.0396505837710808
%!    0.0071313101796869
%!   -0.0200018667208152
%!    0.0022986511991473
%!    0.0139589163577450
%!   -0.0065107266931741
%!    0.0075675123109956
%!    0.0009886964787895
%!   -0.0057449593160819];
%! assert(remez(31, [0,0.2,0.4,0.7,0.85,1], [0,0,1,1,0,0], [10,1,10]), b, 1e-14);

%!test
%! % This test is based on
%! %   "Digital Signal Processing - Principles, Algorithms, and Applications", 
%! %   by John G. Proakis and Dimitri G. Manolakis, Third Edition.
%! %   Example 8.2.5 (pp 653-654)
%! %     "Use the Remez algorithm to design a linear-phase FIR differentiator
%! %     of length M = 60. The Passband edge frequency is 0.1 and the
%! %     stopband edge frequency is 0.15."
%! % Note however that this implementation comes within ~5e-6 of the expected
%! % values provided in the reference, but here we match a sample execution
%! % of the implementation to catch subtle variations of the output of more
%! % than 1e-14 precision.
%! b = [
%!   -0.0012472648414296
%!   -0.0015719461409651
%!    0.0036842028016721
%!    0.0019299333612279
%!    0.0014276887563474
%!   -0.0017621204681462
%!   -0.0043098863358648
%!   -0.0046957468293414
%!   -0.0014116599771143
%!    0.0041683503807290
%!    0.0085738512750088
%!    0.0079829000431730
%!    0.0011829082196982
%!   -0.0087370737841627
%!   -0.0154019205012872
%!   -0.0128826131274226
%!   -0.0001870852213829
%!    0.0166197066586607
%!    0.0267395223235475
%!    0.0208938142163349
%!   -0.0018549506376756
%!   -0.0311063074842088
%!   -0.0488259163700931
%!   -0.0386756813189078
%!    0.0036749605876484
%!    0.0654615861680145
%!    0.1206640310148217
%!    0.1418226280689449
%!    0.1140399346967509
%!    0.0436203905283076
%!   -0.0436203905283076
%!   -0.1140399346967509
%!   -0.1418226280689449
%!   -0.1206640310148217
%!   -0.0654615861680145
%!   -0.0036749605876484
%!    0.0386756813189078
%!    0.0488259163700931
%!    0.0311063074842088
%!    0.0018549506376756
%!   -0.0208938142163349
%!   -0.0267395223235475
%!   -0.0166197066586607
%!    0.0001870852213829
%!    0.0128826131274226
%!    0.0154019205012872
%!    0.0087370737841627
%!   -0.0011829082196982
%!   -0.0079829000431730
%!   -0.0085738512750088
%!   -0.0041683503807290
%!    0.0014116599771143
%!    0.0046957468293414
%!    0.0043098863358648
%!    0.0017621204681462
%!   -0.0014276887563474
%!   -0.0019299333612279
%!   -0.0036842028016721
%!    0.0015719461409651
%!    0.0012472648414296];
%! assert(remez(59, [0,0.2,0.3,1], [0,1,0,0], 'differentiator'), b, 1e-14);

%!test
%! % This test is based on
%! %   "Digital Signal Processing - Principles, Algorithms, and Applications", 
%! %   by John G. Proakis and Dimitri G. Manolakis, Third Edition.
%! %   Example 8.2.6 (pp 660-661)
%! %     "Design a Hilbert transformer with parameters M = 31, f_l = 0.05 and
%! %     f_u = 0.45."
%! % Note however that this implementation comes within ~4e-6 of the expected
%! % values provided in the reference, but here we match a sample execution
%! % of the implementation to catch subtle variations of the output of more
%! % than 1e-14 precision.
%! b = [
%!    0.0041996640945240
%!    0.0000000000000010
%!    0.0092829039647483
%!    0.0000000000000012
%!    0.0188353295177154
%!    0.0000000000000016
%!    0.0344008575273274
%!    0.0000000000000015
%!    0.0595539454724730
%!    0.0000000000000013
%!    0.1030380333867637
%!    0.0000000000000010
%!    0.1968306093300132
%!    0.0000000000000006
%!    0.6313549925597921
%!   -0.0000000000000000
%!   -0.6313549925597921
%!   -0.0000000000000006
%!   -0.1968306093300132
%!   -0.0000000000000010
%!   -0.1030380333867637
%!   -0.0000000000000013
%!   -0.0595539454724730
%!   -0.0000000000000015
%!   -0.0344008575273274
%!   -0.0000000000000016
%!   -0.0188353295177154
%!   -0.0000000000000012
%!   -0.0092829039647483
%!   -0.0000000000000010
%!   -0.0041996640945240];
%! assert(remez(30, [0.1,0.9], [1,1], 'hilbert'), b, 1e-14);

%!test
%! % This test is based on
%! %   http://octave.1599824.n4.nabble.com/digital-differentiator-using-remez-td1607186.html
%! % Note however that this implementation comes within ~7e-4 of the expected
%! % values provided from Matlab, but here we match a sample execution of the
%! % implementation to catch subtle variations of the output of more than
%! % 1e-14 precision.
%! b = [
%!   -0.2724192551394910
%!    0.9429429929421833
%!   -0.0000000000000000
%!   -0.9429429929421833
%!    0.2724192551394910];
%! assert(remez(4, [0,0.75], [0,0.75*pi], 'differentiator'), b, 1e-14);

%!test
%! % This test is based on
%! %   http://savannah.gnu.org/bugs/?38134
%! % With some fixes to remez, we now get a descent (frequency response shows
%! % an equiripple relative error of ~27%) differentiator output for an 
%! % order 4 design. Here we match a sample execution of the implementation
%! % to catch subtle variations of the output of more than 1e-14 precision.
%! b = [
%!   -0.3248287372748795
%!    1.0075846739944025
%!   -0.0000000000000000
%!   -1.0075846739944025
%!    0.3248287372748795];
%! assert(remez(4, [0,0.8], [0,0.8*pi], 'differentiator'), b, 1e-14);
%!
%! % To obtain a result similar to quoted result from MATLAB's firpm,
%! % we need to specify an order 6 design (relative error of ~13%), which is
%! % presumably what was used given the number of coefficients.
%! % Note however that this implementation comes within ~3e-4 of the quoted
%! % results from MATLAB (unless we assume that the quoted 0.1519 coefficient
%! % is a typo and should really be 0.1516, in which case the result is really
%! % within ~8e-5).
%! b = [
%!    0.1516425674336023
%!   -0.4161751757098517
%!    0.9424454379608089
%!   -0.0000000000000000
%!   -0.9424454379608089
%!    0.4161751757098517
%!   -0.1516425674336023];
%! assert(remez(6, [0,0.8], [0,0.8*pi], 'differentiator'), b, 1e-14);
 */

