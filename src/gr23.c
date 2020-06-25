/* hwscrt.f -- translated by f2c (version 20031025).
   You must link the resulting object file with libf2c:
	on Microsoft Windows system, link with libf2c.lib;
	on Linux or Unix systems, link with .../path/to/libf2c.a -lm
	or, if you install libf2c.a in a standard place, with -lf2c -lm
	-- in that order, at the end of the command line, as in
		cc *.o -lf2c -lm
	Source for libf2c is in /netlib/f2c/libf2c.zip, e.g.,

		http://www.netlib.org/f2c/libf2c.zip
*/

#include "f2c.h"

/* DECK HWSCRT */
/* Subroutine */ int hwscrt_(real *a, real *b, integer *m, integer *mbdcnd, 
	real *bda, real *bdb, real *c__, real *d__, integer *n, integer *
	nbdcnd, real *bdc, real *bdd, real *elmbda, real *f, integer *idimf, 
	real *pertrb, integer *ierror, real *w)
{
    /* System generated locals */
    integer f_dim1, f_offset, i__1, i__2;
    real r__1;

    /* Local variables */
    static integer i__, j;
    static real s, a1, a2, s1;
    static integer mp, np, id2, id3, id4, mp1, np1;
    static real st2;
    static integer msp1, nsp1, munk, nunk, ierr1, mstm1, nstm1, mskip, nskip, 
	    mstop, nstop;
    extern /* Subroutine */ int genbun_(integer *, integer *, integer *, 
	    integer *, real *, real *, real *, integer *, real *, integer *, 
	    real *);
    static real deltax, deltay;
    static integer mperod, nperod;
    static real delxsq, delysq, twdelx, twdely;
    static integer nstart, mstart;

/* ***BEGIN PROLOGUE  HWSCRT */
/* ***PURPOSE  Solves the standard five-point finite difference */
/*            approximation to the Helmholtz equation in Cartesian */
/*            coordinates. */
/* ***LIBRARY   SLATEC (FISHPACK) */
/* ***CATEGORY  I2B1A1A */
/* ***TYPE      SINGLE PRECISION (HWSCRT-S) */
/* ***KEYWORDS  CARTESIAN, ELLIPTIC, FISHPACK, HELMHOLTZ, PDE */
/* ***AUTHOR  Adams, J., (NCAR) */
/*           Swarztrauber, P. N., (NCAR) */
/*           Sweet, R., (NCAR) */
/* ***DESCRIPTION */

/*     Subroutine HWSCRT solves the standard five-point finite */
/*     difference approximation to the Helmholtz equation in Cartesian */
/*     coordinates: */

/*          (d/dX)(dU/dX) + (d/dY)(dU/dY) + LAMBDA*U = F(X,Y). */



/*     * * * * * * * *    Parameter Description     * * * * * * * * * * */

/*             * * * * * *   On Input    * * * * * * */

/*     A,B */
/*       The range of X, i.e., A .LE. X .LE. B.  A must be less than B. */

/*     M */
/*       The number of panels into which the interval (A,B) is */
/*       subdivided.  Hence, there will be M+1 grid points in the */
/*       X-direction given by X(I) = A+(I-1)DX for I = 1,2,...,M+1, */
/*       where DX = (B-A)/M is the panel width. M must be greater than 3. */

/*     MBDCND */
/*       Indicates the type of boundary conditions at X = A and X = B. */

/*       = 0  If the solution is periodic in X, i.e., U(I,J) = U(M+I,J). */
/*       = 1  If the solution is specified at X = A and X = B. */
/*       = 2  If the solution is specified at X = A and the derivative of */
/*            the solution with respect to X is specified at X = B. */
/*       = 3  If the derivative of the solution with respect to X is */
/*            specified at X = A and X = B. */
/*       = 4  If the derivative of the solution with respect to X is */
/*            specified at X = A and the solution is specified at X = B. */

/*     BDA */
/*       A one-dimensional array of length N+1 that specifies the values */
/*       of the derivative of the solution with respect to X at X = A. */
/*       When MBDCND = 3 or 4, */

/*            BDA(J) = (d/dX)U(A,Y(J)), J = 1,2,...,N+1  . */

/*       When MBDCND has any other value, BDA is a dummy variable. */

/*     BDB */
/*       A one-dimensional array of length N+1 that specifies the values */
/*       of the derivative of the solution with respect to X at X = B. */
/*       When MBDCND = 2 or 3, */

/*            BDB(J) = (d/dX)U(B,Y(J)), J = 1,2,...,N+1  . */

/*       When MBDCND has any other value BDB is a dummy variable. */

/*     C,D */
/*       The range of Y, i.e., C .LE. Y .LE. D.  C must be less than D. */

/*     N */
/*       The number of panels into which the interval (C,D) is */
/*       subdivided.  Hence, there will be N+1 grid points in the */
/*       Y-direction given by Y(J) = C+(J-1)DY for J = 1,2,...,N+1, where */
/*       DY = (D-C)/N is the panel width.  N must be greater than 3. */

/*     NBDCND */
/*       Indicates the type of boundary conditions at Y = C and Y = D. */

/*       = 0  If the solution is periodic in Y, i.e., U(I,J) = U(I,N+J). */
/*       = 1  If the solution is specified at Y = C and Y = D. */
/*       = 2  If the solution is specified at Y = C and the derivative of */
/*            the solution with respect to Y is specified at Y = D. */
/*       = 3  If the derivative of the solution with respect to Y is */
/*            specified at Y = C and Y = D. */
/*       = 4  If the derivative of the solution with respect to Y is */
/*            specified at Y = C and the solution is specified at Y = D. */

/*     BDC */
/*       A one-dimensional array of length M+1 that specifies the values */
/*       of the derivative of the solution with respect to Y at Y = C. */
/*       When NBDCND = 3 or 4, */

/*            BDC(I) = (d/dY)U(X(I),C), I = 1,2,...,M+1  . */

/*       When NBDCND has any other value, BDC is a dummy variable. */

/*     BDD */
/*       A one-dimensional array of length M+1 that specifies the values */
/*       of the derivative of the solution with respect to Y at Y = D. */
/*       When NBDCND = 2 or 3, */

/*            BDD(I) = (d/dY)U(X(I),D), I = 1,2,...,M+1  . */

/*       When NBDCND has any other value, BDD is a dummy variable. */

/*     ELMBDA */
/*       The constant LAMBDA in the Helmholtz equation.  If */
/*       LAMBDA .GT. 0, a solution may not exist.  However, HWSCRT will */
/*       attempt to find a solution. */

/*     F */
/*       A two-dimensional array which specifies the values of the right */
/*       side of the Helmholtz equation and boundary values (if any). */
/*       For I = 2,3,...,M and J = 2,3,...,N */

/*            F(I,J) = F(X(I),Y(J)). */

/*       On the boundaries F is defined by */

/*            MBDCND     F(1,J)        F(M+1,J) */
/*            ------     ---------     -------- */

/*              0        F(A,Y(J))     F(A,Y(J)) */
/*              1        U(A,Y(J))     U(B,Y(J)) */
/*              2        U(A,Y(J))     F(B,Y(J))     J = 1,2,...,N+1 */
/*              3        F(A,Y(J))     F(B,Y(J)) */
/*              4        F(A,Y(J))     U(B,Y(J)) */


/*            NBDCND     F(I,1)        F(I,N+1) */
/*            ------     ---------     -------- */

/*              0        F(X(I),C)     F(X(I),C) */
/*              1        U(X(I),C)     U(X(I),D) */
/*              2        U(X(I),C)     F(X(I),D)     I = 1,2,...,M+1 */
/*              3        F(X(I),C)     F(X(I),D) */
/*              4        F(X(I),C)     U(X(I),D) */

/*       F must be dimensioned at least (M+1)*(N+1). */

/*       NOTE: */

/*       If the table calls for both the solution U and the right side F */
/*       at a corner then the solution must be specified. */

/*     IDIMF */
/*       The row (or first) dimension of the array F as it appears in the */
/*       program calling HWSCRT.  This parameter is used to specify the */
/*       variable dimension of F.  IDIMF must be at least M+1  . */

/*     W */
/*       A one-dimensional array that must be provided by the user for */
/*       work space.  W may require up to 4*(N+1) + */
/*       (13 + INT(log2(N+1)))*(M+1) locations.  The actual number of */
/*       locations used is computed by HWSCRT and is returned in location */
/*       W(1). */


/*             * * * * * *   On Output     * * * * * * */

/*     F */
/*       Contains the solution U(I,J) of the finite difference */
/*       approximation for the grid point (X(I),Y(J)), I = 1,2,...,M+1, */
/*       J = 1,2,...,N+1  . */

/*     PERTRB */
/*       If a combination of periodic or derivative boundary conditions */
/*       is specified for a Poisson equation (LAMBDA = 0), a solution may */
/*       not exist.  PERTRB is a constant, calculated and subtracted from */
/*       F, which ensures that a solution exists.  HWSCRT then computes */
/*       this solution, which is a least squares solution to the original */
/*       approximation.  This solution plus any constant is also a */
/*       solution.  Hence, the solution is not unique.  The value of */
/*       PERTRB should be small compared to the right side F.  Otherwise, */
/*       a solution is obtained to an essentially different problem. */
/*       This comparison should always be made to insure that a */
/*       meaningful solution has been obtained. */

/*     IERROR */
/*       An error flag that indicates invalid input parameters.  Except */
/*       for numbers 0 and 6, a solution is not attempted. */

/*       = 0  No error. */
/*       = 1  A .GE. B. */
/*       = 2  MBDCND .LT. 0 or MBDCND .GT. 4  . */
/*       = 3  C .GE. D. */
/*       = 4  N .LE. 3 */
/*       = 5  NBDCND .LT. 0 or NBDCND .GT. 4  . */
/*       = 6  LAMBDA .GT. 0  . */
/*       = 7  IDIMF .LT. M+1  . */
/*       = 8  M .LE. 3 */

/*       Since this is the only means of indicating a possibly incorrect */
/*       call to HWSCRT, the user should test IERROR after the call. */

/*     W */
/*       W(1) contains the required length of W. */

/* *Long Description: */

/*     * * * * * * *   Program Specifications    * * * * * * * * * * * * */


/*     Dimension of   BDA(N+1),BDB(N+1),BDC(M+1),BDD(M+1),F(IDIMF,N+1), */
/*     Arguments      W(see argument list) */

/*     Latest         June 1, 1976 */
/*     Revision */

/*     Subprograms    HWSCRT,GENBUN,POISD2,POISN2,POISP2,COSGEN,MERGE, */
/*     Required       TRIX,TRI3,PIMACH */

/*     Special        NONE */
/*     Conditions */

/*     Common         NONE */
/*     Blocks */

/*     I/O            NONE */

/*     Precision      Single */

/*     Specialist     Roland Sweet */

/*     Language       FORTRAN */

/*     History        Standardized September 1, 1973 */
/*                    Revised April 1, 1976 */

/*     Algorithm      The routine defines the finite difference */
/*                    equations, incorporates boundary data, and adjusts */
/*                    the right side of singular systems and then calls */
/*                    GENBUN to solve the system. */

/*     Space          13110(octal) = 5704(decimal) locations on the NCAR */
/*     Required       Control Data 7600 */

/*     Timing and        The execution time T on the NCAR Control Data */
/*     Accuracy       7600 for subroutine HWSCRT is roughly proportional */
/*                    to M*N*log2(N), but also depends on the input */
/*                    parameters NBDCND and MBDCND.  Some typical values */
/*                    are listed in the table below. */
/*                       The solution process employed results in a loss */
/*                    of no more than three significant digits for N and */
/*                    M as large as 64.  More detailed information about */
/*                    accuracy can be found in the documentation for */
/*                    subroutine GENBUN which is the routine that */
/*                    solves the finite difference equations. */


/*                       M(=N)    MBDCND    NBDCND    T(MSECS) */
/*                       -----    ------    ------    -------- */

/*                        32        0         0          31 */
/*                        32        1         1          23 */
/*                        32        3         3          36 */
/*                        64        0         0         128 */
/*                        64        1         1          96 */
/*                        64        3         3         142 */

/*     Portability    American National Standards Institute FORTRAN. */
/*                    The machine dependent constant PI is defined in */
/*                    function PIMACH. */

/*     Reference      Swarztrauber, P. and R. Sweet, 'Efficient FORTRAN */
/*                    Subprograms for The Solution Of Elliptic Equations' */
/*                    NCAR TN/IA-109, July, 1975, 138 pp. */

/*     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/* ***REFERENCES  P. N. Swarztrauber and R. Sweet, Efficient Fortran */
/*                 subprograms for the solution of elliptic equations, */
/*                 NCAR TN/IA-109, July 1975, 138 pp. */
/* ***ROUTINES CALLED  GENBUN */
/* ***REVISION HISTORY  (YYMMDD) */
/*   801001  DATE WRITTEN */
/*   890531  Changed all specific intrinsics to generic.  (WRB) */
/*   890531  REVISION DATE from Version 3.2 */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/*   920501  Reformatted the REFERENCES section.  (WRB) */
/* ***END PROLOGUE  HWSCRT */


/* ***FIRST EXECUTABLE STATEMENT  HWSCRT */
    /* Parameter adjustments */
    --bda;
    --bdb;
    --bdc;
    --bdd;
    f_dim1 = *idimf;
    f_offset = 1 + f_dim1;
    f -= f_offset;
    --w;

    /* Function Body */
    *ierror = 0;
    if (*a >= *b) {
	*ierror = 1;
    }
    if (*mbdcnd < 0 || *mbdcnd > 4) {
	*ierror = 2;
    }
    if (*c__ >= *d__) {
	*ierror = 3;
    }
    if (*n <= 3) {
	*ierror = 4;
    }
    if (*nbdcnd < 0 || *nbdcnd > 4) {
	*ierror = 5;
    }
    if (*idimf < *m + 1) {
	*ierror = 7;
    }
    if (*m <= 3) {
	*ierror = 8;
    }
    if (*ierror != 0) {
	return 0;
    }
    nperod = *nbdcnd;
    mperod = 0;
    if (*mbdcnd > 0) {
	mperod = 1;
    }
    deltax = (*b - *a) / *m;
    twdelx = 2.f / deltax;
/* Computing 2nd power */
    r__1 = deltax;
    delxsq = 1.f / (r__1 * r__1);
    deltay = (*d__ - *c__) / *n;
    twdely = 2.f / deltay;
/* Computing 2nd power */
    r__1 = deltay;
    delysq = 1.f / (r__1 * r__1);
    np = *nbdcnd + 1;
    np1 = *n + 1;
    mp = *mbdcnd + 1;
    mp1 = *m + 1;
    nstart = 1;
    nstop = *n;
    nskip = 1;
    switch (np) {
	case 1:  goto L104;
	case 2:  goto L101;
	case 3:  goto L102;
	case 4:  goto L103;
	case 5:  goto L104;
    }
L101:
    nstart = 2;
    goto L104;
L102:
    nstart = 2;
L103:
    nstop = np1;
    nskip = 2;
L104:
    nunk = nstop - nstart + 1;

/*     ENTER BOUNDARY DATA FOR X-BOUNDARIES. */

    mstart = 1;
    mstop = *m;
    mskip = 1;
    switch (mp) {
	case 1:  goto L117;
	case 2:  goto L105;
	case 3:  goto L106;
	case 4:  goto L109;
	case 5:  goto L110;
    }
L105:
    mstart = 2;
    goto L107;
L106:
    mstart = 2;
    mstop = mp1;
    mskip = 2;
L107:
    i__1 = nstop;
    for (j = nstart; j <= i__1; ++j) {
	f[j * f_dim1 + 2] -= f[j * f_dim1 + 1] * delxsq;
/* L108: */
    }
    goto L112;
L109:
    mstop = mp1;
    mskip = 2;
L110:
    i__1 = nstop;
    for (j = nstart; j <= i__1; ++j) {
	f[j * f_dim1 + 1] += bda[j] * twdelx;
/* L111: */
    }
L112:
    switch (mskip) {
	case 1:  goto L113;
	case 2:  goto L115;
    }
L113:
    i__1 = nstop;
    for (j = nstart; j <= i__1; ++j) {
	f[*m + j * f_dim1] -= f[mp1 + j * f_dim1] * delxsq;
/* L114: */
    }
    goto L117;
L115:
    i__1 = nstop;
    for (j = nstart; j <= i__1; ++j) {
	f[mp1 + j * f_dim1] -= bdb[j] * twdelx;
/* L116: */
    }
L117:
    munk = mstop - mstart + 1;

/*     ENTER BOUNDARY DATA FOR Y-BOUNDARIES. */

    switch (np) {
	case 1:  goto L127;
	case 2:  goto L118;
	case 3:  goto L118;
	case 4:  goto L120;
	case 5:  goto L120;
    }
L118:
    i__1 = mstop;
    for (i__ = mstart; i__ <= i__1; ++i__) {
	f[i__ + (f_dim1 << 1)] -= f[i__ + f_dim1] * delysq;
/* L119: */
    }
    goto L122;
L120:
    i__1 = mstop;
    for (i__ = mstart; i__ <= i__1; ++i__) {
	f[i__ + f_dim1] += bdc[i__] * twdely;
/* L121: */
    }
L122:
    switch (nskip) {
	case 1:  goto L123;
	case 2:  goto L125;
    }
L123:
    i__1 = mstop;
    for (i__ = mstart; i__ <= i__1; ++i__) {
	f[i__ + *n * f_dim1] -= f[i__ + np1 * f_dim1] * delysq;
/* L124: */
    }
    goto L127;
L125:
    i__1 = mstop;
    for (i__ = mstart; i__ <= i__1; ++i__) {
	f[i__ + np1 * f_dim1] -= bdd[i__] * twdely;
/* L126: */
    }

/*    MULTIPLY RIGHT SIDE BY DELTAY**2. */

L127:
    delysq = deltay * deltay;
    i__1 = mstop;
    for (i__ = mstart; i__ <= i__1; ++i__) {
	i__2 = nstop;
	for (j = nstart; j <= i__2; ++j) {
	    f[i__ + j * f_dim1] *= delysq;
/* L128: */
	}
/* L129: */
    }

/*     DEFINE THE A,B,C COEFFICIENTS IN W-ARRAY. */

    id2 = munk;
    id3 = id2 + munk;
    id4 = id3 + munk;
    s = delysq * delxsq;
    st2 = s * 2.f;
    i__1 = munk;
    for (i__ = 1; i__ <= i__1; ++i__) {
	w[i__] = s;
	j = id2 + i__;
	w[j] = -st2 + *elmbda * delysq;
	j = id3 + i__;
	w[j] = s;
/* L130: */
    }
    if (mp == 1) {
	goto L131;
    }
    w[1] = 0.f;
    w[id4] = 0.f;
L131:
    switch (mp) {
	case 1:  goto L135;
	case 2:  goto L135;
	case 3:  goto L132;
	case 4:  goto L133;
	case 5:  goto L134;
    }
L132:
    w[id2] = st2;
    goto L135;
L133:
    w[id2] = st2;
L134:
    w[id3 + 1] = st2;
L135:
    *pertrb = 0.f;
    if (*elmbda < 0.f) {
	goto L144;
    } else if (*elmbda == 0) {
	goto L137;
    } else {
	goto L136;
    }
L136:
    *ierror = 6;
    goto L144;
L137:
    if ((*nbdcnd == 0 || *nbdcnd == 3) && (*mbdcnd == 0 || *mbdcnd == 3)) {
	goto L138;
    }
    goto L144;

/*     FOR SINGULAR PROBLEMS MUST ADJUST DATA TO INSURE THAT A SOLUTION */
/*     WILL EXIST. */

L138:
    a1 = 1.f;
    a2 = 1.f;
    if (*nbdcnd == 3) {
	a2 = 2.f;
    }
    if (*mbdcnd == 3) {
	a1 = 2.f;
    }
    s1 = 0.f;
    msp1 = mstart + 1;
    mstm1 = mstop - 1;
    nsp1 = nstart + 1;
    nstm1 = nstop - 1;
    i__1 = nstm1;
    for (j = nsp1; j <= i__1; ++j) {
	s = 0.f;
	i__2 = mstm1;
	for (i__ = msp1; i__ <= i__2; ++i__) {
	    s += f[i__ + j * f_dim1];
/* L139: */
	}
	s1 = s1 + s * a1 + f[mstart + j * f_dim1] + f[mstop + j * f_dim1];
/* L140: */
    }
    s1 = a2 * s1;
    s = 0.f;
    i__1 = mstm1;
    for (i__ = msp1; i__ <= i__1; ++i__) {
	s = s + f[i__ + nstart * f_dim1] + f[i__ + nstop * f_dim1];
/* L141: */
    }
    s1 = s1 + s * a1 + f[mstart + nstart * f_dim1] + f[mstart + nstop * 
	    f_dim1] + f[mstop + nstart * f_dim1] + f[mstop + nstop * f_dim1];
    s = ((nunk - 2) * a2 + 2.f) * ((munk - 2) * a1 + 2.f);
    *pertrb = s1 / s;
    i__1 = nstop;
    for (j = nstart; j <= i__1; ++j) {
	i__2 = mstop;
	for (i__ = mstart; i__ <= i__2; ++i__) {
	    f[i__ + j * f_dim1] -= *pertrb;
/* L142: */
	}
/* L143: */
    }
    *pertrb /= delysq;

/*     SOLVE THE EQUATION. */

L144:
    genbun_(&nperod, &nunk, &mperod, &munk, &w[1], &w[id2 + 1], &w[id3 + 1], 
	    idimf, &f[mstart + nstart * f_dim1], &ierr1, &w[id4 + 1]);
    w[1] = w[id4 + 1] + munk * 3;

/*     FILL IN IDENTICAL VALUES WHEN HAVE PERIODIC BOUNDARY CONDITIONS. */

    if (*nbdcnd != 0) {
	goto L146;
    }
    i__1 = mstop;
    for (i__ = mstart; i__ <= i__1; ++i__) {
	f[i__ + np1 * f_dim1] = f[i__ + f_dim1];
/* L145: */
    }
L146:
    if (*mbdcnd != 0) {
	goto L148;
    }
    i__1 = nstop;
    for (j = nstart; j <= i__1; ++j) {
	f[mp1 + j * f_dim1] = f[j * f_dim1 + 1];
/* L147: */
    }
    if (*nbdcnd == 0) {
	f[mp1 + np1 * f_dim1] = f[np1 * f_dim1 + 1];
    }
L148:
    return 0;
} /* hwscrt_ */

/* genbun.f -- translated by f2c (version 20031025).
   You must link the resulting object file with libf2c:
	on Microsoft Windows system, link with libf2c.lib;
	on Linux or Unix systems, link with .../path/to/libf2c.a -lm
	or, if you install libf2c.a in a standard place, with -lf2c -lm
	-- in that order, at the end of the command line, as in
		cc *.o -lf2c -lm
	Source for libf2c is in /netlib/f2c/libf2c.zip, e.g.,

		http://www.netlib.org/f2c/libf2c.zip
*/

#include "f2c.h"

/* Table of constant values */

static integer c__1 = 1;
static integer c__2 = 2;

/* DECK GENBUN */
/* Subroutine */ int genbun_(integer *nperod, integer *n, integer *mperod, 
	integer *m, real *a, real *b, real *c__, integer *idimy, real *y, 
	integer *ierror, real *w)
{
    /* System generated locals */
    integer y_dim1, y_offset, i__1, i__2;

    /* Local variables */
    static integer i__, j, k;
    static real a1;
    static integer mh, mp, np, mp1, iwd, iwp, mhm1, iwb2, iwb3, nby2, iww1, 
	    iww2, iww3, iwba, iwbb, iwbc, modd, mhmi, mhpi, irev, mskip;
    extern /* Subroutine */ int poisd2_(integer *, integer *, integer *, real 
	    *, real *, real *, real *, integer *, real *, real *, real *, 
	    real *, real *), poisn2_(integer *, integer *, integer *, integer 
	    *, real *, real *, real *, real *, integer *, real *, real *, 
	    real *, real *, real *, real *, real *, real *, real *), poisp2_(
	    integer *, integer *, real *, real *, real *, real *, integer *, 
	    real *, real *, real *, real *, real *, real *, real *, real *, 
	    real *);
    static integer iwtcos, ipstor;

/* ***BEGIN PROLOGUE  GENBUN */
/* ***PURPOSE  Solve by a cyclic reduction algorithm the linear system */
/*            of equations that results from a finite difference */
/*            approximation to certain 2-d elliptic PDE's on a centered */
/*            grid . */
/* ***LIBRARY   SLATEC (FISHPACK) */
/* ***CATEGORY  I2B4B */
/* ***TYPE      SINGLE PRECISION (GENBUN-S, CMGNBN-C) */
/* ***KEYWORDS  ELLIPTIC, FISHPACK, PDE, TRIDIAGONAL */
/* ***AUTHOR  Adams, J., (NCAR) */
/*           Swarztrauber, P. N., (NCAR) */
/*           Sweet, R., (NCAR) */
/* ***DESCRIPTION */

/*     Subroutine GENBUN solves the linear system of equations */

/*          A(I)*X(I-1,J) + B(I)*X(I,J) + C(I)*X(I+1,J) */

/*          + X(I,J-1) - 2.*X(I,J) + X(I,J+1) = Y(I,J) */

/*               for I = 1,2,...,M  and  J = 1,2,...,N. */

/*     The indices I+1 and I-1 are evaluated modulo M, i.e., */
/*     X(0,J) = X(M,J) and X(M+1,J) = X(1,J), and X(I,0) may be equal to */
/*     0, X(I,2), or X(I,N) and X(I,N+1) may be equal to 0, X(I,N-1), or */
/*     X(I,1) depending on an input parameter. */


/*     * * * * * * * *    Parameter Description     * * * * * * * * * * */

/*             * * * * * *   On Input    * * * * * * */

/*     NPEROD */
/*       Indicates the values that X(I,0) and X(I,N+1) are assumed to */
/*       have. */

/*       = 0  If X(I,0) = X(I,N) and X(I,N+1) = X(I,1). */
/*       = 1  If X(I,0) = X(I,N+1) = 0  . */
/*       = 2  If X(I,0) = 0 and X(I,N+1) = X(I,N-1). */
/*       = 3  If X(I,0) = X(I,2) and X(I,N+1) = X(I,N-1). */
/*       = 4  If X(I,0) = X(I,2) and X(I,N+1) = 0. */

/*     N */
/*       The number of unknowns in the J-direction.  N must be greater */
/*       than 2. */

/*     MPEROD */
/*       = 0 if A(1) and C(M) are not zero. */
/*       = 1 if A(1) = C(M) = 0. */

/*     M */
/*       The number of unknowns in the I-direction.  M must be greater */
/*       than 2. */

/*     A,B,C */
/*       One-dimensional arrays of length M that specify the */
/*       coefficients in the linear equations given above.  If MPEROD = 0 */
/*       the array elements must not depend upon the index I, but must be */
/*       constant.  Specifically, the subroutine checks the following */
/*       condition */

/*             A(I) = C(1) */
/*             C(I) = C(1) */
/*             B(I) = B(1) */

/*       for I=1,2,...,M. */

/*     IDIMY */
/*       The row (or first) dimension of the two-dimensional array Y as */
/*       it appears in the program calling GENBUN.  This parameter is */
/*       used to specify the variable dimension of Y.  IDIMY must be at */
/*       least M. */

/*     Y */
/*       A two-dimensional array that specifies the values of the right */
/*       side of the linear system of equations given above.  Y must be */
/*       dimensioned at least M*N. */

/*     W */
/*       A one-dimensional array that must be provided by the user for */
/*       work space.  W may require up to 4*N + (10 + INT(log2(N)))*M */
/*       locations.  The actual number of locations used is computed by */
/*       GENBUN and is returned in location W(1). */


/*             * * * * * *   On Output     * * * * * * */

/*     Y */
/*       Contains the solution X. */

/*     IERROR */
/*       An error flag that indicates invalid input parameters.  Except */
/*       for number zero, a solution is not attempted. */

/*       = 0  No error. */
/*       = 1  M .LE. 2 */
/*       = 2  N .LE. 2 */
/*       = 3  IDIMY .LT. M */
/*       = 4  NPEROD .LT. 0 or NPEROD .GT. 4 */
/*       = 5  MPEROD .LT. 0 or MPEROD .GT. 1 */
/*       = 6  A(I) .NE. C(1) or C(I) .NE. C(1) or B(I) .NE. B(1) for */
/*            some I=1,2,...,M. */
/*       = 7  A(1) .NE. 0 or C(M) .NE. 0 and MPEROD = 1 */

/*     W */
/*       W(1) contains the required length of W. */

/* *Long Description: */

/*     * * * * * * *   Program Specifications    * * * * * * * * * * * * */

/*     Dimension of   A(M),B(M),C(M),Y(IDIMY,N),W(see parameter list) */
/*     Arguments */

/*     Latest         June 1, 1976 */
/*     Revision */

/*     Subprograms    GENBUN,POISD2,POISN2,POISP2,COSGEN,MERGE,TRIX,TRI3, */
/*     Required       PIMACH */

/*     Special        NONE */
/*     Conditions */

/*     Common         NONE */
/*     Blocks */

/*     I/O            NONE */

/*     Precision      Single */

/*     Specialist     Roland Sweet */

/*     Language       FORTRAN */

/*     History        Standardized April 1, 1973 */
/*                    Revised August 20,1973 */
/*                    Revised January 1, 1976 */

/*     Algorithm      The linear system is solved by a cyclic reduction */
/*                    algorithm described in the reference. */

/*     Space          4944(decimal) = 11520(octal) locations on the NCAR */
/*     Required       Control Data 7600. */

/*     Timing and        The execution time T on the NCAR Control Data */
/*     Accuracy       7600 for subroutine GENBUN is roughly proportional */
/*                    to M*N*log2(N), but also depends on the input */
/*                    parameter NPEROD.  Some typical values are listed */
/*                    in the table below.  More comprehensive timing */
/*                    charts may be found in the reference. */
/*                       To measure the accuracy of the algorithm a */
/*                    uniform random number generator was used to create */
/*                    a solution array X for the system given in the */
/*                    'PURPOSE' with */

/*                       A(I) = C(I) = -0.5*B(I) = 1,       I=1,2,...,M */

/*                    and, when MPEROD = 1 */

/*                       A(1) = C(M) = 0 */
/*                       A(M) = C(1) = 2. */

/*                    The solution X was substituted into the given sys- */
/*                    tem and, using double precision, a right side Y was */
/*                    computed.  Using this array Y subroutine GENBUN was */
/*                    called to produce an approximate solution Z.  Then */
/*                    the relative error, defined as */

/*                       E = MAX(ABS(Z(I,J)-X(I,J)))/MAX(ABS(X(I,J))) */

/*                    where the two maxima are taken over all I=1,2,...,M */
/*                    and J=1,2,...,N, was computed.  The value of E is */
/*                    given in the table below for some typical values of */
/*                    M and N. */


/*                       M (=N)    MPEROD    NPEROD    T(MSECS)    E */
/*                       ------    ------    ------    --------  ------ */

/*                         31        0         0          36     6.E-14 */
/*                         31        1         1          21     4.E-13 */
/*                         31        1         3          41     3.E-13 */
/*                         32        0         0          29     9.E-14 */
/*                         32        1         1          32     3.E-13 */
/*                         32        1         3          48     1.E-13 */
/*                         33        0         0          36     9.E-14 */
/*                         33        1         1          30     4.E-13 */
/*                         33        1         3          34     1.E-13 */
/*                         63        0         0         150     1.E-13 */
/*                         63        1         1          91     1.E-12 */
/*                         63        1         3         173     2.E-13 */
/*                         64        0         0         122     1.E-13 */
/*                         64        1         1         128     1.E-12 */
/*                         64        1         3         199     6.E-13 */
/*                         65        0         0         143     2.E-13 */
/*                         65        1         1         120     1.E-12 */
/*                         65        1         3         138     4.E-13 */

/*     Portability    American National Standards Institute Fortran. */
/*                    The machine dependent constant PI is defined in */
/*                    function PIMACH. */

/*     Required       COS */
/*     Resident */
/*     Routines */

/*     Reference      Sweet, R., 'A Cyclic Reduction Algorithm For */
/*                    Solving Block Tridiagonal Systems Of Arbitrary */
/*                    Dimensions,' SIAM J. on Numer. Anal., */
/*                    14(Sept., 1977), PP. 706-720. */

/*     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/* ***REFERENCES  R. Sweet, A cyclic reduction algorithm for solving */
/*                 block tridiagonal systems of arbitrary dimensions, */
/*                 SIAM Journal on Numerical Analysis 14, (September */
/*                 1977), pp. 706-720. */
/* ***ROUTINES CALLED  POISD2, POISN2, POISP2 */
/* ***REVISION HISTORY  (YYMMDD) */
/*   801001  DATE WRITTEN */
/*   861211  REVISION DATE from Version 3.2 */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/*   920501  Reformatted the REFERENCES section.  (WRB) */
/* ***END PROLOGUE  GENBUN */


/* ***FIRST EXECUTABLE STATEMENT  GENBUN */
    /* Parameter adjustments */
    --a;
    --b;
    --c__;
    y_dim1 = *idimy;
    y_offset = 1 + y_dim1;
    y -= y_offset;
    --w;

    /* Function Body */
    *ierror = 0;
    if (*m <= 2) {
	*ierror = 1;
    }
    if (*n <= 2) {
	*ierror = 2;
    }
    if (*idimy < *m) {
	*ierror = 3;
    }
    if (*nperod < 0 || *nperod > 4) {
	*ierror = 4;
    }
    if (*mperod < 0 || *mperod > 1) {
	*ierror = 5;
    }
    if (*mperod == 1) {
	goto L102;
    }
    i__1 = *m;
    for (i__ = 2; i__ <= i__1; ++i__) {
	if (a[i__] != c__[1]) {
	    goto L103;
	}
	if (c__[i__] != c__[1]) {
	    goto L103;
	}
	if (b[i__] != b[1]) {
	    goto L103;
	}
/* L101: */
    }
    goto L104;
L102:
    if (a[1] != 0.f || c__[*m] != 0.f) {
	*ierror = 7;
    }
    goto L104;
L103:
    *ierror = 6;
L104:
    if (*ierror != 0) {
	return 0;
    }
    mp1 = *m + 1;
    iwba = mp1;
    iwbb = iwba + *m;
    iwbc = iwbb + *m;
    iwb2 = iwbc + *m;
    iwb3 = iwb2 + *m;
    iww1 = iwb3 + *m;
    iww2 = iww1 + *m;
    iww3 = iww2 + *m;
    iwd = iww3 + *m;
    iwtcos = iwd + *m;
    iwp = iwtcos + (*n << 2);
    i__1 = *m;
    for (i__ = 1; i__ <= i__1; ++i__) {
	k = iwba + i__ - 1;
	w[k] = -a[i__];
	k = iwbc + i__ - 1;
	w[k] = -c__[i__];
	k = iwbb + i__ - 1;
	w[k] = 2.f - b[i__];
	i__2 = *n;
	for (j = 1; j <= i__2; ++j) {
	    y[i__ + j * y_dim1] = -y[i__ + j * y_dim1];
/* L105: */
	}
/* L106: */
    }
    mp = *mperod + 1;
    np = *nperod + 1;
    switch (mp) {
	case 1:  goto L114;
	case 2:  goto L107;
    }
L107:
    switch (np) {
	case 1:  goto L108;
	case 2:  goto L109;
	case 3:  goto L110;
	case 4:  goto L111;
	case 5:  goto L123;
    }
L108:
    poisp2_(m, n, &w[iwba], &w[iwbb], &w[iwbc], &y[y_offset], idimy, &w[1], &
	    w[iwb2], &w[iwb3], &w[iww1], &w[iww2], &w[iww3], &w[iwd], &w[
	    iwtcos], &w[iwp]);
    goto L112;
L109:
    poisd2_(m, n, &c__1, &w[iwba], &w[iwbb], &w[iwbc], &y[y_offset], idimy, &
	    w[1], &w[iww1], &w[iwd], &w[iwtcos], &w[iwp]);
    goto L112;
L110:
    poisn2_(m, n, &c__1, &c__2, &w[iwba], &w[iwbb], &w[iwbc], &y[y_offset], 
	    idimy, &w[1], &w[iwb2], &w[iwb3], &w[iww1], &w[iww2], &w[iww3], &
	    w[iwd], &w[iwtcos], &w[iwp]);
    goto L112;
L111:
    poisn2_(m, n, &c__1, &c__1, &w[iwba], &w[iwbb], &w[iwbc], &y[y_offset], 
	    idimy, &w[1], &w[iwb2], &w[iwb3], &w[iww1], &w[iww2], &w[iww3], &
	    w[iwd], &w[iwtcos], &w[iwp]);
L112:
    ipstor = w[iww1];
    irev = 2;
    if (*nperod == 4) {
	goto L124;
    }
L113:
    switch (mp) {
	case 1:  goto L127;
	case 2:  goto L133;
    }
L114:

/*     REORDER UNKNOWNS WHEN MP =0 */

    mh = (*m + 1) / 2;
    mhm1 = mh - 1;
    modd = 1;
    if (mh << 1 == *m) {
	modd = 2;
    }
    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
	i__2 = mhm1;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    mhpi = mh + i__;
	    mhmi = mh - i__;
	    w[i__] = y[mhmi + j * y_dim1] - y[mhpi + j * y_dim1];
	    w[mhpi] = y[mhmi + j * y_dim1] + y[mhpi + j * y_dim1];
/* L115: */
	}
	w[mh] = y[mh + j * y_dim1] * 2.f;
	switch (modd) {
	    case 1:  goto L117;
	    case 2:  goto L116;
	}
L116:
	w[*m] = y[*m + j * y_dim1] * 2.f;
L117:
	i__2 = *m;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    y[i__ + j * y_dim1] = w[i__];
/* L118: */
	}
/* L119: */
    }
    k = iwbc + mhm1 - 1;
    i__ = iwba + mhm1;
    w[k] = 0.f;
    w[i__] = 0.f;
    w[k + 1] *= 2.f;
    switch (modd) {
	case 1:  goto L120;
	case 2:  goto L121;
    }
L120:
    k = iwbb + mhm1 - 1;
    w[k] -= w[i__ - 1];
    w[iwbc - 1] += w[iwbb - 1];
    goto L122;
L121:
    w[iwbb - 1] = w[k + 1];
L122:
    goto L107;

/*     REVERSE COLUMNS WHEN NPEROD = 4. */

L123:
    irev = 1;
    nby2 = *n / 2;
L124:
    i__1 = nby2;
    for (j = 1; j <= i__1; ++j) {
	mskip = *n + 1 - j;
	i__2 = *m;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    a1 = y[i__ + j * y_dim1];
	    y[i__ + j * y_dim1] = y[i__ + mskip * y_dim1];
	    y[i__ + mskip * y_dim1] = a1;
/* L125: */
	}
/* L126: */
    }
    switch (irev) {
	case 1:  goto L110;
	case 2:  goto L113;
    }
L127:
    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
	i__2 = mhm1;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    mhmi = mh - i__;
	    mhpi = mh + i__;
	    w[mhmi] = (y[mhpi + j * y_dim1] + y[i__ + j * y_dim1]) * .5f;
	    w[mhpi] = (y[mhpi + j * y_dim1] - y[i__ + j * y_dim1]) * .5f;
/* L128: */
	}
	w[mh] = y[mh + j * y_dim1] * .5f;
	switch (modd) {
	    case 1:  goto L130;
	    case 2:  goto L129;
	}
L129:
	w[*m] = y[*m + j * y_dim1] * .5f;
L130:
	i__2 = *m;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    y[i__ + j * y_dim1] = w[i__];
/* L131: */
	}
/* L132: */
    }
L133:

/*     RETURN STORAGE REQUIREMENTS FOR W ARRAY. */

    w[1] = (real) (ipstor + iwp - 1);
    return 0;
} /* genbun_ */

/* poisd2.f -- translated by f2c (version 20031025).
   You must link the resulting object file with libf2c:
	on Microsoft Windows system, link with libf2c.lib;
	on Linux or Unix systems, link with .../path/to/libf2c.a -lm
	or, if you install libf2c.a in a standard place, with -lf2c -lm
	-- in that order, at the end of the command line, as in
		cc *.o -lf2c -lm
	Source for libf2c is in /netlib/f2c/libf2c.zip, e.g.,

		http://www.netlib.org/f2c/libf2c.zip
*/

#include "f2c.h"

/* Table of constant values */

static integer c__0 = 0;
static real c_b17 = .5f;
static real c_b18 = 0.f;

/* DECK POISD2 */
/* Subroutine */ int poisd2_(integer *mr, integer *nr, integer *istag, real *
	ba, real *bb, real *bc, real *q, integer *idimq, real *b, real *w, 
	real *d__, real *tcos, real *p)
{
    /* System generated locals */
    integer q_dim1, q_offset, i__1, i__2, i__3;

    /* Local variables */
    static integer i__, j, l, m, n;
    static real t, fi;
    static integer ip, kr, lr, jm1, jm2, jm3, jp1, jp2, jp3, ip1, jsh, jsp, 
	    nun, jst, ideg, jdeg, nodd, krpi;
    extern /* Subroutine */ int trix_(integer *, integer *, integer *, real *,
	     real *, real *, real *, real *, real *, real *);
    static integer irreg;
    extern /* Subroutine */ int s1merg_(real *, integer *, integer *, integer 
	    *, integer *, integer *), cosgen_(integer *, integer *, real *, 
	    real *, real *);
    static integer noddpr, jstsav, ipstor;

/* ***BEGIN PROLOGUE  POISD2 */
/* ***SUBSIDIARY */
/* ***PURPOSE  Subsidiary to GENBUN */
/* ***LIBRARY   SLATEC */
/* ***TYPE      SINGLE PRECISION (POISD2-S, CMPOSD-C) */
/* ***AUTHOR  (UNKNOWN) */
/* ***DESCRIPTION */

/*     Subroutine to solve Poisson's equation for Dirichlet boundary */
/*     conditions. */

/*     ISTAG = 1 if the last diagonal block is the matrix A. */
/*     ISTAG = 2 if the last diagonal block is the matrix A+I. */

/* ***SEE ALSO  GENBUN */
/* ***ROUTINES CALLED  COSGEN, S1MERG, TRIX */
/* ***REVISION HISTORY  (YYMMDD) */
/*   801001  DATE WRITTEN */
/*   890531  Changed all specific intrinsics to generic.  (WRB) */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/*   900402  Added TYPE section.  (WRB) */
/*   920130  Modified to use merge routine S1MERG rather than deleted */
/*           routine MERGE.  (WRB) */
/* ***END PROLOGUE  POISD2 */

/* ***FIRST EXECUTABLE STATEMENT  POISD2 */
    /* Parameter adjustments */
    --ba;
    --bb;
    --bc;
    q_dim1 = *idimq;
    q_offset = 1 + q_dim1;
    q -= q_offset;
    --b;
    --w;
    --d__;
    --tcos;
    --p;

    /* Function Body */
    m = *mr;
    n = *nr;
    jsh = 0;
    fi = 1.f / *istag;
    ip = -m;
    ipstor = 0;
    switch (*istag) {
	case 1:  goto L101;
	case 2:  goto L102;
    }
L101:
    kr = 0;
    irreg = 1;
    if (n > 1) {
	goto L106;
    }
    tcos[1] = 0.f;
    goto L103;
L102:
    kr = 1;
    jstsav = 1;
    irreg = 2;
    if (n > 1) {
	goto L106;
    }
    tcos[1] = -1.f;
L103:
    i__1 = m;
    for (i__ = 1; i__ <= i__1; ++i__) {
	b[i__] = q[i__ + q_dim1];
/* L104: */
    }
    trix_(&c__1, &c__0, &m, &ba[1], &bb[1], &bc[1], &b[1], &tcos[1], &d__[1], 
	    &w[1]);
    i__1 = m;
    for (i__ = 1; i__ <= i__1; ++i__) {
	q[i__ + q_dim1] = b[i__];
/* L105: */
    }
    goto L183;
L106:
    lr = 0;
    i__1 = m;
    for (i__ = 1; i__ <= i__1; ++i__) {
	p[i__] = 0.f;
/* L107: */
    }
    nun = n;
    jst = 1;
    jsp = n;

/*     IRREG = 1 WHEN NO IRREGULARITIES HAVE OCCURRED, OTHERWISE IT IS 2. */

L108:
    l = jst << 1;
    nodd = 2 - ((nun + 1) / 2 << 1) + nun;

/*     NODD = 1 WHEN NUN IS ODD, OTHERWISE IT IS 2. */

    switch (nodd) {
	case 1:  goto L110;
	case 2:  goto L109;
    }
L109:
    jsp -= l;
    goto L111;
L110:
    jsp -= jst;
    if (irreg != 1) {
	jsp -= l;
    }
L111:

/*     REGULAR REDUCTION */

    cosgen_(&jst, &c__1, &c_b17, &c_b18, &tcos[1]);
    if (l > jsp) {
	goto L118;
    }
    i__1 = jsp;
    i__2 = l;
    for (j = l; i__2 < 0 ? j >= i__1 : j <= i__1; j += i__2) {
	jm1 = j - jsh;
	jp1 = j + jsh;
	jm2 = j - jst;
	jp2 = j + jst;
	jm3 = jm2 - jsh;
	jp3 = jp2 + jsh;
	if (jst != 1) {
	    goto L113;
	}
	i__3 = m;
	for (i__ = 1; i__ <= i__3; ++i__) {
	    b[i__] = q[i__ + j * q_dim1] * 2.f;
	    q[i__ + j * q_dim1] = q[i__ + jm2 * q_dim1] + q[i__ + jp2 * 
		    q_dim1];
/* L112: */
	}
	goto L115;
L113:
	i__3 = m;
	for (i__ = 1; i__ <= i__3; ++i__) {
	    t = q[i__ + j * q_dim1] - q[i__ + jm1 * q_dim1] - q[i__ + jp1 * 
		    q_dim1] + q[i__ + jm2 * q_dim1] + q[i__ + jp2 * q_dim1];
	    b[i__] = t + q[i__ + j * q_dim1] - q[i__ + jm3 * q_dim1] - q[i__ 
		    + jp3 * q_dim1];
	    q[i__ + j * q_dim1] = t;
/* L114: */
	}
L115:
	trix_(&jst, &c__0, &m, &ba[1], &bb[1], &bc[1], &b[1], &tcos[1], &d__[
		1], &w[1]);
	i__3 = m;
	for (i__ = 1; i__ <= i__3; ++i__) {
	    q[i__ + j * q_dim1] += b[i__];
/* L116: */
	}
/* L117: */
    }

/*     REDUCTION FOR LAST UNKNOWN */

L118:
    switch (nodd) {
	case 1:  goto L119;
	case 2:  goto L136;
    }
L119:
    switch (irreg) {
	case 1:  goto L152;
	case 2:  goto L120;
    }

/*     ODD NUMBER OF UNKNOWNS */

L120:
    jsp += l;
    j = jsp;
    jm1 = j - jsh;
    jp1 = j + jsh;
    jm2 = j - jst;
    jp2 = j + jst;
    jm3 = jm2 - jsh;
    switch (*istag) {
	case 1:  goto L123;
	case 2:  goto L121;
    }
L121:
    if (jst != 1) {
	goto L123;
    }
    i__2 = m;
    for (i__ = 1; i__ <= i__2; ++i__) {
	b[i__] = q[i__ + j * q_dim1];
	q[i__ + j * q_dim1] = 0.f;
/* L122: */
    }
    goto L130;
L123:
    switch (noddpr) {
	case 1:  goto L124;
	case 2:  goto L126;
    }
L124:
    i__2 = m;
    for (i__ = 1; i__ <= i__2; ++i__) {
	ip1 = ip + i__;
	b[i__] = (q[i__ + jm2 * q_dim1] - q[i__ + jm1 * q_dim1] - q[i__ + jm3 
		* q_dim1]) * .5f + p[ip1] + q[i__ + j * q_dim1];
/* L125: */
    }
    goto L128;
L126:
    i__2 = m;
    for (i__ = 1; i__ <= i__2; ++i__) {
	b[i__] = (q[i__ + jm2 * q_dim1] - q[i__ + jm1 * q_dim1] - q[i__ + jm3 
		* q_dim1]) * .5f + q[i__ + jp2 * q_dim1] - q[i__ + jp1 * 
		q_dim1] + q[i__ + j * q_dim1];
/* L127: */
    }
L128:
    i__2 = m;
    for (i__ = 1; i__ <= i__2; ++i__) {
	q[i__ + j * q_dim1] = (q[i__ + j * q_dim1] - q[i__ + jm1 * q_dim1] - 
		q[i__ + jp1 * q_dim1]) * .5f;
/* L129: */
    }
L130:
    trix_(&jst, &c__0, &m, &ba[1], &bb[1], &bc[1], &b[1], &tcos[1], &d__[1], &
	    w[1]);
    ip += m;
/* Computing MAX */
    i__2 = ipstor, i__1 = ip + m;
    ipstor = max(i__2,i__1);
    i__2 = m;
    for (i__ = 1; i__ <= i__2; ++i__) {
	ip1 = ip + i__;
	p[ip1] = q[i__ + j * q_dim1] + b[i__];
	b[i__] = q[i__ + jp2 * q_dim1] + p[ip1];
/* L131: */
    }
    if (lr != 0) {
	goto L133;
    }
    i__2 = jst;
    for (i__ = 1; i__ <= i__2; ++i__) {
	krpi = kr + i__;
	tcos[krpi] = tcos[i__];
/* L132: */
    }
    goto L134;
L133:
    cosgen_(&lr, &jstsav, &c_b18, &fi, &tcos[jst + 1]);
    s1merg_(&tcos[1], &c__0, &jst, &jst, &lr, &kr);
L134:
    cosgen_(&kr, &jstsav, &c_b18, &fi, &tcos[1]);
    trix_(&kr, &kr, &m, &ba[1], &bb[1], &bc[1], &b[1], &tcos[1], &d__[1], &w[
	    1]);
    i__2 = m;
    for (i__ = 1; i__ <= i__2; ++i__) {
	ip1 = ip + i__;
	q[i__ + j * q_dim1] = q[i__ + jm2 * q_dim1] + b[i__] + p[ip1];
/* L135: */
    }
    lr = kr;
    kr += l;
    goto L152;

/*     EVEN NUMBER OF UNKNOWNS */

L136:
    jsp += l;
    j = jsp;
    jm1 = j - jsh;
    jp1 = j + jsh;
    jm2 = j - jst;
    jp2 = j + jst;
    jm3 = jm2 - jsh;
    switch (irreg) {
	case 1:  goto L137;
	case 2:  goto L138;
    }
L137:
    jstsav = jst;
    ideg = jst;
    kr = l;
    goto L139;
L138:
    cosgen_(&kr, &jstsav, &c_b18, &fi, &tcos[1]);
    cosgen_(&lr, &jstsav, &c_b18, &fi, &tcos[kr + 1]);
    ideg = kr;
    kr += jst;
L139:
    if (jst != 1) {
	goto L141;
    }
    irreg = 2;
    i__2 = m;
    for (i__ = 1; i__ <= i__2; ++i__) {
	b[i__] = q[i__ + j * q_dim1];
	q[i__ + j * q_dim1] = q[i__ + jm2 * q_dim1];
/* L140: */
    }
    goto L150;
L141:
    i__2 = m;
    for (i__ = 1; i__ <= i__2; ++i__) {
	b[i__] = q[i__ + j * q_dim1] + (q[i__ + jm2 * q_dim1] - q[i__ + jm1 * 
		q_dim1] - q[i__ + jm3 * q_dim1]) * .5f;
/* L142: */
    }
    switch (irreg) {
	case 1:  goto L143;
	case 2:  goto L145;
    }
L143:
    i__2 = m;
    for (i__ = 1; i__ <= i__2; ++i__) {
	q[i__ + j * q_dim1] = q[i__ + jm2 * q_dim1] + (q[i__ + j * q_dim1] - 
		q[i__ + jm1 * q_dim1] - q[i__ + jp1 * q_dim1]) * .5f;
/* L144: */
    }
    irreg = 2;
    goto L150;
L145:
    switch (noddpr) {
	case 1:  goto L146;
	case 2:  goto L148;
    }
L146:
    i__2 = m;
    for (i__ = 1; i__ <= i__2; ++i__) {
	ip1 = ip + i__;
	q[i__ + j * q_dim1] = q[i__ + jm2 * q_dim1] + p[ip1];
/* L147: */
    }
    ip -= m;
    goto L150;
L148:
    i__2 = m;
    for (i__ = 1; i__ <= i__2; ++i__) {
	q[i__ + j * q_dim1] = q[i__ + jm2 * q_dim1] + q[i__ + j * q_dim1] - q[
		i__ + jm1 * q_dim1];
/* L149: */
    }
L150:
    trix_(&ideg, &lr, &m, &ba[1], &bb[1], &bc[1], &b[1], &tcos[1], &d__[1], &
	    w[1]);
    i__2 = m;
    for (i__ = 1; i__ <= i__2; ++i__) {
	q[i__ + j * q_dim1] += b[i__];
/* L151: */
    }
L152:
    nun /= 2;
    noddpr = nodd;
    jsh = jst;
    jst <<= 1;
    if (nun >= 2) {
	goto L108;
    }

/*     START SOLUTION. */

    j = jsp;
    i__2 = m;
    for (i__ = 1; i__ <= i__2; ++i__) {
	b[i__] = q[i__ + j * q_dim1];
/* L153: */
    }
    switch (irreg) {
	case 1:  goto L154;
	case 2:  goto L155;
    }
L154:
    cosgen_(&jst, &c__1, &c_b17, &c_b18, &tcos[1]);
    ideg = jst;
    goto L156;
L155:
    kr = lr + jst;
    cosgen_(&kr, &jstsav, &c_b18, &fi, &tcos[1]);
    cosgen_(&lr, &jstsav, &c_b18, &fi, &tcos[kr + 1]);
    ideg = kr;
L156:
    trix_(&ideg, &lr, &m, &ba[1], &bb[1], &bc[1], &b[1], &tcos[1], &d__[1], &
	    w[1]);
    jm1 = j - jsh;
    jp1 = j + jsh;
    switch (irreg) {
	case 1:  goto L157;
	case 2:  goto L159;
    }
L157:
    i__2 = m;
    for (i__ = 1; i__ <= i__2; ++i__) {
	q[i__ + j * q_dim1] = (q[i__ + j * q_dim1] - q[i__ + jm1 * q_dim1] - 
		q[i__ + jp1 * q_dim1]) * .5f + b[i__];
/* L158: */
    }
    goto L164;
L159:
    switch (noddpr) {
	case 1:  goto L160;
	case 2:  goto L162;
    }
L160:
    i__2 = m;
    for (i__ = 1; i__ <= i__2; ++i__) {
	ip1 = ip + i__;
	q[i__ + j * q_dim1] = p[ip1] + b[i__];
/* L161: */
    }
    ip -= m;
    goto L164;
L162:
    i__2 = m;
    for (i__ = 1; i__ <= i__2; ++i__) {
	q[i__ + j * q_dim1] = q[i__ + j * q_dim1] - q[i__ + jm1 * q_dim1] + b[
		i__];
/* L163: */
    }
L164:

/*     START BACK SUBSTITUTION. */

    jst /= 2;
    jsh = jst / 2;
    nun <<= 1;
    if (nun > n) {
	goto L183;
    }
    i__2 = n;
    i__1 = l;
    for (j = jst; i__1 < 0 ? j >= i__2 : j <= i__2; j += i__1) {
	jm1 = j - jsh;
	jp1 = j + jsh;
	jm2 = j - jst;
	jp2 = j + jst;
	if (j > jst) {
	    goto L166;
	}
	i__3 = m;
	for (i__ = 1; i__ <= i__3; ++i__) {
	    b[i__] = q[i__ + j * q_dim1] + q[i__ + jp2 * q_dim1];
/* L165: */
	}
	goto L170;
L166:
	if (jp2 <= n) {
	    goto L168;
	}
	i__3 = m;
	for (i__ = 1; i__ <= i__3; ++i__) {
	    b[i__] = q[i__ + j * q_dim1] + q[i__ + jm2 * q_dim1];
/* L167: */
	}
	if (jst < jstsav) {
	    irreg = 1;
	}
	switch (irreg) {
	    case 1:  goto L170;
	    case 2:  goto L171;
	}
L168:
	i__3 = m;
	for (i__ = 1; i__ <= i__3; ++i__) {
	    b[i__] = q[i__ + j * q_dim1] + q[i__ + jm2 * q_dim1] + q[i__ + 
		    jp2 * q_dim1];
/* L169: */
	}
L170:
	cosgen_(&jst, &c__1, &c_b17, &c_b18, &tcos[1]);
	ideg = jst;
	jdeg = 0;
	goto L172;
L171:
	if (j + l > n) {
	    lr -= jst;
	}
	kr = jst + lr;
	cosgen_(&kr, &jstsav, &c_b18, &fi, &tcos[1]);
	cosgen_(&lr, &jstsav, &c_b18, &fi, &tcos[kr + 1]);
	ideg = kr;
	jdeg = lr;
L172:
	trix_(&ideg, &jdeg, &m, &ba[1], &bb[1], &bc[1], &b[1], &tcos[1], &d__[
		1], &w[1]);
	if (jst > 1) {
	    goto L174;
	}
	i__3 = m;
	for (i__ = 1; i__ <= i__3; ++i__) {
	    q[i__ + j * q_dim1] = b[i__];
/* L173: */
	}
	goto L182;
L174:
	if (jp2 > n) {
	    goto L177;
	}
L175:
	i__3 = m;
	for (i__ = 1; i__ <= i__3; ++i__) {
	    q[i__ + j * q_dim1] = (q[i__ + j * q_dim1] - q[i__ + jm1 * q_dim1]
		     - q[i__ + jp1 * q_dim1]) * .5f + b[i__];
/* L176: */
	}
	goto L182;
L177:
	switch (irreg) {
	    case 1:  goto L175;
	    case 2:  goto L178;
	}
L178:
	if (j + jsh > n) {
	    goto L180;
	}
	i__3 = m;
	for (i__ = 1; i__ <= i__3; ++i__) {
	    ip1 = ip + i__;
	    q[i__ + j * q_dim1] = b[i__] + p[ip1];
/* L179: */
	}
	ip -= m;
	goto L182;
L180:
	i__3 = m;
	for (i__ = 1; i__ <= i__3; ++i__) {
	    q[i__ + j * q_dim1] = b[i__] + q[i__ + j * q_dim1] - q[i__ + jm1 *
		     q_dim1];
/* L181: */
	}
L182:
	;
    }
    l /= 2;
    goto L164;
L183:

/*     RETURN STORAGE REQUIREMENTS FOR P VECTORS. */

    w[1] = (real) ipstor;
    return 0;
} /* poisd2_ */

/* poisn2.f -- translated by f2c (version 20031025).
   You must link the resulting object file with libf2c:
	on Microsoft Windows system, link with libf2c.lib;
	on Linux or Unix systems, link with .../path/to/libf2c.a -lm
	or, if you install libf2c.a in a standard place, with -lf2c -lm
	-- in that order, at the end of the command line, as in
		cc *.o -lf2c -lm
	Source for libf2c is in /netlib/f2c/libf2c.zip, e.g.,

		http://www.netlib.org/f2c/libf2c.zip
*/

#include "f2c.h"

/* Table of constant values */

static real c_b11 = .5f;
static real c_b12 = 0.f;
static real c_b120 = 1.f;

/* DECK POISN2 */
/* Subroutine */ int poisn2_(integer *m, integer *n, integer *istag, integer *
	mixbnd, real *a, real *bb, real *c__, real *q, integer *idimq, real *
	b, real *b2, real *b3, real *w, real *w2, real *w3, real *d__, real *
	tcos, real *p)
{
    /* System generated locals */
    integer q_dim1, q_offset, i__1, i__2, i__3;
    static integer equiv_3[4];

    /* Local variables */
    static integer i__, j;
#define k (equiv_3)
    static real t;
    static integer i1, i2;
#define k1 (equiv_3)
#define k2 (equiv_3 + 1)
#define k3 (equiv_3 + 2)
#define k4 (equiv_3 + 3)
    static real fi;
    static integer ii, ip, jr, kr, lr, mr, nr, jm1, jm2, jm3, jp1, jp2, i2r, 
	    jp3, jr2;
    extern /* Subroutine */ int tri3_(integer *, real *, real *, real *, 
	    integer *, real *, real *, real *, real *, real *, real *, real *,
	     real *);
    static real fden;
    static integer nrod;
    static real fnum;
    extern /* Subroutine */ int trix_(integer *, integer *, integer *, real *,
	     real *, real *, real *, real *, real *, real *);
    static integer i2rby2, nlast, jstep, jstop;
    extern /* Subroutine */ int s1merg_(real *, integer *, integer *, integer 
	    *, integer *, integer *);
    static real fistag;
    extern /* Subroutine */ int cosgen_(integer *, integer *, real *, real *, 
	    real *);
    static integer nlastp, nrodpr, jstart, ipstor;

/* ***BEGIN PROLOGUE  POISN2 */
/* ***SUBSIDIARY */
/* ***PURPOSE  Subsidiary to GENBUN */
/* ***LIBRARY   SLATEC */
/* ***TYPE      SINGLE PRECISION (POISN2-S, CMPOSN-C) */
/* ***AUTHOR  (UNKNOWN) */
/* ***DESCRIPTION */

/*     Subroutine to solve Poisson's equation with Neumann boundary */
/*     conditions. */

/*     ISTAG = 1 if the last diagonal block is A. */
/*     ISTAG = 2 if the last diagonal block is A-I. */
/*     MIXBND = 1 if have Neumann boundary conditions at both boundaries. */
/*     MIXBND = 2 if have Neumann boundary conditions at bottom and */
/*     Dirichlet condition at top.  (for this case, must have ISTAG = 1.) */

/* ***SEE ALSO  GENBUN */
/* ***ROUTINES CALLED  COSGEN, S1MERG, TRI3, TRIX */
/* ***REVISION HISTORY  (YYMMDD) */
/*   801001  DATE WRITTEN */
/*   890531  Changed all specific intrinsics to generic.  (WRB) */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/*   900402  Added TYPE section.  (WRB) */
/*   920130  Modified to use merge routine S1MERG rather than deleted */
/*           routine MERGE.  (WRB) */
/* ***END PROLOGUE  POISN2 */

/* ***FIRST EXECUTABLE STATEMENT  POISN2 */
    /* Parameter adjustments */
    --a;
    --bb;
    --c__;
    q_dim1 = *idimq;
    q_offset = 1 + q_dim1;
    q -= q_offset;
    --b;
    --b2;
    --b3;
    --w;
    --w2;
    --w3;
    --d__;
    --tcos;
    --p;

    /* Function Body */
    fistag = (real) (3 - *istag);
    fnum = 1.f / *istag;
    fden = (*istag - 1) * .5f;
    mr = *m;
    ip = -mr;
    ipstor = 0;
    i2r = 1;
    jr = 2;
    nr = *n;
    nlast = *n;
    kr = 1;
    lr = 0;
    switch (*istag) {
	case 1:  goto L101;
	case 2:  goto L103;
    }
L101:
    i__1 = mr;
    for (i__ = 1; i__ <= i__1; ++i__) {
	q[i__ + *n * q_dim1] *= .5f;
/* L102: */
    }
    switch (*mixbnd) {
	case 1:  goto L103;
	case 2:  goto L104;
    }
L103:
    if (*n <= 3) {
	goto L155;
    }
L104:
    jr = i2r << 1;
    nrod = 1;
    if (nr / 2 << 1 == nr) {
	nrod = 0;
    }
    switch (*mixbnd) {
	case 1:  goto L105;
	case 2:  goto L106;
    }
L105:
    jstart = 1;
    goto L107;
L106:
    jstart = jr;
    nrod = 1 - nrod;
L107:
    jstop = nlast - jr;
    if (nrod == 0) {
	jstop -= i2r;
    }
    cosgen_(&i2r, &c__1, &c_b11, &c_b12, &tcos[1]);
    i2rby2 = i2r / 2;
    if (jstop >= jstart) {
	goto L108;
    }
    j = jr;
    goto L116;
L108:

/*     REGULAR REDUCTION. */

    i__1 = jstop;
    i__2 = jr;
    for (j = jstart; i__2 < 0 ? j >= i__1 : j <= i__1; j += i__2) {
	jp1 = j + i2rby2;
	jp2 = j + i2r;
	jp3 = jp2 + i2rby2;
	jm1 = j - i2rby2;
	jm2 = j - i2r;
	jm3 = jm2 - i2rby2;
	if (j != 1) {
	    goto L109;
	}
	jm1 = jp1;
	jm2 = jp2;
	jm3 = jp3;
L109:
	if (i2r != 1) {
	    goto L111;
	}
	if (j == 1) {
	    jm2 = jp2;
	}
	i__3 = mr;
	for (i__ = 1; i__ <= i__3; ++i__) {
	    b[i__] = q[i__ + j * q_dim1] * 2.f;
	    q[i__ + j * q_dim1] = q[i__ + jm2 * q_dim1] + q[i__ + jp2 * 
		    q_dim1];
/* L110: */
	}
	goto L113;
L111:
	i__3 = mr;
	for (i__ = 1; i__ <= i__3; ++i__) {
	    fi = q[i__ + j * q_dim1];
	    q[i__ + j * q_dim1] = q[i__ + j * q_dim1] - q[i__ + jm1 * q_dim1] 
		    - q[i__ + jp1 * q_dim1] + q[i__ + jm2 * q_dim1] + q[i__ + 
		    jp2 * q_dim1];
	    b[i__] = fi + q[i__ + j * q_dim1] - q[i__ + jm3 * q_dim1] - q[i__ 
		    + jp3 * q_dim1];
/* L112: */
	}
L113:
	trix_(&i2r, &c__0, &mr, &a[1], &bb[1], &c__[1], &b[1], &tcos[1], &d__[
		1], &w[1]);
	i__3 = mr;
	for (i__ = 1; i__ <= i__3; ++i__) {
	    q[i__ + j * q_dim1] += b[i__];
/* L114: */
	}

/*     END OF REDUCTION FOR REGULAR UNKNOWNS. */

/* L115: */
    }

/*     BEGIN SPECIAL REDUCTION FOR LAST UNKNOWN. */

    j = jstop + jr;
L116:
    nlast = j;
    jm1 = j - i2rby2;
    jm2 = j - i2r;
    jm3 = jm2 - i2rby2;
    if (nrod == 0) {
	goto L128;
    }

/*     ODD NUMBER OF UNKNOWNS */

    if (i2r != 1) {
	goto L118;
    }
    i__2 = mr;
    for (i__ = 1; i__ <= i__2; ++i__) {
	b[i__] = fistag * q[i__ + j * q_dim1];
	q[i__ + j * q_dim1] = q[i__ + jm2 * q_dim1];
/* L117: */
    }
    goto L126;
L118:
    i__2 = mr;
    for (i__ = 1; i__ <= i__2; ++i__) {
	b[i__] = q[i__ + j * q_dim1] + (q[i__ + jm2 * q_dim1] - q[i__ + jm1 * 
		q_dim1] - q[i__ + jm3 * q_dim1]) * .5f;
/* L119: */
    }
    if (nrodpr != 0) {
	goto L121;
    }
    i__2 = mr;
    for (i__ = 1; i__ <= i__2; ++i__) {
	ii = ip + i__;
	q[i__ + j * q_dim1] = q[i__ + jm2 * q_dim1] + p[ii];
/* L120: */
    }
    ip -= mr;
    goto L123;
L121:
    i__2 = mr;
    for (i__ = 1; i__ <= i__2; ++i__) {
	q[i__ + j * q_dim1] = q[i__ + j * q_dim1] - q[i__ + jm1 * q_dim1] + q[
		i__ + jm2 * q_dim1];
/* L122: */
    }
L123:
    if (lr == 0) {
	goto L124;
    }
    cosgen_(&lr, &c__1, &c_b11, &fden, &tcos[kr + 1]);
    goto L126;
L124:
    i__2 = mr;
    for (i__ = 1; i__ <= i__2; ++i__) {
	b[i__] = fistag * b[i__];
/* L125: */
    }
L126:
    cosgen_(&kr, &c__1, &c_b11, &fden, &tcos[1]);
    trix_(&kr, &lr, &mr, &a[1], &bb[1], &c__[1], &b[1], &tcos[1], &d__[1], &w[
	    1]);
    i__2 = mr;
    for (i__ = 1; i__ <= i__2; ++i__) {
	q[i__ + j * q_dim1] += b[i__];
/* L127: */
    }
    kr += i2r;
    goto L151;
L128:

/*     EVEN NUMBER OF UNKNOWNS */

    jp1 = j + i2rby2;
    jp2 = j + i2r;
    if (i2r != 1) {
	goto L135;
    }
    i__2 = mr;
    for (i__ = 1; i__ <= i__2; ++i__) {
	b[i__] = q[i__ + j * q_dim1];
/* L129: */
    }
    trix_(&c__1, &c__0, &mr, &a[1], &bb[1], &c__[1], &b[1], &tcos[1], &d__[1],
	     &w[1]);
    ip = 0;
    ipstor = mr;
    switch (*istag) {
	case 1:  goto L133;
	case 2:  goto L130;
    }
L130:
    i__2 = mr;
    for (i__ = 1; i__ <= i__2; ++i__) {
	p[i__] = b[i__];
	b[i__] += q[i__ + *n * q_dim1];
/* L131: */
    }
    tcos[1] = 1.f;
    tcos[2] = 0.f;
    trix_(&c__1, &c__1, &mr, &a[1], &bb[1], &c__[1], &b[1], &tcos[1], &d__[1],
	     &w[1]);
    i__2 = mr;
    for (i__ = 1; i__ <= i__2; ++i__) {
	q[i__ + j * q_dim1] = q[i__ + jm2 * q_dim1] + p[i__] + b[i__];
/* L132: */
    }
    goto L150;
L133:
    i__2 = mr;
    for (i__ = 1; i__ <= i__2; ++i__) {
	p[i__] = b[i__];
	q[i__ + j * q_dim1] = q[i__ + jm2 * q_dim1] + q[i__ + jp2 * q_dim1] * 
		2.f + b[i__] * 3.f;
/* L134: */
    }
    goto L150;
L135:
    i__2 = mr;
    for (i__ = 1; i__ <= i__2; ++i__) {
	b[i__] = q[i__ + j * q_dim1] + (q[i__ + jm2 * q_dim1] - q[i__ + jm1 * 
		q_dim1] - q[i__ + jm3 * q_dim1]) * .5f;
/* L136: */
    }
    if (nrodpr != 0) {
	goto L138;
    }
    i__2 = mr;
    for (i__ = 1; i__ <= i__2; ++i__) {
	ii = ip + i__;
	b[i__] += p[ii];
/* L137: */
    }
    goto L140;
L138:
    i__2 = mr;
    for (i__ = 1; i__ <= i__2; ++i__) {
	b[i__] = b[i__] + q[i__ + jp2 * q_dim1] - q[i__ + jp1 * q_dim1];
/* L139: */
    }
L140:
    trix_(&i2r, &c__0, &mr, &a[1], &bb[1], &c__[1], &b[1], &tcos[1], &d__[1], 
	    &w[1]);
    ip += mr;
/* Computing MAX */
    i__2 = ipstor, i__1 = ip + mr;
    ipstor = max(i__2,i__1);
    i__2 = mr;
    for (i__ = 1; i__ <= i__2; ++i__) {
	ii = ip + i__;
	p[ii] = b[i__] + (q[i__ + j * q_dim1] - q[i__ + jm1 * q_dim1] - q[i__ 
		+ jp1 * q_dim1]) * .5f;
	b[i__] = p[ii] + q[i__ + jp2 * q_dim1];
/* L141: */
    }
    if (lr == 0) {
	goto L142;
    }
    cosgen_(&lr, &c__1, &c_b11, &fden, &tcos[i2r + 1]);
    s1merg_(&tcos[1], &c__0, &i2r, &i2r, &lr, &kr);
    goto L144;
L142:
    i__2 = i2r;
    for (i__ = 1; i__ <= i__2; ++i__) {
	ii = kr + i__;
	tcos[ii] = tcos[i__];
/* L143: */
    }
L144:
    cosgen_(&kr, &c__1, &c_b11, &fden, &tcos[1]);
    if (lr != 0) {
	goto L145;
    }
    switch (*istag) {
	case 1:  goto L146;
	case 2:  goto L145;
    }
L145:
    trix_(&kr, &kr, &mr, &a[1], &bb[1], &c__[1], &b[1], &tcos[1], &d__[1], &w[
	    1]);
    goto L148;
L146:
    i__2 = mr;
    for (i__ = 1; i__ <= i__2; ++i__) {
	b[i__] = fistag * b[i__];
/* L147: */
    }
L148:
    i__2 = mr;
    for (i__ = 1; i__ <= i__2; ++i__) {
	ii = ip + i__;
	q[i__ + j * q_dim1] = q[i__ + jm2 * q_dim1] + p[ii] + b[i__];
/* L149: */
    }
L150:
    lr = kr;
    kr += jr;
L151:
    switch (*mixbnd) {
	case 1:  goto L152;
	case 2:  goto L153;
    }
L152:
    nr = (nlast - 1) / jr + 1;
    if (nr <= 3) {
	goto L155;
    }
    goto L154;
L153:
    nr = nlast / jr;
    if (nr <= 1) {
	goto L192;
    }
L154:
    i2r = jr;
    nrodpr = nrod;
    goto L104;
L155:

/*      BEGIN SOLUTION */

    j = jr + 1;
    jm1 = j - i2r;
    jp1 = j + i2r;
    jm2 = nlast - i2r;
    if (nr == 2) {
	goto L184;
    }
    if (lr != 0) {
	goto L170;
    }
    if (*n != 3) {
	goto L161;
    }

/*     CASE N = 3. */

    switch (*istag) {
	case 1:  goto L156;
	case 2:  goto L168;
    }
L156:
    i__2 = mr;
    for (i__ = 1; i__ <= i__2; ++i__) {
	b[i__] = q[i__ + (q_dim1 << 1)];
/* L157: */
    }
    tcos[1] = 0.f;
    trix_(&c__1, &c__0, &mr, &a[1], &bb[1], &c__[1], &b[1], &tcos[1], &d__[1],
	     &w[1]);
    i__2 = mr;
    for (i__ = 1; i__ <= i__2; ++i__) {
	q[i__ + (q_dim1 << 1)] = b[i__];
	b[i__] = b[i__] * 4.f + q[i__ + q_dim1] + q[i__ + q_dim1 * 3] * 2.f;
/* L158: */
    }
    tcos[1] = -2.f;
    tcos[2] = 2.f;
    i1 = 2;
    i2 = 0;
    trix_(&i1, &i2, &mr, &a[1], &bb[1], &c__[1], &b[1], &tcos[1], &d__[1], &w[
	    1]);
    i__2 = mr;
    for (i__ = 1; i__ <= i__2; ++i__) {
	q[i__ + (q_dim1 << 1)] += b[i__];
	b[i__] = q[i__ + q_dim1] + q[i__ + (q_dim1 << 1)] * 2.f;
/* L159: */
    }
    tcos[1] = 0.f;
    trix_(&c__1, &c__0, &mr, &a[1], &bb[1], &c__[1], &b[1], &tcos[1], &d__[1],
	     &w[1]);
    i__2 = mr;
    for (i__ = 1; i__ <= i__2; ++i__) {
	q[i__ + q_dim1] = b[i__];
/* L160: */
    }
    jr = 1;
    i2r = 0;
    goto L194;

/*     CASE N = 2**P+1 */

L161:
    switch (*istag) {
	case 1:  goto L162;
	case 2:  goto L170;
    }
L162:
    i__2 = mr;
    for (i__ = 1; i__ <= i__2; ++i__) {
	b[i__] = q[i__ + j * q_dim1] + q[i__ + q_dim1] * .5f - q[i__ + jm1 * 
		q_dim1] + q[i__ + nlast * q_dim1] - q[i__ + jm2 * q_dim1];
/* L163: */
    }
    cosgen_(&jr, &c__1, &c_b11, &c_b12, &tcos[1]);
    trix_(&jr, &c__0, &mr, &a[1], &bb[1], &c__[1], &b[1], &tcos[1], &d__[1], &
	    w[1]);
    i__2 = mr;
    for (i__ = 1; i__ <= i__2; ++i__) {
	q[i__ + j * q_dim1] = (q[i__ + j * q_dim1] - q[i__ + jm1 * q_dim1] - 
		q[i__ + jp1 * q_dim1]) * .5f + b[i__];
	b[i__] = q[i__ + q_dim1] + q[i__ + nlast * q_dim1] * 2.f + q[i__ + j *
		 q_dim1] * 4.f;
/* L164: */
    }
    jr2 = jr << 1;
    cosgen_(&jr, &c__1, &c_b12, &c_b12, &tcos[1]);
    i__2 = jr;
    for (i__ = 1; i__ <= i__2; ++i__) {
	i1 = jr + i__;
	i2 = jr + 1 - i__;
	tcos[i1] = -tcos[i2];
/* L165: */
    }
    trix_(&jr2, &c__0, &mr, &a[1], &bb[1], &c__[1], &b[1], &tcos[1], &d__[1], 
	    &w[1]);
    i__2 = mr;
    for (i__ = 1; i__ <= i__2; ++i__) {
	q[i__ + j * q_dim1] += b[i__];
	b[i__] = q[i__ + q_dim1] + q[i__ + j * q_dim1] * 2.f;
/* L166: */
    }
    cosgen_(&jr, &c__1, &c_b11, &c_b12, &tcos[1]);
    trix_(&jr, &c__0, &mr, &a[1], &bb[1], &c__[1], &b[1], &tcos[1], &d__[1], &
	    w[1]);
    i__2 = mr;
    for (i__ = 1; i__ <= i__2; ++i__) {
	q[i__ + q_dim1] = q[i__ + q_dim1] * .5f - q[i__ + jm1 * q_dim1] + b[
		i__];
/* L167: */
    }
    goto L194;

/*     CASE OF GENERAL N WITH NR = 3 . */

L168:
    i__2 = mr;
    for (i__ = 1; i__ <= i__2; ++i__) {
	b[i__] = q[i__ + (q_dim1 << 1)];
	q[i__ + (q_dim1 << 1)] = 0.f;
	b2[i__] = q[i__ + q_dim1 * 3];
	b3[i__] = q[i__ + q_dim1];
/* L169: */
    }
    jr = 1;
    i2r = 0;
    j = 2;
    goto L177;
L170:
    i__2 = mr;
    for (i__ = 1; i__ <= i__2; ++i__) {
	b[i__] = q[i__ + q_dim1] * .5f - q[i__ + jm1 * q_dim1] + q[i__ + j * 
		q_dim1];
/* L171: */
    }
    if (nrod != 0) {
	goto L173;
    }
    i__2 = mr;
    for (i__ = 1; i__ <= i__2; ++i__) {
	ii = ip + i__;
	b[i__] += p[ii];
/* L172: */
    }
    goto L175;
L173:
    i__2 = mr;
    for (i__ = 1; i__ <= i__2; ++i__) {
	b[i__] = b[i__] + q[i__ + nlast * q_dim1] - q[i__ + jm2 * q_dim1];
/* L174: */
    }
L175:
    i__2 = mr;
    for (i__ = 1; i__ <= i__2; ++i__) {
	t = (q[i__ + j * q_dim1] - q[i__ + jm1 * q_dim1] - q[i__ + jp1 * 
		q_dim1]) * .5f;
	q[i__ + j * q_dim1] = t;
	b2[i__] = q[i__ + nlast * q_dim1] + t;
	b3[i__] = q[i__ + q_dim1] + t * 2.f;
/* L176: */
    }
L177:
    *k1 = kr + (jr << 1) - 1;
    *k2 = kr + jr;
    tcos[*k1 + 1] = -2.f;
    *k4 = *k1 + 3 - *istag;
    i__2 = *k2 + *istag - 2;
    cosgen_(&i__2, &c__1, &c_b12, &fnum, &tcos[*k4]);
    *k4 = *k1 + *k2 + 1;
    i__2 = jr - 1;
    cosgen_(&i__2, &c__1, &c_b12, &c_b120, &tcos[*k4]);
    i__2 = *k1 + *k2;
    i__1 = jr - 1;
    s1merg_(&tcos[1], k1, k2, &i__2, &i__1, &c__0);
    *k3 = *k1 + *k2 + lr;
    cosgen_(&jr, &c__1, &c_b11, &c_b12, &tcos[*k3 + 1]);
    *k4 = *k3 + jr + 1;
    cosgen_(&kr, &c__1, &c_b11, &fden, &tcos[*k4]);
    i__2 = *k3 + jr;
    s1merg_(&tcos[1], k3, &jr, &i__2, &kr, k1);
    if (lr == 0) {
	goto L178;
    }
    cosgen_(&lr, &c__1, &c_b11, &fden, &tcos[*k4]);
    i__2 = *k3 + jr;
    i__1 = *k3 - lr;
    s1merg_(&tcos[1], k3, &jr, &i__2, &lr, &i__1);
    cosgen_(&kr, &c__1, &c_b11, &fden, &tcos[*k4]);
L178:
    *k3 = kr;
    *k4 = kr;
    tri3_(&mr, &a[1], &bb[1], &c__[1], k, &b[1], &b2[1], &b3[1], &tcos[1], &
	    d__[1], &w[1], &w2[1], &w3[1]);
    i__2 = mr;
    for (i__ = 1; i__ <= i__2; ++i__) {
	b[i__] = b[i__] + b2[i__] + b3[i__];
/* L179: */
    }
    tcos[1] = 2.f;
    trix_(&c__1, &c__0, &mr, &a[1], &bb[1], &c__[1], &b[1], &tcos[1], &d__[1],
	     &w[1]);
    i__2 = mr;
    for (i__ = 1; i__ <= i__2; ++i__) {
	q[i__ + j * q_dim1] += b[i__];
	b[i__] = q[i__ + q_dim1] + q[i__ + j * q_dim1] * 2.f;
/* L180: */
    }
    cosgen_(&jr, &c__1, &c_b11, &c_b12, &tcos[1]);
    trix_(&jr, &c__0, &mr, &a[1], &bb[1], &c__[1], &b[1], &tcos[1], &d__[1], &
	    w[1]);
    if (jr != 1) {
	goto L182;
    }
    i__2 = mr;
    for (i__ = 1; i__ <= i__2; ++i__) {
	q[i__ + q_dim1] = b[i__];
/* L181: */
    }
    goto L194;
L182:
    i__2 = mr;
    for (i__ = 1; i__ <= i__2; ++i__) {
	q[i__ + q_dim1] = q[i__ + q_dim1] * .5f - q[i__ + jm1 * q_dim1] + b[
		i__];
/* L183: */
    }
    goto L194;
L184:
    if (*n != 2) {
	goto L188;
    }

/*     CASE  N = 2 */

    i__2 = mr;
    for (i__ = 1; i__ <= i__2; ++i__) {
	b[i__] = q[i__ + q_dim1];
/* L185: */
    }
    tcos[1] = 0.f;
    trix_(&c__1, &c__0, &mr, &a[1], &bb[1], &c__[1], &b[1], &tcos[1], &d__[1],
	     &w[1]);
    i__2 = mr;
    for (i__ = 1; i__ <= i__2; ++i__) {
	q[i__ + q_dim1] = b[i__];
	b[i__] = (q[i__ + (q_dim1 << 1)] + b[i__]) * 2.f * fistag;
/* L186: */
    }
    tcos[1] = -fistag;
    tcos[2] = 2.f;
    trix_(&c__2, &c__0, &mr, &a[1], &bb[1], &c__[1], &b[1], &tcos[1], &d__[1],
	     &w[1]);
    i__2 = mr;
    for (i__ = 1; i__ <= i__2; ++i__) {
	q[i__ + q_dim1] += b[i__];
/* L187: */
    }
    jr = 1;
    i2r = 0;
    goto L194;
L188:

/*     CASE OF GENERAL N AND NR = 2 . */

    i__2 = mr;
    for (i__ = 1; i__ <= i__2; ++i__) {
	ii = ip + i__;
	b3[i__] = 0.f;
	b[i__] = q[i__ + q_dim1] + p[ii] * 2.f;
	q[i__ + q_dim1] = q[i__ + q_dim1] * .5f - q[i__ + jm1 * q_dim1];
	b2[i__] = (q[i__ + q_dim1] + q[i__ + nlast * q_dim1]) * 2.f;
/* L189: */
    }
    *k1 = kr + jr - 1;
    tcos[*k1 + 1] = -2.f;
    *k4 = *k1 + 3 - *istag;
    i__2 = kr + *istag - 2;
    cosgen_(&i__2, &c__1, &c_b12, &fnum, &tcos[*k4]);
    *k4 = *k1 + kr + 1;
    i__2 = jr - 1;
    cosgen_(&i__2, &c__1, &c_b12, &c_b120, &tcos[*k4]);
    i__2 = *k1 + kr;
    i__1 = jr - 1;
    s1merg_(&tcos[1], k1, &kr, &i__2, &i__1, &c__0);
    cosgen_(&kr, &c__1, &c_b11, &fden, &tcos[*k1 + 1]);
    *k2 = kr;
    *k4 = *k1 + *k2 + 1;
    cosgen_(&lr, &c__1, &c_b11, &fden, &tcos[*k4]);
    *k3 = lr;
    *k4 = 0;
    tri3_(&mr, &a[1], &bb[1], &c__[1], k, &b[1], &b2[1], &b3[1], &tcos[1], &
	    d__[1], &w[1], &w2[1], &w3[1]);
    i__2 = mr;
    for (i__ = 1; i__ <= i__2; ++i__) {
	b[i__] += b2[i__];
/* L190: */
    }
    tcos[1] = 2.f;
    trix_(&c__1, &c__0, &mr, &a[1], &bb[1], &c__[1], &b[1], &tcos[1], &d__[1],
	     &w[1]);
    i__2 = mr;
    for (i__ = 1; i__ <= i__2; ++i__) {
	q[i__ + q_dim1] += b[i__];
/* L191: */
    }
    goto L194;
L192:
    i__2 = mr;
    for (i__ = 1; i__ <= i__2; ++i__) {
	b[i__] = q[i__ + nlast * q_dim1];
/* L193: */
    }
    goto L196;
L194:

/*     START BACK SUBSTITUTION. */

    j = nlast - jr;
    i__2 = mr;
    for (i__ = 1; i__ <= i__2; ++i__) {
	b[i__] = q[i__ + nlast * q_dim1] + q[i__ + j * q_dim1];
/* L195: */
    }
L196:
    jm2 = nlast - i2r;
    if (jr != 1) {
	goto L198;
    }
    i__2 = mr;
    for (i__ = 1; i__ <= i__2; ++i__) {
	q[i__ + nlast * q_dim1] = 0.f;
/* L197: */
    }
    goto L202;
L198:
    if (nrod != 0) {
	goto L200;
    }
    i__2 = mr;
    for (i__ = 1; i__ <= i__2; ++i__) {
	ii = ip + i__;
	q[i__ + nlast * q_dim1] = p[ii];
/* L199: */
    }
    ip -= mr;
    goto L202;
L200:
    i__2 = mr;
    for (i__ = 1; i__ <= i__2; ++i__) {
	q[i__ + nlast * q_dim1] -= q[i__ + jm2 * q_dim1];
/* L201: */
    }
L202:
    cosgen_(&kr, &c__1, &c_b11, &fden, &tcos[1]);
    cosgen_(&lr, &c__1, &c_b11, &fden, &tcos[kr + 1]);
    if (lr != 0) {
	goto L204;
    }
    i__2 = mr;
    for (i__ = 1; i__ <= i__2; ++i__) {
	b[i__] = fistag * b[i__];
/* L203: */
    }
L204:
    trix_(&kr, &lr, &mr, &a[1], &bb[1], &c__[1], &b[1], &tcos[1], &d__[1], &w[
	    1]);
    i__2 = mr;
    for (i__ = 1; i__ <= i__2; ++i__) {
	q[i__ + nlast * q_dim1] += b[i__];
/* L205: */
    }
    nlastp = nlast;
L206:
    jstep = jr;
    jr = i2r;
    i2r /= 2;
    if (jr == 0) {
	goto L222;
    }
    switch (*mixbnd) {
	case 1:  goto L207;
	case 2:  goto L208;
    }
L207:
    jstart = jr + 1;
    goto L209;
L208:
    jstart = jr;
L209:
    kr -= jr;
    if (nlast + jr > *n) {
	goto L210;
    }
    kr -= jr;
    nlast += jr;
    jstop = nlast - jstep;
    goto L211;
L210:
    jstop = nlast - jr;
L211:
    lr = kr - jr;
    cosgen_(&jr, &c__1, &c_b11, &c_b12, &tcos[1]);
    i__2 = jstop;
    i__1 = jstep;
    for (j = jstart; i__1 < 0 ? j >= i__2 : j <= i__2; j += i__1) {
	jm2 = j - jr;
	jp2 = j + jr;
	if (j != jr) {
	    goto L213;
	}
	i__3 = mr;
	for (i__ = 1; i__ <= i__3; ++i__) {
	    b[i__] = q[i__ + j * q_dim1] + q[i__ + jp2 * q_dim1];
/* L212: */
	}
	goto L215;
L213:
	i__3 = mr;
	for (i__ = 1; i__ <= i__3; ++i__) {
	    b[i__] = q[i__ + j * q_dim1] + q[i__ + jm2 * q_dim1] + q[i__ + 
		    jp2 * q_dim1];
/* L214: */
	}
L215:
	if (jr != 1) {
	    goto L217;
	}
	i__3 = mr;
	for (i__ = 1; i__ <= i__3; ++i__) {
	    q[i__ + j * q_dim1] = 0.f;
/* L216: */
	}
	goto L219;
L217:
	jm1 = j - i2r;
	jp1 = j + i2r;
	i__3 = mr;
	for (i__ = 1; i__ <= i__3; ++i__) {
	    q[i__ + j * q_dim1] = (q[i__ + j * q_dim1] - q[i__ + jm1 * q_dim1]
		     - q[i__ + jp1 * q_dim1]) * .5f;
/* L218: */
	}
L219:
	trix_(&jr, &c__0, &mr, &a[1], &bb[1], &c__[1], &b[1], &tcos[1], &d__[
		1], &w[1]);
	i__3 = mr;
	for (i__ = 1; i__ <= i__3; ++i__) {
	    q[i__ + j * q_dim1] += b[i__];
/* L220: */
	}
/* L221: */
    }
    nrod = 1;
    if (nlast + i2r <= *n) {
	nrod = 0;
    }
    if (nlastp != nlast) {
	goto L194;
    }
    goto L206;
L222:

/*     RETURN STORAGE REQUIREMENTS FOR P VECTORS. */

    w[1] = (real) ipstor;
    return 0;
} /* poisn2_ */

#undef k4
#undef k3
#undef k2
#undef k1
#undef k


/* poisp2.f -- translated by f2c (version 20031025).
   You must link the resulting object file with libf2c:
	on Microsoft Windows system, link with libf2c.lib;
	on Linux or Unix systems, link with .../path/to/libf2c.a -lm
	or, if you install libf2c.a in a standard place, with -lf2c -lm
	-- in that order, at the end of the command line, as in
		cc *.o -lf2c -lm
	Source for libf2c is in /netlib/f2c/libf2c.zip, e.g.,

		http://www.netlib.org/f2c/libf2c.zip
*/

#include "f2c.h"

/* Table of constant values */

/* DECK POISP2 */
/* Subroutine */ int poisp2_(integer *m, integer *n, real *a, real *bb, real *
	c__, real *q, integer *idimq, real *b, real *b2, real *b3, real *w, 
	real *w2, real *w3, real *d__, real *tcos, real *p)
{
    /* System generated locals */
    integer q_dim1, q_offset, i__1, i__2;

    /* Local variables */
    static integer i__, j;
    static real s, t;
    static integer lh, mr, nr, nrm1, nrmj, nrpj;
    extern /* Subroutine */ int poisd2_(integer *, integer *, integer *, real 
	    *, real *, real *, real *, integer *, real *, real *, real *, 
	    real *, real *), poisn2_(integer *, integer *, integer *, integer 
	    *, real *, real *, real *, real *, integer *, real *, real *, 
	    real *, real *, real *, real *, real *, real *, real *);
    static integer ipstor;

/* ***BEGIN PROLOGUE  POISP2 */
/* ***SUBSIDIARY */
/* ***PURPOSE  Subsidiary to GENBUN */
/* ***LIBRARY   SLATEC */
/* ***TYPE      SINGLE PRECISION (POISP2-S, CMPOSP-C) */
/* ***AUTHOR  (UNKNOWN) */
/* ***DESCRIPTION */

/*     Subroutine to solve Poisson equation with periodic boundary */
/*     conditions. */

/* ***SEE ALSO  GENBUN */
/* ***ROUTINES CALLED  POISD2, POISN2 */
/* ***REVISION HISTORY  (YYMMDD) */
/*   801001  DATE WRITTEN */
/*   890531  Changed all specific intrinsics to generic.  (WRB) */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/*   900402  Added TYPE section.  (WRB) */
/* ***END PROLOGUE  POISP2 */

/* ***FIRST EXECUTABLE STATEMENT  POISP2 */
    /* Parameter adjustments */
    --a;
    --bb;
    --c__;
    q_dim1 = *idimq;
    q_offset = 1 + q_dim1;
    q -= q_offset;
    --b;
    --b2;
    --b3;
    --w;
    --w2;
    --w3;
    --d__;
    --tcos;
    --p;

    /* Function Body */
    mr = *m;
    nr = (*n + 1) / 2;
    nrm1 = nr - 1;
    if (nr << 1 != *n) {
	goto L107;
    }

/*     EVEN NUMBER OF UNKNOWNS */

    i__1 = nrm1;
    for (j = 1; j <= i__1; ++j) {
	nrmj = nr - j;
	nrpj = nr + j;
	i__2 = mr;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    s = q[i__ + nrmj * q_dim1] - q[i__ + nrpj * q_dim1];
	    t = q[i__ + nrmj * q_dim1] + q[i__ + nrpj * q_dim1];
	    q[i__ + nrmj * q_dim1] = s;
	    q[i__ + nrpj * q_dim1] = t;
/* L101: */
	}
/* L102: */
    }
    i__1 = mr;
    for (i__ = 1; i__ <= i__1; ++i__) {
	q[i__ + nr * q_dim1] *= 2.f;
	q[i__ + *n * q_dim1] *= 2.f;
/* L103: */
    }
    poisd2_(&mr, &nrm1, &c__1, &a[1], &bb[1], &c__[1], &q[q_offset], idimq, &
	    b[1], &w[1], &d__[1], &tcos[1], &p[1]);
    ipstor = w[1];
    i__1 = nr + 1;
    poisn2_(&mr, &i__1, &c__1, &c__1, &a[1], &bb[1], &c__[1], &q[nr * q_dim1 
	    + 1], idimq, &b[1], &b2[1], &b3[1], &w[1], &w2[1], &w3[1], &d__[1]
	    , &tcos[1], &p[1]);
/* Computing MAX */
    i__1 = ipstor, i__2 = (integer) w[1];
    ipstor = max(i__1,i__2);
    i__1 = nrm1;
    for (j = 1; j <= i__1; ++j) {
	nrmj = nr - j;
	nrpj = nr + j;
	i__2 = mr;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    s = (q[i__ + nrpj * q_dim1] + q[i__ + nrmj * q_dim1]) * .5f;
	    t = (q[i__ + nrpj * q_dim1] - q[i__ + nrmj * q_dim1]) * .5f;
	    q[i__ + nrmj * q_dim1] = s;
	    q[i__ + nrpj * q_dim1] = t;
/* L104: */
	}
/* L105: */
    }
    i__1 = mr;
    for (i__ = 1; i__ <= i__1; ++i__) {
	q[i__ + nr * q_dim1] *= .5f;
	q[i__ + *n * q_dim1] *= .5f;
/* L106: */
    }
    goto L118;
L107:

/*     ODD  NUMBER OF UNKNOWNS */

    i__1 = nrm1;
    for (j = 1; j <= i__1; ++j) {
	nrpj = *n + 1 - j;
	i__2 = mr;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    s = q[i__ + j * q_dim1] - q[i__ + nrpj * q_dim1];
	    t = q[i__ + j * q_dim1] + q[i__ + nrpj * q_dim1];
	    q[i__ + j * q_dim1] = s;
	    q[i__ + nrpj * q_dim1] = t;
/* L108: */
	}
/* L109: */
    }
    i__1 = mr;
    for (i__ = 1; i__ <= i__1; ++i__) {
	q[i__ + nr * q_dim1] *= 2.f;
/* L110: */
    }
    lh = nrm1 / 2;
    i__1 = lh;
    for (j = 1; j <= i__1; ++j) {
	nrmj = nr - j;
	i__2 = mr;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    s = q[i__ + j * q_dim1];
	    q[i__ + j * q_dim1] = q[i__ + nrmj * q_dim1];
	    q[i__ + nrmj * q_dim1] = s;
/* L111: */
	}
/* L112: */
    }
    poisd2_(&mr, &nrm1, &c__2, &a[1], &bb[1], &c__[1], &q[q_offset], idimq, &
	    b[1], &w[1], &d__[1], &tcos[1], &p[1]);
    ipstor = w[1];
    poisn2_(&mr, &nr, &c__2, &c__1, &a[1], &bb[1], &c__[1], &q[nr * q_dim1 + 
	    1], idimq, &b[1], &b2[1], &b3[1], &w[1], &w2[1], &w3[1], &d__[1], 
	    &tcos[1], &p[1]);
/* Computing MAX */
    i__1 = ipstor, i__2 = (integer) w[1];
    ipstor = max(i__1,i__2);
    i__1 = nrm1;
    for (j = 1; j <= i__1; ++j) {
	nrpj = nr + j;
	i__2 = mr;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    s = (q[i__ + nrpj * q_dim1] + q[i__ + j * q_dim1]) * .5f;
	    t = (q[i__ + nrpj * q_dim1] - q[i__ + j * q_dim1]) * .5f;
	    q[i__ + nrpj * q_dim1] = t;
	    q[i__ + j * q_dim1] = s;
/* L113: */
	}
/* L114: */
    }
    i__1 = mr;
    for (i__ = 1; i__ <= i__1; ++i__) {
	q[i__ + nr * q_dim1] *= .5f;
/* L115: */
    }
    i__1 = lh;
    for (j = 1; j <= i__1; ++j) {
	nrmj = nr - j;
	i__2 = mr;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    s = q[i__ + j * q_dim1];
	    q[i__ + j * q_dim1] = q[i__ + nrmj * q_dim1];
	    q[i__ + nrmj * q_dim1] = s;
/* L116: */
	}
/* L117: */
    }
L118:

/*     RETURN STORAGE REQUIREMENTS FOR P VECTORS. */

    w[1] = (real) ipstor;
    return 0;
} /* poisp2_ */

/* cosgen.f -- translated by f2c (version 20031025).
   You must link the resulting object file with libf2c:
	on Microsoft Windows system, link with libf2c.lib;
	on Linux or Unix systems, link with .../path/to/libf2c.a -lm
	or, if you install libf2c.a in a standard place, with -lf2c -lm
	-- in that order, at the end of the command line, as in
		cc *.o -lf2c -lm
	Source for libf2c is in /netlib/f2c/libf2c.zip, e.g.,

		http://www.netlib.org/f2c/libf2c.zip
*/

#include "f2c.h"

/* DECK COSGEN */
/* Subroutine */ int cosgen_(integer *n, integer *ijump, real *fnum, real *
	fden, real *a)
{
    /* System generated locals */
    integer i__1, i__2;

    /* Builtin functions */
    double cos(doublereal);

    /* Local variables */
    static integer i__, k;
    static real x, y;
    static integer k1, k2, k3, k4, k5;
    static real pi;
    static integer np1;
    static real dum, pibyn;
    extern doublereal pimach_(real *);

/* ***BEGIN PROLOGUE  COSGEN */
/* ***SUBSIDIARY */
/* ***PURPOSE  Subsidiary to GENBUN */
/* ***LIBRARY   SLATEC */
/* ***TYPE      SINGLE PRECISION (COSGEN-S, CMPCSG-C) */
/* ***AUTHOR  (UNKNOWN) */
/* ***DESCRIPTION */

/*     This subroutine computes required cosine values in ascending */
/*     order.  When IJUMP .GT. 1 the routine computes values */

/*        2*COS(J*PI/L) , J=1,2,...,L and J .NE. 0(MOD N/IJUMP+1) */

/*     where L = IJUMP*(N/IJUMP+1). */


/*     when IJUMP = 1 it computes */

/*            2*COS((J-FNUM)*PI/(N+FDEN)) ,  J=1, 2, ... ,N */

/*     where */
/*        FNUM = 0.5, FDEN = 0.0, for regular reduction values. */
/*        FNUM = 0.0, FDEN = 1.0, for B-R and C-R when ISTAG = 1 */
/*        FNUM = 0.0, FDEN = 0.5, for B-R and C-R when ISTAG = 2 */
/*        FNUM = 0.5, FDEN = 0.5, for B-R and C-R when ISTAG = 2 */
/*                                in POISN2 only. */

/* ***SEE ALSO  GENBUN */
/* ***ROUTINES CALLED  PIMACH */
/* ***REVISION HISTORY  (YYMMDD) */
/*   801001  DATE WRITTEN */
/*   890531  Changed all specific intrinsics to generic.  (WRB) */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/*   900402  Added TYPE section.  (WRB) */
/* ***END PROLOGUE  COSGEN */


/* ***FIRST EXECUTABLE STATEMENT  COSGEN */
    /* Parameter adjustments */
    --a;

    /* Function Body */
    pi = pimach_(&dum);
    if (*n == 0) {
	goto L105;
    }
    if (*ijump == 1) {
	goto L103;
    }
    k3 = *n / *ijump + 1;
    k4 = k3 - 1;
    pibyn = pi / (*n + *ijump);
    i__1 = *ijump;
    for (k = 1; k <= i__1; ++k) {
	k1 = (k - 1) * k3;
	k5 = (k - 1) * k4;
	i__2 = k4;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    x = (real) (k1 + i__);
	    k2 = k5 + i__;
	    a[k2] = cos(x * pibyn) * -2.f;
/* L101: */
	}
/* L102: */
    }
    goto L105;
L103:
    np1 = *n + 1;
    y = pi / (*n + *fden);
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	x = np1 - i__ - *fnum;
	a[i__] = cos(x * y) * 2.f;
/* L104: */
    }
L105:
    return 0;
} /* cosgen_ */

/* pimach.f -- translated by f2c (version 20031025).
   You must link the resulting object file with libf2c:
	on Microsoft Windows system, link with libf2c.lib;
	on Linux or Unix systems, link with .../path/to/libf2c.a -lm
	or, if you install libf2c.a in a standard place, with -lf2c -lm
	-- in that order, at the end of the command line, as in
		cc *.o -lf2c -lm
	Source for libf2c is in /netlib/f2c/libf2c.zip, e.g.,

		http://www.netlib.org/f2c/libf2c.zip
*/

#include "f2c.h"

/* DECK PIMACH */
doublereal pimach_(real *dum)
{
    /* System generated locals */
    real ret_val;

/* ***BEGIN PROLOGUE  PIMACH */
/* ***SUBSIDIARY */
/* ***PURPOSE  Subsidiary to HSTCSP, HSTSSP and HWSCSP */
/* ***LIBRARY   SLATEC */
/* ***TYPE      SINGLE PRECISION (PIMACH-S) */
/* ***AUTHOR  (UNKNOWN) */
/* ***DESCRIPTION */

/*     This subprogram supplies the value of the constant PI correct to */
/*     machine precision where */

/*     PI=3.1415926535897932384626433832795028841971693993751058209749446 */

/* ***SEE ALSO  HSTCSP, HSTSSP, HWSCSP */
/* ***ROUTINES CALLED  (NONE) */
/* ***REVISION HISTORY  (YYMMDD) */
/*   801001  DATE WRITTEN */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/*   900402  Added TYPE section.  (WRB) */
/* ***END PROLOGUE  PIMACH */

/* ***FIRST EXECUTABLE STATEMENT  PIMACH */
    ret_val = 3.14159265358979f;
    return ret_val;
} /* pimach_ */

/* s1merg.f -- translated by f2c (version 20031025).
   You must link the resulting object file with libf2c:
	on Microsoft Windows system, link with libf2c.lib;
	on Linux or Unix systems, link with .../path/to/libf2c.a -lm
	or, if you install libf2c.a in a standard place, with -lf2c -lm
	-- in that order, at the end of the command line, as in
		cc *.o -lf2c -lm
	Source for libf2c is in /netlib/f2c/libf2c.zip, e.g.,

		http://www.netlib.org/f2c/libf2c.zip
*/

#include "f2c.h"

/* Table of constant values */

/* DECK S1MERG */
/* Subroutine */ int s1merg_(real *tcos, integer *i1, integer *m1, integer *
	i2, integer *m2, integer *i3)
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    static integer j1, j2, j3;
    extern /* Subroutine */ int scopy_(integer *, real *, integer *, real *, 
	    integer *);

/* ***BEGIN PROLOGUE  S1MERG */
/* ***SUBSIDIARY */
/* ***PURPOSE  Merge two strings of ascending real numbers. */
/* ***LIBRARY   SLATEC */
/* ***TYPE      SINGLE PRECISION (S1MERG-S, D1MERG-D, C1MERG-C, I1MERG-I) */
/* ***AUTHOR  (UNKNOWN) */
/* ***DESCRIPTION */

/*   This subroutine merges two ascending strings of numbers in the */
/*   array TCOS.  The first string is of length M1 and starts at */
/*   TCOS(I1+1).  The second string is of length M2 and starts at */
/*   TCOS(I2+1).  The merged string goes into TCOS(I3+1). */

/* ***SEE ALSO  GENBUN */
/* ***ROUTINES CALLED  SCOPY */
/* ***REVISION HISTORY  (YYMMDD) */
/*   801001  DATE WRITTEN */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/*   900402  Added TYPE section.  (WRB) */
/*   901120  Modified to use IF-THEN-ELSE.  Previous spaghetti code did */
/*           not compile correctly with optimization on the IBM RS6000. */
/*           (RWC) */
/*   920130  Code name changed from MERGE to S1MERG.  (WRB) */
/* ***END PROLOGUE  S1MERG */


/* ***FIRST EXECUTABLE STATEMENT  S1MERG */
    /* Parameter adjustments */
    --tcos;

    /* Function Body */
    if (*m1 == 0 && *m2 == 0) {
	return 0;
    }

    if (*m1 == 0 && *m2 != 0) {
	scopy_(m2, &tcos[*i2 + 1], &c__1, &tcos[*i3 + 1], &c__1);
	return 0;
    }

    if (*m1 != 0 && *m2 == 0) {
	scopy_(m1, &tcos[*i1 + 1], &c__1, &tcos[*i3 + 1], &c__1);
	return 0;
    }

    j1 = 1;
    j2 = 1;
    j3 = 1;

L10:
    if (tcos[*i1 + j1] <= tcos[*i2 + j2]) {
	tcos[*i3 + j3] = tcos[*i1 + j1];
	++j1;
	if (j1 > *m1) {
	    i__1 = *m2 - j2 + 1;
	    scopy_(&i__1, &tcos[*i2 + j2], &c__1, &tcos[*i3 + j3 + 1], &c__1);
	    return 0;
	}
    } else {
	tcos[*i3 + j3] = tcos[*i2 + j2];
	++j2;
	if (j2 > *m2) {
	    i__1 = *m1 - j1 + 1;
	    scopy_(&i__1, &tcos[*i1 + j1], &c__1, &tcos[*i3 + j3 + 1], &c__1);
	    return 0;
	}
    }
    ++j3;
    goto L10;
} /* s1merg_ */

/* trix.f -- translated by f2c (version 20031025).
   You must link the resulting object file with libf2c:
	on Microsoft Windows system, link with libf2c.lib;
	on Linux or Unix systems, link with .../path/to/libf2c.a -lm
	or, if you install libf2c.a in a standard place, with -lf2c -lm
	-- in that order, at the end of the command line, as in
		cc *.o -lf2c -lm
	Source for libf2c is in /netlib/f2c/libf2c.zip, e.g.,

		http://www.netlib.org/f2c/libf2c.zip
*/

#include "f2c.h"

/* DECK TRIX */
/* Subroutine */ int trix_(integer *idegbr, integer *idegcr, integer *m, real 
	*a, real *b, real *c__, real *y, real *tcos, real *d__, real *w)
{
    /* System generated locals */
    integer i__1, i__2;

    /* Local variables */
    static integer i__, k, l;
    static real x, z__;
    static integer kb, kc, ip;
    static real xx;
    static integer mm1, lint;

/* ***BEGIN PROLOGUE  TRIX */
/* ***SUBSIDIARY */
/* ***PURPOSE  Subsidiary to GENBUN */
/* ***LIBRARY   SLATEC */
/* ***TYPE      SINGLE PRECISION (TRIX-S, CMPTRX-C) */
/* ***AUTHOR  (UNKNOWN) */
/* ***DESCRIPTION */

/*     Subroutine to solve a system of linear equations where the */
/*     coefficient matrix is a rational function in the matrix given by */
/*     TRIDIAGONAL  ( . . . , A(I), B(I), C(I), . . . ). */

/* ***SEE ALSO  GENBUN */
/* ***ROUTINES CALLED  (NONE) */
/* ***REVISION HISTORY  (YYMMDD) */
/*   801001  DATE WRITTEN */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/*   900402  Added TYPE section.  (WRB) */
/* ***END PROLOGUE  TRIX */

/* ***FIRST EXECUTABLE STATEMENT  TRIX */
    /* Parameter adjustments */
    --w;
    --d__;
    --tcos;
    --y;
    --c__;
    --b;
    --a;

    /* Function Body */
    mm1 = *m - 1;
    kb = *idegbr + 1;
    kc = *idegcr + 1;
    l = (*idegbr + 1) / (*idegcr + 1);
    lint = 1;
    i__1 = *idegbr;
    for (k = 1; k <= i__1; ++k) {
	x = tcos[k];
	if (k != l) {
	    goto L102;
	}
	i__ = *idegbr + lint;
	xx = x - tcos[i__];
	i__2 = *m;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    w[i__] = y[i__];
	    y[i__] = xx * y[i__];
/* L101: */
	}
L102:
	z__ = 1.f / (b[1] - x);
	d__[1] = c__[1] * z__;
	y[1] *= z__;
	i__2 = mm1;
	for (i__ = 2; i__ <= i__2; ++i__) {
	    z__ = 1.f / (b[i__] - x - a[i__] * d__[i__ - 1]);
	    d__[i__] = c__[i__] * z__;
	    y[i__] = (y[i__] - a[i__] * y[i__ - 1]) * z__;
/* L103: */
	}
	z__ = b[*m] - x - a[*m] * d__[mm1];
	if (z__ != 0.f) {
	    goto L104;
	}
	y[*m] = 0.f;
	goto L105;
L104:
	y[*m] = (y[*m] - a[*m] * y[mm1]) / z__;
L105:
	i__2 = mm1;
	for (ip = 1; ip <= i__2; ++ip) {
	    i__ = *m - ip;
	    y[i__] -= d__[i__] * y[i__ + 1];
/* L106: */
	}
	if (k != l) {
	    goto L108;
	}
	i__2 = *m;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    y[i__] += w[i__];
/* L107: */
	}
	++lint;
	l = lint * kb / kc;
L108:
	;
    }
    return 0;
} /* trix_ */

/* tri3.f -- translated by f2c (version 20031025).
   You must link the resulting object file with libf2c:
	on Microsoft Windows system, link with libf2c.lib;
	on Linux or Unix systems, link with .../path/to/libf2c.a -lm
	or, if you install libf2c.a in a standard place, with -lf2c -lm
	-- in that order, at the end of the command line, as in
		cc *.o -lf2c -lm
	Source for libf2c is in /netlib/f2c/libf2c.zip, e.g.,

		http://www.netlib.org/f2c/libf2c.zip
*/

#include "f2c.h"

/* DECK TRI3 */
/* Subroutine */ int tri3_(integer *m, real *a, real *b, real *c__, integer *
	k, real *y1, real *y2, real *y3, real *tcos, real *d__, real *w1, 
	real *w2, real *w3)
{
    /* System generated locals */
    integer i__1, i__2;

    /* Local variables */
    static integer i__, n;
    static real x, z__;
    static integer k1, k2, k3, k4, l1, l2, l3, ip;
    static real xx;
    static integer mm1, k1p1, k2p1, k3p1, k4p1, k2k3k4, kint1, lint1, lint2, 
	    lint3, kint2, kint3;

/* ***BEGIN PROLOGUE  TRI3 */
/* ***SUBSIDIARY */
/* ***PURPOSE  Subsidiary to GENBUN */
/* ***LIBRARY   SLATEC */
/* ***TYPE      SINGLE PRECISION (TRI3-S, CMPTR3-C) */
/* ***AUTHOR  (UNKNOWN) */
/* ***DESCRIPTION */

/*     Subroutine to solve three linear systems whose common coefficient */
/*     matrix is a rational function in the matrix given by */

/*                  TRIDIAGONAL (...,A(I),B(I),C(I),...) */

/* ***SEE ALSO  GENBUN */
/* ***ROUTINES CALLED  (NONE) */
/* ***REVISION HISTORY  (YYMMDD) */
/*   801001  DATE WRITTEN */
/*   890206  REVISION DATE from Version 3.2 */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/*   900402  Added TYPE section.  (WRB) */
/* ***END PROLOGUE  TRI3 */

/* ***FIRST EXECUTABLE STATEMENT  TRI3 */
    /* Parameter adjustments */
    --w3;
    --w2;
    --w1;
    --d__;
    --tcos;
    --y3;
    --y2;
    --y1;
    --k;
    --c__;
    --b;
    --a;

    /* Function Body */
    mm1 = *m - 1;
    k1 = k[1];
    k2 = k[2];
    k3 = k[3];
    k4 = k[4];
    k1p1 = k1 + 1;
    k2p1 = k2 + 1;
    k3p1 = k3 + 1;
    k4p1 = k4 + 1;
    k2k3k4 = k2 + k3 + k4;
    if (k2k3k4 == 0) {
	goto L101;
    }
    l1 = (k1 + 1) / (k2 + 1);
    l2 = (k1 + 1) / (k3 + 1);
    l3 = (k1 + 1) / (k4 + 1);
    lint1 = 1;
    lint2 = 1;
    lint3 = 1;
    kint1 = k1;
    kint2 = kint1 + k2;
    kint3 = kint2 + k3;
L101:
    i__1 = k1;
    for (n = 1; n <= i__1; ++n) {
	x = tcos[n];
	if (k2k3k4 == 0) {
	    goto L107;
	}
	if (n != l1) {
	    goto L103;
	}
	i__2 = *m;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    w1[i__] = y1[i__];
/* L102: */
	}
L103:
	if (n != l2) {
	    goto L105;
	}
	i__2 = *m;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    w2[i__] = y2[i__];
/* L104: */
	}
L105:
	if (n != l3) {
	    goto L107;
	}
	i__2 = *m;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    w3[i__] = y3[i__];
/* L106: */
	}
L107:
	z__ = 1.f / (b[1] - x);
	d__[1] = c__[1] * z__;
	y1[1] *= z__;
	y2[1] *= z__;
	y3[1] *= z__;
	i__2 = *m;
	for (i__ = 2; i__ <= i__2; ++i__) {
	    z__ = 1.f / (b[i__] - x - a[i__] * d__[i__ - 1]);
	    d__[i__] = c__[i__] * z__;
	    y1[i__] = (y1[i__] - a[i__] * y1[i__ - 1]) * z__;
	    y2[i__] = (y2[i__] - a[i__] * y2[i__ - 1]) * z__;
	    y3[i__] = (y3[i__] - a[i__] * y3[i__ - 1]) * z__;
/* L108: */
	}
	i__2 = mm1;
	for (ip = 1; ip <= i__2; ++ip) {
	    i__ = *m - ip;
	    y1[i__] -= d__[i__] * y1[i__ + 1];
	    y2[i__] -= d__[i__] * y2[i__ + 1];
	    y3[i__] -= d__[i__] * y3[i__ + 1];
/* L109: */
	}
	if (k2k3k4 == 0) {
	    goto L115;
	}
	if (n != l1) {
	    goto L111;
	}
	i__ = lint1 + kint1;
	xx = x - tcos[i__];
	i__2 = *m;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    y1[i__] = xx * y1[i__] + w1[i__];
/* L110: */
	}
	++lint1;
	l1 = lint1 * k1p1 / k2p1;
L111:
	if (n != l2) {
	    goto L113;
	}
	i__ = lint2 + kint2;
	xx = x - tcos[i__];
	i__2 = *m;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    y2[i__] = xx * y2[i__] + w2[i__];
/* L112: */
	}
	++lint2;
	l2 = lint2 * k1p1 / k3p1;
L113:
	if (n != l3) {
	    goto L115;
	}
	i__ = lint3 + kint3;
	xx = x - tcos[i__];
	i__2 = *m;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    y3[i__] = xx * y3[i__] + w3[i__];
/* L114: */
	}
	++lint3;
	l3 = lint3 * k1p1 / k4p1;
L115:
	;
    }
    return 0;
} /* tri3_ */

/* scopy.f -- translated by f2c (version 20031025).
   You must link the resulting object file with libf2c:
	on Microsoft Windows system, link with libf2c.lib;
	on Linux or Unix systems, link with .../path/to/libf2c.a -lm
	or, if you install libf2c.a in a standard place, with -lf2c -lm
	-- in that order, at the end of the command line, as in
		cc *.o -lf2c -lm
	Source for libf2c is in /netlib/f2c/libf2c.zip, e.g.,

		http://www.netlib.org/f2c/libf2c.zip
*/

#include "f2c.h"

/* Subroutine */ int scopy_(integer *n, real *sx, integer *incx, real *sy, 
	integer *incy)
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    static integer i__, m, ix, iy, mp1;


/*     copies a vector, x, to a vector, y. */
/*     uses unrolled loops for increments equal to 1. */
/*     jack dongarra, linpack, 3/11/78. */
/*     modified 12/3/93, array(1) declarations changed to array(*) */


    /* Parameter adjustments */
    --sy;
    --sx;

    /* Function Body */
    if (*n <= 0) {
	return 0;
    }
    if (*incx == 1 && *incy == 1) {
	goto L20;
    }

/*        code for unequal increments or equal increments */
/*          not equal to 1 */

    ix = 1;
    iy = 1;
    if (*incx < 0) {
	ix = (-(*n) + 1) * *incx + 1;
    }
    if (*incy < 0) {
	iy = (-(*n) + 1) * *incy + 1;
    }
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	sy[iy] = sx[ix];
	ix += *incx;
	iy += *incy;
/* L10: */
    }
    return 0;

/*        code for both increments equal to 1 */


/*        clean-up loop */

L20:
    m = *n % 7;
    if (m == 0) {
	goto L40;
    }
    i__1 = m;
    for (i__ = 1; i__ <= i__1; ++i__) {
	sy[i__] = sx[i__];
/* L30: */
    }
    if (*n < 7) {
	return 0;
    }
L40:
    mp1 = m + 1;
    i__1 = *n;
    for (i__ = mp1; i__ <= i__1; i__ += 7) {
	sy[i__] = sx[i__];
	sy[i__ + 1] = sx[i__ + 1];
	sy[i__ + 2] = sx[i__ + 2];
	sy[i__ + 3] = sx[i__ + 3];
	sy[i__ + 4] = sx[i__ + 4];
	sy[i__ + 5] = sx[i__ + 5];
	sy[i__ + 6] = sx[i__ + 6];
/* L50: */
    }
    return 0;
} /* scopy_ */

/* gr3.f -- translated by f2c (version 20000121).
   You must link the resulting object file with the libraries:
	-lf2c -lm   (in that order)
*/

#include "f2c.h"

/*   Fishpack calls routines from fftpack.  If you don't already have */
/*   those routines, submit your request in the form */
/*        send blktri from fishpack fftpack */
/* aveat receptor.  (Jack) dongarra@anl-mcs.arpa, (Eric Grosse) research!ehg */
/* areful! Anything free comes with no guarantee. */
/* *** from netlib, Fri Jun 10 16:46:03 CDT 1988 *** */
/* fishpk12 from portlib                                  12/30/83 */
/* Subroutine */ int hw3crt_(real *xs, real *xf, integer *l, integer *lbdcnd, 
	real *bdxs, real *bdxf, real *ys, real *yf, integer *m, integer *
	mbdcnd, real *bdys, real *bdyf, real *zs, real *zf, integer *n, 
	integer *nbdcnd, real *bdzs, real *bdzf, real *elmbda, integer *ldimf,
	 integer *mdimf, real *f, real *pertrb, integer *ierror, real *w)
{
    /* System generated locals */
    integer bdxs_dim1, bdxs_offset, bdxf_dim1, bdxf_offset, bdys_dim1, 
	    bdys_offset, bdyf_dim1, bdyf_offset, bdzs_dim1, bdzs_offset, 
	    bdzf_dim1, bdzf_offset, f_dim1, f_dim2, f_offset, i__1, i__2, 
	    i__3;
    real r__1;

    /* Local variables */
    static integer lunk, munk, nunk, i__, j, k;
    static real s, c1, c2, c3;
    static integer lstop, mstop, nstop;
    static real s1, s2;
    extern /* Subroutine */ int pois3d_(integer *, integer *, real *, integer 
	    *, integer *, real *, integer *, integer *, real *, real *, real *
	    , integer *, integer *, real *, integer *, real *);
    static integer lstpm1, mstpm1, nstpm1, ir;
    static real dx, dy;
    static integer mp;
    static real dz;
    static integer np, lp, nperod, lstart, mstart, nstart, lp1, mp1, np1;
    static real twbydx, twbydy, twbydz;
    static integer iwb, iwc;
    static real xlp, ylp, zlp;
    static integer iww;



/*     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*     *                                                               * */
/*     *                        f i s h p a k                          * */
/*     *                                                               * */
/*     *                                                               * */
/*     *     a package of fortran subprograms for the solution of      * */
/*     *                                                               * */
/*     *      separable elliptic partial differential equations        * */
/*     *                                                               * */
/*     *                  (version 3.1 , october 1980)                  * */
/*     *                                                               * */
/*     *                             by                                * */
/*     *                                                               * */
/*     *        john adams, paul swarztrauber and roland sweet         * */
/*     *                                                               * */
/*     *                             of                                * */
/*     *                                                               * */
/*     *         the national center for atmospheric research          * */
/*     *                                                               * */
/*     *                boulder, colorado  (80307)  u.s.a.             * */
/*     *                                                               * */
/*     *                   which is sponsored by                       * */
/*     *                                                               * */
/*     *              the national science foundation                  * */
/*     *                                                               * */
/*     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */


/*    * * * * * * * * *  purpose    * * * * * * * * * * * * * * * * * * */

/*          subroutine hw3crt solves the standard seven-point finite */
/*     difference approximation to the helmholtz equation in cartesian */
/*     coordinates: */

/*         (d/dx)(du/dx) + (d/dy)(du/dy) + (d/dz)(du/dz) */

/*                    + lambda*u = f(x,y,z) . */

/*    * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */


/*    * * * * * * * *    parameter description     * * * * * * * * * * */


/*            * * * * * *   on input    * * * * * * */

/*     xs,xf */
/*        the range of x, i.e. xs .le. x .le. xf . */
/*        xs must be less than xf. */

/*     l */
/*        the number of panels into which the interval (xs,xf) is */
/*        subdivided.  hence, there will be l+1 grid points in the */
/*        x-direction given by x(i) = xs+(i-1)dx for i=1,2,...,l+1, */
/*        where dx = (xf-xs)/l is the panel width.  l must be at */
/*        least 5 . */

/*     lbdcnd */
/*        indicates the type of boundary conditions at x = xs and x = xf. */

/*        = 0  if the solution is periodic in x, i.e. */
/*             u(l+i,j,k) = u(i,j,k). */
/*        = 1  if the solution is specified at x = xs and x = xf. */
/*        = 2  if the solution is specified at x = xs and the derivative */
/*             of the solution with respect to x is specified at x = xf. */
/*        = 3  if the derivative of the solution with respect to x is */
/*             specified at x = xs and x = xf. */
/*        = 4  if the derivative of the solution with respect to x is */
/*             specified at x = xs and the solution is specified at x=xf. */

/*     bdxs */
/*        a two-dimensional array that specifies the values of the */
/*        derivative of the solution with respect to x at x = xs. */
/*        when lbdcnd = 3 or 4, */

/*             bdxs(j,k) = (d/dx)u(xs,y(j),z(k)), j=1,2,...,m+1, */
/*                                                k=1,2,...,n+1. */

/*        when lbdcnd has any other value, bdxs is a dummy variable. */
/*        bdxs must be dimensioned at least (m+1)*(n+1). */

/*     bdxf */
/*        a two-dimensional array that specifies the values of the */
/*        derivative of the solution with respect to x at x = xf. */
/*        when lbdcnd = 2 or 3, */

/*             bdxf(j,k) = (d/dx)u(xf,y(j),z(k)), j=1,2,...,m+1, */
/*                                                k=1,2,...,n+1. */

/*        when lbdcnd has any other value, bdxf is a dummy variable. */
/*        bdxf must be dimensioned at least (m+1)*(n+1). */

/*     ys,yf */
/*        the range of y, i.e. ys .le. y .le. yf. */
/*        ys must be less than yf. */

/*     m */
/*        the number of panels into which the interval (ys,yf) is */
/*        subdivided.  hence, there will be m+1 grid points in the */
/*        y-direction given by y(j) = ys+(j-1)dy for j=1,2,...,m+1, */
/*        where dy = (yf-ys)/m is the panel width.  m must be at */
/*        least 5 . */

/*     mbdcnd */
/*        indicates the type of boundary conditions at y = ys and y = yf. */

/*        = 0  if the solution is periodic in y, i.e. */
/*             u(i,m+j,k) = u(i,j,k). */
/*        = 1  if the solution is specified at y = ys and y = yf. */
/*        = 2  if the solution is specified at y = ys and the derivative */
/*             of the solution with respect to y is specified at y = yf. */
/*        = 3  if the derivative of the solution with respect to y is */
/*             specified at y = ys and y = yf. */
/*        = 4  if the derivative of the solution with respect to y is */
/*             specified at y = ys and the solution is specified at y=yf. */

/*     bdys */
/*        a two-dimensional array that specifies the values of the */
/*        derivative of the solution with respect to y at y = ys. */
/*        when mbdcnd = 3 or 4, */

/*             bdys(i,k) = (d/dy)u(x(i),ys,z(k)), i=1,2,...,l+1, */
/*                                                k=1,2,...,n+1. */

/*        when mbdcnd has any other value, bdys is a dummy variable. */
/*        bdys must be dimensioned at least (l+1)*(n+1). */

/*     bdyf */
/*        a two-dimensional array that specifies the values of the */
/*        derivative of the solution with respect to y at y = yf. */
/*        when mbdcnd = 2 or 3, */

/*             bdyf(i,k) = (d/dy)u(x(i),yf,z(k)), i=1,2,...,l+1, */
/*                                                k=1,2,...,n+1. */

/*        when mbdcnd has any other value, bdyf is a dummy variable. */
/*        bdyf must be dimensioned at least (l+1)*(n+1). */

/*     zs,zf */
/*        the range of z, i.e. zs .le. z .le. zf. */
/*        zs must be less than zf. */

/*     n */
/*        the number of panels into which the interval (zs,zf) is */
/*        subdivided.  hence, there will be n+1 grid points in the */
/*        z-direction given by z(k) = zs+(k-1)dz for k=1,2,...,n+1, */
/*        where dz = (zf-zs)/n is the panel width.  n must be at least 5. */

/*     nbdcnd */
/*        indicates the type of boundary conditions at z = zs and z = zf. */

/*        = 0  if the solution is periodic in z, i.e. */
/*             u(i,j,n+k) = u(i,j,k). */
/*        = 1  if the solution is specified at z = zs and z = zf. */
/*        = 2  if the solution is specified at z = zs and the derivative */
/*             of the solution with respect to z is specified at z = zf. */
/*        = 3  if the derivative of the solution with respect to z is */
/*             specified at z = zs and z = zf. */
/*        = 4  if the derivative of the solution with respect to z is */
/*             specified at z = zs and the solution is specified at z=zf. */

/*     bdzs */
/*        a two-dimensional array that specifies the values of the */
/*        derivative of the solution with respect to z at z = zs. */
/*        when nbdcnd = 3 or 4, */

/*             bdzs(i,j) = (d/dz)u(x(i),y(j),zs), i=1,2,...,l+1, */
/*                                                j=1,2,...,m+1. */

/*        when nbdcnd has any other value, bdzs is a dummy variable. */
/*        bdzs must be dimensioned at least (l+1)*(m+1). */

/*     bdzf */
/*        a two-dimensional array that specifies the values of the */
/*        derivative of the solution with respect to z at z = zf. */
/*        when nbdcnd = 2 or 3, */

/*             bdzf(i,j) = (d/dz)u(x(i),y(j),zf), i=1,2,...,l+1, */
/*                                                j=1,2,...,m+1. */

/*        when nbdcnd has any other value, bdzf is a dummy variable. */
/*        bdzf must be dimensioned at least (l+1)*(m+1). */

/*     elmbda */
/*        the constant lambda in the helmholtz equation. if */
/*        lambda .gt. 0, a solution may not exist.  however, hw3crt will */
/*        attempt to find a solution. */

/*     f */
/*        a three-dimensional array that specifies the values of the */
/*        right side of the helmholtz equation and boundary values (if */
/*        any).  for i=2,3,...,l, j=2,3,...,m, and k=2,3,...,n */

/*                   f(i,j,k) = f(x(i),y(j),z(k)). */

/*        on the boundaries f is defined by */

/*        lbdcnd      f(1,j,k)         f(l+1,j,k) */
/*        ------   ---------------   --------------- */

/*          0      f(xs,y(j),z(k))   f(xs,y(j),z(k)) */
/*          1      u(xs,y(j),z(k))   u(xf,y(j),z(k)) */
/*          2      u(xs,y(j),z(k))   f(xf,y(j),z(k))   j=1,2,...,m+1 */
/*          3      f(xs,y(j),z(k))   f(xf,y(j),z(k))   k=1,2,...,n+1 */
/*          4      f(xs,y(j),z(k))   u(xf,y(j),z(k)) */

/*        mbdcnd      f(i,1,k)         f(i,m+1,k) */
/*        ------   ---------------   --------------- */

/*          0      f(x(i),ys,z(k))   f(x(i),ys,z(k)) */
/*          1      u(x(i),ys,z(k))   u(x(i),yf,z(k)) */
/*          2      u(x(i),ys,z(k))   f(x(i),yf,z(k))   i=1,2,...,l+1 */
/*          3      f(x(i),ys,z(k))   f(x(i),yf,z(k))   k=1,2,...,n+1 */
/*          4      f(x(i),ys,z(k))   u(x(i),yf,z(k)) */

/*        nbdcnd      f(i,j,1)         f(i,j,n+1) */
/*        ------   ---------------   --------------- */

/*          0      f(x(i),y(j),zs)   f(x(i),y(j),zs) */
/*          1      u(x(i),y(j),zs)   u(x(i),y(j),zf) */
/*          2      u(x(i),y(j),zs)   f(x(i),y(j),zf)   i=1,2,...,l+1 */
/*          3      f(x(i),y(j),zs)   f(x(i),y(j),zf)   j=1,2,...,m+1 */
/*          4      f(x(i),y(j),zs)   u(x(i),y(j),zf) */

/*        f must be dimensioned at least (l+1)*(m+1)*(n+1). */

/*        note: */

/*        if the table calls for both the solution u and the right side f */
/*        on a boundary, then the solution must be specified. */

/*     ldimf */
/*        the row (or first) dimension of the arrays f,bdys,bdyf,bdzs, */
/*        and bdzf as it appears in the program calling hw3crt. this */
/*        parameter is used to specify the variable dimension of these */
/*        arrays.  ldimf must be at least l+1. */

/*     mdimf */
/*        the column (or second) dimension of the array f and the row (or */
/*        first) dimension of the arrays bdxs and bdxf as it appears in */
/*        the program calling hw3crt.  this parameter is used to specify */
/*        the variable dimension of these arrays. */
/*        mdimf must be at least m+1. */

/*     w */
/*        a one-dimensional array that must be provided by the user for */
/*        work space.  the length of w must be at least 30 + l + m + 5*n */
/*        + max(l,m,n) + 7*(int((l+1)/2) + int((m+1)/2)) */


/*            * * * * * *   on output   * * * * * * */

/*     f */
/*        contains the solution u(i,j,k) of the finite difference */
/*        approximation for the grid point (x(i),y(j),z(k)) for */
/*        i=1,2,...,l+1, j=1,2,...,m+1, and k=1,2,...,n+1. */

/*     pertrb */
/*        if a combination of periodic or derivative boundary conditions */
/*        is specified for a poisson equation (lambda = 0), a solution */
/*        may not exist.  pertrb is a constant, calculated and subtracted */
/*        from f, which ensures that a solution exists.  pwscrt then */
/*        computes this solution, which is a least squares solution to */
/*        the original approximation.  this solution is not unique and is */
/*        unnormalized.  the value of pertrb should be small compared to */
/*        the right side f.  otherwise, a solution is obtained to an */
/*        essentially different problem.  this comparison should always */
/*        be made to insure that a meaningful solution has been obtained. */

/*     ierror */
/*        an error flag that indicates invalid input parameters.  except */
/*        for numbers 0 and 12, a solution is not attempted. */

/*        =  0  no error */
/*        =  1  xs .ge. xf */
/*        =  2  l .lt. 5 */
/*        =  3  lbdcnd .lt. 0 .or. lbdcnd .gt. 4 */
/*        =  4  ys .ge. yf */
/*        =  5  m .lt. 5 */
/*        =  6  mbdcnd .lt. 0 .or. mbdcnd .gt. 4 */
/*        =  7  zs .ge. zf */
/*        =  8  n .lt. 5 */
/*        =  9  nbdcnd .lt. 0 .or. nbdcnd .gt. 4 */
/*        = 10  ldimf .lt. l+1 */
/*        = 11  mdimf .lt. m+1 */
/*        = 12  lambda .gt. 0 */

/*        since this is the only means of indicating a possibly incorrect */
/*        call to hw3crt, the user should test ierror after the call. */


/*    * * * * * * *   program specifications    * * * * * * * * * * * * */

/*     dimension of   bdxs(mdimf,n+1),bdxf(mdimf,n+1),bdys(ldimf,n+1), */
/*     arguments      bdyf(ldimf,n+1),bdzs(ldimf,m+1),bdzf(ldimf,m+1), */
/*                    f(ldimf,mdimf,n+1),w(see argument list) */

/*     latest         december 1, 1978 */
/*     revision */

/*     subprograms    hw3crt,pois3d,pos3d1,trid,rffti,rfftf,rfftf1, */
/*     required       rfftb,rfftb1,costi,cost,sinti,sint,cosqi,cosqf, */
/*                    cosqf1,cosqb,cosqb1,sinqi,sinqf,sinqb,cffti, */
/*                    cffti1,cfftb,cfftb1,passb2,passb3,passb4,passb, */
/*                    cfftf,cfftf1,passf1,passf2,passf3,passf4,passf, */
/*                    pimach */

/*     special        none */
/*     conditions */

/*     common         value */
/*     blocks */

/*     i/o            none */

/*     precision      single */

/*     specialist     roland sweet */

/*     language       fortran */

/*     history        written by roland sweet at ncar in july,1977 */

/*     algorithm      this subroutine defines the finite difference */
/*                    equations, incorporates boundary data, and */
/*                    adjusts the right side of singular systems and */
/*                    then calls pois3d to solve the system. */

/*     space          7862(decimal) = 17300(octal) locations on the */
/*     required       ncar control data 7600 */

/*     timing and        the execution time t on the ncar control data */
/*     accuracy       7600 for subroutine hw3crt is roughly proportional */
/*                    to l*m*n*(log2(l)+log2(m)+5), but also depends on */
/*                    input parameters lbdcnd and mbdcnd.  some typical */
/*                    values are listed in the table below. */
/*                       the solution process employed results in a loss */
/*                    of no more than three significant digits for l,m an */
/*                    n as large as 32.  more detailed information about */
/*                    accuracy can be found in the documentation for */
/*                    subroutine pois3d which is the routine that actuall */
/*                    solves the finite difference equations. */


/*                       l(=m=n)     lbdcnd(=mbdcnd=nbdcnd)      t(msecs) */
/*                       -------     ----------------------      -------- */

/*                         16                  0                    300 */
/*                         16                  1                    302 */
/*                         16                  3                    348 */
/*                         32                  0                   1925 */
/*                         32                  1                   1929 */
/*                         32                  3                   2109 */

/*     portability    american national standards institute fortran. */
/*                    the machine dependent constant pi is defined in */
/*                    function pimach. */

/*     required       cos,sin,atan */
/*     resident */
/*     routines */

/*     reference      none */

/*     required         cos,sin,atan */
/*     resident */
/*     routines */

/*     reference        none */

/*    * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */


/*     check for invalid input. */

    /* Parameter adjustments */
    bdzf_dim1 = *ldimf;
    bdzf_offset = 1 + bdzf_dim1 * 1;
    bdzf -= bdzf_offset;
    bdzs_dim1 = *ldimf;
    bdzs_offset = 1 + bdzs_dim1 * 1;
    bdzs -= bdzs_offset;
    bdyf_dim1 = *ldimf;
    bdyf_offset = 1 + bdyf_dim1 * 1;
    bdyf -= bdyf_offset;
    bdys_dim1 = *ldimf;
    bdys_offset = 1 + bdys_dim1 * 1;
    bdys -= bdys_offset;
    f_dim1 = *ldimf;
    f_dim2 = *mdimf;
    f_offset = 1 + f_dim1 * (1 + f_dim2 * 1);
    f -= f_offset;
    bdxf_dim1 = *mdimf;
    bdxf_offset = 1 + bdxf_dim1 * 1;
    bdxf -= bdxf_offset;
    bdxs_dim1 = *mdimf;
    bdxs_offset = 1 + bdxs_dim1 * 1;
    bdxs -= bdxs_offset;
    --w;

    /* Function Body */

    // printf("xs and xf are %lf %lf\n",*xs,*xf);
    // printf("ys and yf are %lf %lf\n",*ys,*yf);
    // printf("zs and zf are %lf %lf\n",*zs,*zf);
    // printf("f_dim1 %d, f_dim2 %d, f_offset %d\n",f_dim1,f_dim2,f_offset);
    // printf("f(0,0,0) is f[0] is %7.1f\n",f[0]);
    // printf("f(0,0,1) is f[1] is %7.1f\n",f[1]);
    // printf("f(0,1,0) is f[33] is %7.1f\n",f[33]);
    // printf("f(0,1,1) is f[34] is %7.1f\n",f[34]);
    // printf("f(1,0,0) is f[1089] is %7.1f\n",f[1089]);
    // printf("f(1,0,1) is f[1090] is %7.1f\n",f[1090]);
    // printf("f(1,1,1) is f[1123] is %7.1f\n",f[1123]);
    // printf("f(1,1,2) is f[1124] is %7.1f\n",f[1124]);
    // printf("f(2,3,4) is f[2280] is %7.1f\n",f[2280]);
    // printf("f(5,11,27) is f[5834] is %7.1f\n",f[5834]);


    *ierror = 0;
    if (*xf <= *xs) {
	*ierror = 1;
    }
    if (*l < 5) {
	*ierror = 2;
    }
    if (*lbdcnd < 0 || *lbdcnd > 4) {
	*ierror = 3;
    }
    if (*yf <= *ys) {
	*ierror = 4;
    }
    if (*m < 5) {
	*ierror = 5;
    }
    if (*mbdcnd < 0 || *mbdcnd > 4) {
	*ierror = 6;
    }
    if (*zf <= *zs) {
	*ierror = 7;
    }
    if (*n < 5) {
	*ierror = 8;
    }
    if (*nbdcnd < 0 || *nbdcnd > 4) {
	*ierror = 9;
    }
    if (*ldimf < *l + 1) {
	*ierror = 10;
    }
    if (*mdimf < *m + 1) {
	*ierror = 11;
    }
    if (*ierror != 0) {
	goto L188;
    }
    dy = (*yf - *ys) / (real) (*m);
    twbydy = 2.f / dy;
/* Computing 2nd power */
    r__1 = dy;
    c2 = 1.f / (r__1 * r__1);
    mstart = 1;
    mstop = *m;
    mp1 = *m + 1;
    mp = *mbdcnd + 1;
    switch (mp) {
	case 1:  goto L104;
	case 2:  goto L101;
	case 3:  goto L101;
	case 4:  goto L102;
	case 5:  goto L102;
    }
L101:
    mstart = 2;
L102:
    switch (mp) {
	case 1:  goto L104;
	case 2:  goto L104;
	case 3:  goto L103;
	case 4:  goto L103;
	case 5:  goto L104;
    }
L103:
    mstop = mp1;
L104:
    munk = mstop - mstart + 1;
    dz = (*zf - *zs) / (real) (*n);
    twbydz = 2.f / dz;
    np = *nbdcnd + 1;
/* Computing 2nd power */
    r__1 = dz;
    c3 = 1.f / (r__1 * r__1);
    np1 = *n + 1;
    nstart = 1;
    nstop = *n;
    switch (np) {
	case 1:  goto L108;
	case 2:  goto L105;
	case 3:  goto L105;
	case 4:  goto L106;
	case 5:  goto L106;
    }
L105:
    nstart = 2;
L106:
    switch (np) {
	case 1:  goto L108;
	case 2:  goto L108;
	case 3:  goto L107;
	case 4:  goto L107;
	case 5:  goto L108;
    }
L107:
    nstop = np1;
L108:
    nunk = nstop - nstart + 1;
    lp1 = *l + 1;
    dx = (*xf - *xs) / (real) (*l);
/* Computing 2nd power */
    r__1 = dx;
    c1 = 1.f / (r__1 * r__1);
    twbydx = 2.f / dx;
    lp = *lbdcnd + 1;
    lstart = 1;
    lstop = *l;

/*     enter boundary data for x-boundaries. */

    switch (lp) {
	case 1:  goto L122;
	case 2:  goto L109;
	case 3:  goto L109;
	case 4:  goto L112;
	case 5:  goto L112;
    }
L109:
    lstart = 2;
    i__1 = mstop;
    for (j = mstart; j <= i__1; ++j) {
	i__2 = nstop;
	for (k = nstart; k <= i__2; ++k) {
	    f[(j + k * f_dim2) * f_dim1 + 2] -= c1 * f[(j + k * f_dim2) * 
		    f_dim1 + 1];
/* L110: */
	}
/* L111: */
    }
    goto L115;
L112:
    i__1 = mstop;
    for (j = mstart; j <= i__1; ++j) {
	i__2 = nstop;
	for (k = nstart; k <= i__2; ++k) {
	    f[(j + k * f_dim2) * f_dim1 + 1] += twbydx * bdxs[j + k * 
		    bdxs_dim1];
/* L113: */
	}
/* L114: */
    }
L115:
    switch (lp) {
	case 1:  goto L122;
	case 2:  goto L116;
	case 3:  goto L119;
	case 4:  goto L119;
	case 5:  goto L116;
    }
L116:
    i__1 = mstop;
    for (j = mstart; j <= i__1; ++j) {
	i__2 = nstop;
	for (k = nstart; k <= i__2; ++k) {
	    f[*l + (j + k * f_dim2) * f_dim1] -= c1 * f[lp1 + (j + k * f_dim2)
		     * f_dim1];
/* L117: */
	}
/* L118: */
    }
    goto L122;
L119:
    lstop = lp1;
    i__1 = mstop;
    for (j = mstart; j <= i__1; ++j) {
	i__2 = nstop;
	for (k = nstart; k <= i__2; ++k) {
	    f[lp1 + (j + k * f_dim2) * f_dim1] -= twbydx * bdxf[j + k * 
		    bdxf_dim1];
/* L120: */
	}
/* L121: */
    }
L122:
    lunk = lstop - lstart + 1;

/*     enter boundary data for y-boundaries. */

    switch (mp) {
	case 1:  goto L136;
	case 2:  goto L123;
	case 3:  goto L123;
	case 4:  goto L126;
	case 5:  goto L126;
    }
L123:
    i__1 = lstop;
    for (i__ = lstart; i__ <= i__1; ++i__) {
	i__2 = nstop;
	for (k = nstart; k <= i__2; ++k) {
	    f[i__ + (k * f_dim2 + 2) * f_dim1] -= c2 * f[i__ + (k * f_dim2 + 
		    1) * f_dim1];
/* L124: */
	}
/* L125: */
    }
    goto L129;
L126:
    i__1 = lstop;
    for (i__ = lstart; i__ <= i__1; ++i__) {
	i__2 = nstop;
	for (k = nstart; k <= i__2; ++k) {
	    f[i__ + (k * f_dim2 + 1) * f_dim1] += twbydy * bdys[i__ + k * 
		    bdys_dim1];
/* L127: */
	}
/* L128: */
    }
L129:
    switch (mp) {
	case 1:  goto L136;
	case 2:  goto L130;
	case 3:  goto L133;
	case 4:  goto L133;
	case 5:  goto L130;
    }
L130:
    i__1 = lstop;
    for (i__ = lstart; i__ <= i__1; ++i__) {
	i__2 = nstop;
	for (k = nstart; k <= i__2; ++k) {
	    f[i__ + (*m + k * f_dim2) * f_dim1] -= c2 * f[i__ + (mp1 + k * 
		    f_dim2) * f_dim1];
/* L131: */
	}
/* L132: */
    }
    goto L136;
L133:
    i__1 = lstop;
    for (i__ = lstart; i__ <= i__1; ++i__) {
	i__2 = nstop;
	for (k = nstart; k <= i__2; ++k) {
	    f[i__ + (mp1 + k * f_dim2) * f_dim1] -= twbydy * bdyf[i__ + k * 
		    bdyf_dim1];
/* L134: */
	}
/* L135: */
    }
L136:

/*     enter boundary data for z-boundaries. */

    switch (np) {
	case 1:  goto L150;
	case 2:  goto L137;
	case 3:  goto L137;
	case 4:  goto L140;
	case 5:  goto L140;
    }
L137:
    i__1 = lstop;
    for (i__ = lstart; i__ <= i__1; ++i__) {
	i__2 = mstop;
	for (j = mstart; j <= i__2; ++j) {
	    f[i__ + (j + (f_dim2 << 1)) * f_dim1] -= c3 * f[i__ + (j + f_dim2)
		     * f_dim1];
/* L138: */
	}
/* L139: */
    }
    goto L143;
L140:
    i__1 = lstop;
    for (i__ = lstart; i__ <= i__1; ++i__) {
	i__2 = mstop;
	for (j = mstart; j <= i__2; ++j) {
	    f[i__ + (j + f_dim2) * f_dim1] += twbydz * bdzs[i__ + j * 
		    bdzs_dim1];
/* L141: */
	}
/* L142: */
    }
L143:
    switch (np) {
	case 1:  goto L150;
	case 2:  goto L144;
	case 3:  goto L147;
	case 4:  goto L147;
	case 5:  goto L144;
    }
L144:
    i__1 = lstop;
    for (i__ = lstart; i__ <= i__1; ++i__) {
	i__2 = mstop;
	for (j = mstart; j <= i__2; ++j) {
	    f[i__ + (j + *n * f_dim2) * f_dim1] -= c3 * f[i__ + (j + np1 * 
		    f_dim2) * f_dim1];
/* L145: */
	}
/* L146: */
    }
    goto L150;
L147:
    i__1 = lstop;
    for (i__ = lstart; i__ <= i__1; ++i__) {
	i__2 = mstop;
	for (j = mstart; j <= i__2; ++j) {
	    f[i__ + (j + np1 * f_dim2) * f_dim1] -= twbydz * bdzf[i__ + j * 
		    bdzf_dim1];
/* L148: */
	}
/* L149: */
    }

/*     define a,b,c coefficients in w-array. */

L150:
    iwb = nunk + 1;
    iwc = iwb + nunk;
    iww = iwc + nunk;
    i__1 = nunk;
    for (k = 1; k <= i__1; ++k) {
	i__ = iwc + k - 1;
	w[k] = c3;
	w[i__] = c3;
	i__ = iwb + k - 1;
	w[i__] = c3 * -2.f + *elmbda;
/* L151: */
    }
    switch (np) {
	case 1:  goto L155;
	case 2:  goto L155;
	case 3:  goto L153;
	case 4:  goto L152;
	case 5:  goto L152;
    }
L152:
    w[iwc] = c3 * 2.f;
L153:
    switch (np) {
	case 1:  goto L155;
	case 2:  goto L155;
	case 3:  goto L154;
	case 4:  goto L154;
	case 5:  goto L155;
    }
L154:
    w[iwb - 1] = c3 * 2.f;
L155:
    *pertrb = 0.f;

/*     for singular problems adjust data to insure a solution will exist. */

    switch (lp) {
	case 1:  goto L156;
	case 2:  goto L172;
	case 3:  goto L172;
	case 4:  goto L156;
	case 5:  goto L172;
    }
L156:
    switch (mp) {
	case 1:  goto L157;
	case 2:  goto L172;
	case 3:  goto L172;
	case 4:  goto L157;
	case 5:  goto L172;
    }
L157:
    switch (np) {
	case 1:  goto L158;
	case 2:  goto L172;
	case 3:  goto L172;
	case 4:  goto L158;
	case 5:  goto L172;
    }
L158:
    if (*elmbda < 0.f) {
	goto L172;
    } else if (*elmbda == 0) {
	goto L160;
    } else {
	goto L159;
    }
L159:
    *ierror = 12;
    goto L172;
L160:
    mstpm1 = mstop - 1;
    lstpm1 = lstop - 1;
    nstpm1 = nstop - 1;
    xlp = (real) ((lp + 2) / 3);
    ylp = (real) ((mp + 2) / 3);
    zlp = (real) ((np + 2) / 3);
    s1 = 0.f;
    i__1 = nstpm1;
    for (k = 2; k <= i__1; ++k) {
	i__2 = mstpm1;
	for (j = 2; j <= i__2; ++j) {
	    i__3 = lstpm1;
	    for (i__ = 2; i__ <= i__3; ++i__) {
		s1 += f[i__ + (j + k * f_dim2) * f_dim1];
/* L161: */
	    }
	    s1 += (f[(j + k * f_dim2) * f_dim1 + 1] + f[lstop + (j + k * 
		    f_dim2) * f_dim1]) / xlp;
/* L162: */
	}
	s2 = 0.f;
	i__2 = lstpm1;
	for (i__ = 2; i__ <= i__2; ++i__) {
	    s2 = s2 + f[i__ + (k * f_dim2 + 1) * f_dim1] + f[i__ + (mstop + k 
		    * f_dim2) * f_dim1];
/* L163: */
	}
	s2 = (s2 + (f[(k * f_dim2 + 1) * f_dim1 + 1] + f[(mstop + k * f_dim2) 
		* f_dim1 + 1] + f[lstop + (k * f_dim2 + 1) * f_dim1] + f[
		lstop + (mstop + k * f_dim2) * f_dim1]) / xlp) / ylp;
	s1 += s2;
/* L164: */
    }
    s = (f[(f_dim2 + 1) * f_dim1 + 1] + f[lstop + (f_dim2 + 1) * f_dim1] + f[(
	    nstop * f_dim2 + 1) * f_dim1 + 1] + f[lstop + (nstop * f_dim2 + 1)
	     * f_dim1] + f[(mstop + f_dim2) * f_dim1 + 1] + f[lstop + (mstop 
	    + f_dim2) * f_dim1] + f[(mstop + nstop * f_dim2) * f_dim1 + 1] + 
	    f[lstop + (mstop + nstop * f_dim2) * f_dim1]) / (xlp * ylp);
    i__1 = mstpm1;
    for (j = 2; j <= i__1; ++j) {
	i__2 = lstpm1;
	for (i__ = 2; i__ <= i__2; ++i__) {
	    s = s + f[i__ + (j + f_dim2) * f_dim1] + f[i__ + (j + nstop * 
		    f_dim2) * f_dim1];
/* L165: */
	}
/* L166: */
    }
    s2 = 0.f;
    i__1 = lstpm1;
    for (i__ = 2; i__ <= i__1; ++i__) {
	s2 = s2 + f[i__ + (f_dim2 + 1) * f_dim1] + f[i__ + (nstop * f_dim2 + 
		1) * f_dim1] + f[i__ + (mstop + f_dim2) * f_dim1] + f[i__ + (
		mstop + nstop * f_dim2) * f_dim1];
/* L167: */
    }
    s = s2 / ylp + s;
    s2 = 0.f;
    i__1 = mstpm1;
    for (j = 2; j <= i__1; ++j) {
	s2 = s2 + f[(j + f_dim2) * f_dim1 + 1] + f[(j + nstop * f_dim2) * 
		f_dim1 + 1] + f[lstop + (j + f_dim2) * f_dim1] + f[lstop + (j 
		+ nstop * f_dim2) * f_dim1];
/* L168: */
    }
    s = s2 / xlp + s;
    *pertrb = (s / zlp + s1) / (((real) (lunk + 1) - xlp) * ((real) (munk + 1)
	     - ylp) * ((real) (nunk + 1) - zlp));
    i__1 = lunk;
    for (i__ = 1; i__ <= i__1; ++i__) {
	i__2 = munk;
	for (j = 1; j <= i__2; ++j) {
	    i__3 = nunk;
	    for (k = 1; k <= i__3; ++k) {
		f[i__ + (j + k * f_dim2) * f_dim1] -= *pertrb;
/* L169: */
	    }
/* L170: */
	}
/* L171: */
    }
L172:
    nperod = 0;
    if (*nbdcnd == 0) {
	goto L173;
    }
    nperod = 1;
    w[1] = 0.f;
    w[iww - 1] = 0.f;
L173:
    pois3d_(lbdcnd, &lunk, &c1, mbdcnd, &munk, &c2, &nperod, &nunk, &w[1], &w[
	    iwb], &w[iwc], ldimf, mdimf, &f[lstart + (mstart + nstart * 
	    f_dim2) * f_dim1], &ir, &w[iww]);

/*     fill in sides for periodic boundary conditions. */

    if (lp != 1) {
	goto L180;
    }
    if (mp != 1) {
	goto L175;
    }
    i__1 = nstop;
    for (k = nstart; k <= i__1; ++k) {
	f[(mp1 + k * f_dim2) * f_dim1 + 1] = f[(k * f_dim2 + 1) * f_dim1 + 1];
/* L174: */
    }
    mstop = mp1;
L175:
    if (np != 1) {
	goto L177;
    }
    i__1 = mstop;
    for (j = mstart; j <= i__1; ++j) {
	f[(j + np1 * f_dim2) * f_dim1 + 1] = f[(j + f_dim2) * f_dim1 + 1];
/* L176: */
    }
    nstop = np1;
L177:
    i__1 = mstop;
    for (j = mstart; j <= i__1; ++j) {
	i__2 = nstop;
	for (k = nstart; k <= i__2; ++k) {
	    f[lp1 + (j + k * f_dim2) * f_dim1] = f[(j + k * f_dim2) * f_dim1 
		    + 1];
/* L178: */
	}
/* L179: */
    }
L180:
    if (mp != 1) {
	goto L185;
    }
    if (np != 1) {
	goto L182;
    }
    i__1 = lstop;
    for (i__ = lstart; i__ <= i__1; ++i__) {
	f[i__ + (np1 * f_dim2 + 1) * f_dim1] = f[i__ + (f_dim2 + 1) * f_dim1];
/* L181: */
    }
    nstop = np1;
L182:
    i__1 = lstop;
    for (i__ = lstart; i__ <= i__1; ++i__) {
	i__2 = nstop;
	for (k = nstart; k <= i__2; ++k) {
	    f[i__ + (mp1 + k * f_dim2) * f_dim1] = f[i__ + (k * f_dim2 + 1) * 
		    f_dim1];
/* L183: */
	}
/* L184: */
    }
L185:
    if (np != 1) {
	goto L188;
    }
    i__1 = lstop;
    for (i__ = lstart; i__ <= i__1; ++i__) {
	i__2 = mstop;
	for (j = mstart; j <= i__2; ++j) {
	    f[i__ + (j + np1 * f_dim2) * f_dim1] = f[i__ + (j + f_dim2) * 
		    f_dim1];
/* L186: */
	}
/* L187: */
    }
L188:
    return 0;
} /* hw3crt_ */

/* Subroutine */ int pois3d_(integer *lperod, integer *l, real *c1, integer *
	mperod, integer *m, real *c2, integer *nperod, integer *n, real *a, 
	real *b, real *c__, integer *ldimf, integer *mdimf, real *f, integer *
	ierror, real *w)
{
    /* System generated locals */
    integer f_dim1, f_dim2, f_offset, i__1, i__2, i__3;

    /* Local variables */
    static integer iwbb, nodd, nhmk;
    static real save[6];
    static integer nhpk, i__, j, k;
    extern /* Subroutine */ int pos3d1_(integer *, integer *, integer *, 
	    integer *, integer *, real *, real *, real *, integer *, integer *
	    , real *, real *, real *, real *, real *, real *, real *, real *, 
	    real *, real *);
    static integer iwyrt, nh, lp, mp, np, iwd, iwt, iwx, iwy, nhm1;



/*     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*     *                                                               * */
/*     *                        f i s h p a k                          * */
/*     *                                                               * */
/*     *                                                               * */
/*     *     a package of fortran subprograms for the solution of      * */
/*     *                                                               * */
/*     *      separable elliptic partial differential equations        * */
/*     *                                                               * */
/*     *                  (version 3.1 , october 1980)                  * */
/*     *                                                               * */
/*     *                             by                                * */
/*     *                                                               * */
/*     *        john adams, paul swarztrauber and roland sweet         * */
/*     *                                                               * */
/*     *                             of                                * */
/*     *                                                               * */
/*     *         the national center for atmospheric research          * */
/*     *                                                               * */
/*     *                boulder, colorado  (80307)  u.s.a.             * */
/*     *                                                               * */
/*     *                   which is sponsored by                       * */
/*     *                                                               * */
/*     *              the national science foundation                  * */
/*     *                                                               * */
/*     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */


/*    * * * * * * * * *  purpose    * * * * * * * * * * * * * * * * * * */

/*     subroutine pois3d solves the linear system of equations */

/*       c1*(x(i-1,j,k)-2.*x(i,j,k)+x(i+1,j,k)) */
/*     + c2*(x(i,j-1,k)-2.*x(i,j,k)+x(i,j+1,k)) */
/*     + a(k)*x(i,j,k-1)+b(k)*x(i,j,k)+c(k)*x(i,j,k+1) = f(i,j,k) */

/*     for  i=1,2,...,l , j=1,2,...,m , and k=1,2,...,n . */

/*     the indices k-1 and k+1 are evaluated modulo n, i.e. */
/*     x(i,j,0) = x(i,j,n) and x(i,j,n+1) = x(i,j,1). the unknowns */
/*     x(0,j,k), x(l+1,j,k), x(i,0,k), and x(i,m+1,k) are assumed to take */
/*     on certain prescribed values described below. */

/*    * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */


/*    * * * * * * * *    parameter description     * * * * * * * * * * */


/*            * * * * * *   on input    * * * * * * */

/*     lperod   indicates the values that x(0,j,k) and x(l+1,j,k) are */
/*              assumed to have. */

/*              = 0  if x(0,j,k) = x(l,j,k) and x(l+1,j,k) = x(1,j,k). */
/*              = 1  if x(0,j,k) = x(l+1,j,k) = 0. */
/*              = 2  if x(0,j,k) = 0  and x(l+1,j,k) = x(l-1,j,k). */
/*              = 3  if x(0,j,k) = x(2,j,k) and x(l+1,j,k) = x(l-1,j,k). */
/*              = 4  if x(0,j,k) = x(2,j,k) and x(l+1,j,k) = 0. */

/*     l        the number of unknowns in the i-direction. l must be at */
/*              least 3. */

/*     c1       the real constant that appears in the above equation. */

/*     mperod   indicates the values that x(i,0,k) and x(i,m+1,k) are */
/*              assumed to have. */

/*              = 0  if x(i,0,k) = x(i,m,k) and x(i,m+1,k) = x(i,1,k). */
/*              = 1  if x(i,0,k) = x(i,m+1,k) = 0. */
/*              = 2  if x(i,0,k) = 0 and x(i,m+1,k) = x(i,m-1,k). */
/*              = 3  if x(i,0,k) = x(i,2,k) and x(i,m+1,k) = x(i,m-1,k). */
/*              = 4  if x(i,0,k) = x(i,2,k) and x(i,m+1,k) = 0. */

/*     m        the number of unknowns in the j-direction. m must be at */
/*              least 3. */

/*     c2       the real constant which appears in the above equation. */

/*     nperod   = 0  if a(1) and c(n) are not zero. */
/*              = 1  if a(1) = c(n) = 0. */

/*     n        the number of unknowns in the k-direction. n must be at */
/*              least 3. */


/*     a,b,c    one-dimensional arrays of length n that specify the */
/*              coefficients in the linear equations given above. */

/*              if nperod = 0 the array elements must not depend upon the */
/*              index k, but must be constant.  specifically,the */
/*              subroutine checks the following condition */

/*                          a(k) = c(1) */
/*                          c(k) = c(1) */
/*                          b(k) = b(1) */

/*                  for k=1,2,...,n. */

/*     ldimf    the row (or first) dimension of the three-dimensional */
/*              array f as it appears in the program calling pois3d. */
/*              this parameter is used to specify the variable dimension */
/*              of f.  ldimf must be at least l. */

/*     mdimf    the column (or second) dimension of the three-dimensional */
/*              array f as it appears in the program calling pois3d. */
/*              this parameter is used to specify the variable dimension */
/*              of f.  mdimf must be at least m. */

/*     f        a three-dimensional array that specifies the values of */
/*              the right side of the linear system of equations given */
/*              above.  f must be dimensioned at least l x m x n. */

/*     w        a one-dimensional array that must be provided by the */
/*              user for work space.  the length of w must be at least */
/*              30 + l + m + 2*n + max(l,m,n) + */
/*              7*(int((l+1)/2) + int((m+1)/2)). */


/*            * * * * * *   on output   * * * * * * */

/*     f        contains the solution x. */

/*     ierror   an error flag that indicates invalid input parameters. */
/*              except for number zero, a solution is not attempted. */
/*              = 0  no error */
/*              = 1  if lperod .lt. 0 or .gt. 4 */
/*              = 2  if l .lt. 3 */
/*              = 3  if mperod .lt. 0 or .gt. 4 */
/*              = 4  if m .lt. 3 */
/*              = 5  if nperod .lt. 0 or .gt. 1 */
/*              = 6  if n .lt. 3 */
/*              = 7  if ldimf .lt. l */
/*              = 8  if mdimf .lt. m */
/*              = 9  if a(k) .ne. c(1) or c(k) .ne. c(1) or b(i) .ne.b(1) */
/*                      for some k=1,2,...,n. */
/*              = 10 if nperod = 1 and a(1) .ne. 0 or c(n) .ne. 0 */

/*              since this is the only means of indicating a possibly */
/*              incorrect call to pois3d, the user should test ierror */
/*              after the call. */


/*    * * * * * * *   program specifications    * * * * * * * * * * * * */

/*     dimension of   a(n),b(n),c(n),f(ldimf,mdimf,n), */
/*     arguments      w(see argument list) */

/*     latest         december 1, 1978 */
/*     revision */

/*     subprograms    pois3d,pos3d1,trid,rffti,rfftf,rfftf1,rfftb, */
/*     required       rfftb1,costi,cost,sinti,sint,cosqi,cosqf,cosqf1 */
/*                    cosqb,cosqb1,sinqi,sinqf,sinqb,cffti,cffti1, */
/*                    cfftb,cfftb1,passb2,passb3,passb4,passb,cfftf, */
/*                    cfftf1,passf1,passf2,passf3,passf4,passf,pimach, */

/*     special        none */
/*     conditions */

/*     common         value */
/*     blocks */

/*     i/o            none */

/*     precision      single */

/*     specialist     roland sweet */

/*     language       fortran */

/*     history        written by roland sweet at ncar in july,1977 */

/*     algorithm      this subroutine solves three-dimensional block */
/*                    tridiagonal linear systems arising from finite */
/*                    difference approximations to three-dimensional */
/*                    poisson equations using the fourier transform */
/*                    package sclrfftpak written by paul swarztrauber. */

/*     space          6561(decimal) = 14641(octal) locations on the */
/*     required       ncar control data 7600 */

/*     timing and        the execution time t on the ncar control data */
/*     accuracy       7600 for subroutine pois3d is roughly proportional */
/*                    to l*m*n*(log2(l)+log2(m)+5), but also depends on */
/*                    input parameters lperod and mperod.  some typical */
/*                    values are listed in the table below when nperod=0. */
/*                       to measure the accuracy of the algorithm a */
/*                    uniform random number generator was used to create */
/*                    a solution array x for the system given in the */
/*                    'purpose' with */

/*                       a(k) = c(k) = -0.5*b(k) = 1,       k=1,2,...,n */

/*                    and, when nperod = 1 */

/*                       a(1) = c(n) = 0 */
/*                       a(n) = c(1) = 2. */

/*                    the solution x was substituted into the given sys- */
/*                    tem and, using double precision, a right side y was */
/*                    computed.  using this array y subroutine pois was */
/*                    called to produce an approximate solution z.  then */
/*                    the relative error, defined as */

/*                    e = max(abs(z(i,j,k)-x(i,j,k)))/max(abs(x(i,j,k))) */

/*                    where the two maxima are taken over i=1,2,...,l, */
/*                    j=1,2,...,m and k=1,2,...,n, was computed.  the */
/*                    value of e is given in the table below for some */
/*                    typical values of l,m and n. */


/*                       l(=m=n)   lperod    mperod    t(msecs)    e */
/*                       ------    ------    ------    --------  ------ */

/*                         16        0         0         272     1.e-13 */
/*                         15        1         1         287     4.e-13 */
/*                         17        3         3         338     2.e-13 */
/*                         32        0         0        1755     2.e-13 */
/*                         31        1         1        1894     2.e-12 */
/*                         33        3         3        2042     7.e-13 */


/*     portability    american national standards institute fortran. */
/*                    the machine dependent constant pi is defined in */
/*                    function pimach. */

/*     required       cos,sin,atan */
/*     resident */
/*     routines */

/*     reference      none */

/*    * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

    /* Parameter adjustments */
    --a;
    --b;
    --c__;
    f_dim1 = *ldimf;
    f_dim2 = *mdimf;
    f_offset = 1 + f_dim1 * (1 + f_dim2 * 1);
    f -= f_offset;
    --w;

    /* Function Body */
    lp = *lperod + 1;
    mp = *mperod + 1;
    np = *nperod + 1;

/*     check for invalid input. */

    *ierror = 0;
    if (lp < 1 || lp > 5) {
	*ierror = 1;
    }
    if (*l < 3) {
	*ierror = 2;
    }
    if (mp < 1 || mp > 5) {
	*ierror = 3;
    }
    if (*m < 3) {
	*ierror = 4;
    }
    if (np < 1 || np > 2) {
	*ierror = 5;
    }
    if (*n < 3) {
	*ierror = 6;
    }
    if (*ldimf < *l) {
	*ierror = 7;
    }
    if (*mdimf < *m) {
	*ierror = 8;
    }
    if (np != 1) {
	goto L103;
    }
    i__1 = *n;
    for (k = 1; k <= i__1; ++k) {
	if (a[k] != c__[1]) {
	    goto L102;
	}
	if (c__[k] != c__[1]) {
	    goto L102;
	}
	if (b[k] != b[1]) {
	    goto L102;
	}
/* L101: */
    }
    goto L104;
L102:
    *ierror = 9;
L103:
    if (*nperod == 1 && (a[1] != 0.f || c__[*n] != 0.f)) {
	*ierror = 10;
    }
L104:
    if (*ierror != 0) {
	goto L122;
    }
    iwyrt = *l + 1;
    iwt = iwyrt + *m;
/* Computing MAX */
    i__1 = max(*l,*m);
    iwd = iwt + max(i__1,*n) + 1;
    iwbb = iwd + *n;
    iwx = iwbb + *n;
    iwy = iwx + (*l + 1) / 2 * 7 + 15;
    switch (np) {
	case 1:  goto L105;
	case 2:  goto L114;
    }

/*     reorder unknowns when nperod = 0. */

L105:
    nh = (*n + 1) / 2;
    nhm1 = nh - 1;
    nodd = 1;
    if (nh << 1 == *n) {
	nodd = 2;
    }
    i__1 = *l;
    for (i__ = 1; i__ <= i__1; ++i__) {
	i__2 = *m;
	for (j = 1; j <= i__2; ++j) {
	    i__3 = nhm1;
	    for (k = 1; k <= i__3; ++k) {
		nhpk = nh + k;
		nhmk = nh - k;
		w[k] = f[i__ + (j + nhmk * f_dim2) * f_dim1] - f[i__ + (j + 
			nhpk * f_dim2) * f_dim1];
		w[nhpk] = f[i__ + (j + nhmk * f_dim2) * f_dim1] + f[i__ + (j 
			+ nhpk * f_dim2) * f_dim1];
/* L106: */
	    }
	    w[nh] = f[i__ + (j + nh * f_dim2) * f_dim1] * 2.f;
	    switch (nodd) {
		case 1:  goto L108;
		case 2:  goto L107;
	    }
L107:
	    w[*n] = f[i__ + (j + *n * f_dim2) * f_dim1] * 2.f;
L108:
	    i__3 = *n;
	    for (k = 1; k <= i__3; ++k) {
		f[i__ + (j + k * f_dim2) * f_dim1] = w[k];
/* L109: */
	    }
/* L110: */
	}
/* L111: */
    }
    save[0] = c__[nhm1];
    save[1] = a[nh];
    save[2] = c__[nh];
    save[3] = b[nhm1];
    save[4] = b[*n];
    save[5] = a[*n];
    c__[nhm1] = 0.f;
    a[nh] = 0.f;
    c__[nh] *= 2.f;
    switch (nodd) {
	case 1:  goto L112;
	case 2:  goto L113;
    }
L112:
    b[nhm1] -= a[nh - 1];
    b[*n] += a[*n];
    goto L114;
L113:
    a[*n] = c__[nh];
L114:
    pos3d1_(&lp, l, &mp, m, n, &a[1], &b[1], &c__[1], ldimf, mdimf, &f[
	    f_offset], &w[1], &w[iwyrt], &w[iwt], &w[iwd], &w[iwx], &w[iwy], 
	    c1, c2, &w[iwbb]);
    switch (np) {
	case 1:  goto L115;
	case 2:  goto L122;
    }
L115:
    i__1 = *l;
    for (i__ = 1; i__ <= i__1; ++i__) {
	i__2 = *m;
	for (j = 1; j <= i__2; ++j) {
	    i__3 = nhm1;
	    for (k = 1; k <= i__3; ++k) {
		nhmk = nh - k;
		nhpk = nh + k;
		w[nhmk] = (f[i__ + (j + nhpk * f_dim2) * f_dim1] + f[i__ + (j 
			+ k * f_dim2) * f_dim1]) * .5f;
		w[nhpk] = (f[i__ + (j + nhpk * f_dim2) * f_dim1] - f[i__ + (j 
			+ k * f_dim2) * f_dim1]) * .5f;
/* L116: */
	    }
	    w[nh] = f[i__ + (j + nh * f_dim2) * f_dim1] * .5f;
	    switch (nodd) {
		case 1:  goto L118;
		case 2:  goto L117;
	    }
L117:
	    w[*n] = f[i__ + (j + *n * f_dim2) * f_dim1] * .5f;
L118:
	    i__3 = *n;
	    for (k = 1; k <= i__3; ++k) {
		f[i__ + (j + k * f_dim2) * f_dim1] = w[k];
/* L119: */
	    }
/* L120: */
	}
/* L121: */
    }
    c__[nhm1] = save[0];
    a[nh] = save[1];
    c__[nh] = save[2];
    b[nhm1] = save[3];
    b[*n] = save[4];
    a[*n] = save[5];
L122:
    return 0;
} /* pois3d_ */

/* Subroutine */ int pos3d1_(integer *lp, integer *l, integer *mp, integer *m,
	 integer *n, real *a, real *b, real *c__, integer *ldimf, integer *
	mdimf, real *f, real *xrt, real *yrt, real *t, real *d__, real *wx, 
	real *wy, real *c1, real *c2, real *bb)
{
    /* System generated locals */
    integer f_dim1, f_dim2, f_offset, i__1, i__2, i__3;
    real r__1;

    /* Builtin functions */
    double sin(doublereal);

    /* Local variables */
    extern /* Subroutine */ int trid_(integer *, real *, real *, real *, real 
	    *, real *), cost_(integer *, real *, real *), sint_(integer *, 
	    real *, real *);
    static integer i__, j, k, lrdel, mrdel;
    extern /* Subroutine */ int rfftb_(integer *, real *, real *), rfftf_(
	    integer *, real *, real *), cosqb_(integer *, real *, real *);
    static real scalx;
    extern /* Subroutine */ int rffti_(integer *, real *);
    static real scaly;
    static integer ifwrd;
    extern /* Subroutine */ int cosqi_(integer *, real *), sinqb_(integer *, 
	    real *, real *), sinqf_(integer *, real *, real *), costi_(
	    integer *, real *), cosqf_(integer *, real *, real *), sinqi_(
	    integer *, real *), sinti_(integer *, real *);
    static real di, dj, pi;
    static integer lr, mr, nr;
    extern doublereal pimach_(real *);
    static real dx, dy;
    static integer is;
    static real dum;

    /* Parameter adjustments */
    --a;
    --b;
    --c__;
    f_dim1 = *ldimf;
    f_dim2 = *mdimf;
    f_offset = 1 + f_dim1 * (1 + f_dim2 * 1);
    f -= f_offset;
    --xrt;
    --yrt;
    --t;
    --d__;
    --wx;
    --wy;
    --bb;

    /* Function Body */
    pi = pimach_(&dum);
    lr = *l;
    mr = *m;
    nr = *n;

/*     generate transform roots */

    lrdel = (*lp - 1) * (*lp - 3) * (*lp - 5) / 3;
    scalx = (real) (lr + lrdel);
    dx = pi / (scalx * 2.f);
    switch (*lp) {
	case 1:  goto L108;
	case 2:  goto L103;
	case 3:  goto L101;
	case 4:  goto L102;
	case 5:  goto L101;
    }
L101:
    di = .5f;
    scalx *= 2.f;
    goto L104;
L102:
    di = 1.f;
    goto L104;
L103:
    di = 0.f;
L104:
    i__1 = lr;
    for (i__ = 1; i__ <= i__1; ++i__) {
/* Computing 2nd power */
	r__1 = sin(((real) i__ - di) * dx);
	xrt[i__] = *c1 * -4.f * (r__1 * r__1);
/* L105: */
    }
    scalx *= 2.f;
    switch (*lp) {
	case 1:  goto L112;
	case 2:  goto L106;
	case 3:  goto L110;
	case 4:  goto L107;
	case 5:  goto L111;
    }
L106:
    sinti_(&lr, &wx[1]);
    goto L112;
L107:
    costi_(&lr, &wx[1]);
    goto L112;
L108:
    xrt[1] = 0.f;
    xrt[lr] = *c1 * -4.f;
    i__1 = lr;
    for (i__ = 3; i__ <= i__1; i__ += 2) {
/* Computing 2nd power */
	r__1 = sin((real) (i__ - 1) * dx);
	xrt[i__ - 1] = *c1 * -4.f * (r__1 * r__1);
	xrt[i__] = xrt[i__ - 1];
/* L109: */
    }
    rffti_(&lr, &wx[1]);
    goto L112;
L110:
    sinqi_(&lr, &wx[1]);
    goto L112;
L111:
    cosqi_(&lr, &wx[1]);
L112:
    mrdel = (*mp - 1) * (*mp - 3) * (*mp - 5) / 3;
    scaly = (real) (mr + mrdel);
    dy = pi / (scaly * 2.f);
    switch (*mp) {
	case 1:  goto L120;
	case 2:  goto L115;
	case 3:  goto L113;
	case 4:  goto L114;
	case 5:  goto L113;
    }
L113:
    dj = .5f;
    scaly *= 2.f;
    goto L116;
L114:
    dj = 1.f;
    goto L116;
L115:
    dj = 0.f;
L116:
    i__1 = mr;
    for (j = 1; j <= i__1; ++j) {
/* Computing 2nd power */
	r__1 = sin(((real) j - dj) * dy);
	yrt[j] = *c2 * -4.f * (r__1 * r__1);
/* L117: */
    }
    scaly *= 2.f;
    switch (*mp) {
	case 1:  goto L124;
	case 2:  goto L118;
	case 3:  goto L122;
	case 4:  goto L119;
	case 5:  goto L123;
    }
L118:
    sinti_(&mr, &wy[1]);
    goto L124;
L119:
    costi_(&mr, &wy[1]);
    goto L124;
L120:
    yrt[1] = 0.f;
    yrt[mr] = *c2 * -4.f;
    i__1 = mr;
    for (j = 3; j <= i__1; j += 2) {
/* Computing 2nd power */
	r__1 = sin((real) (j - 1) * dy);
	yrt[j - 1] = *c2 * -4.f * (r__1 * r__1);
	yrt[j] = yrt[j - 1];
/* L121: */
    }
    rffti_(&mr, &wy[1]);
    goto L124;
L122:
    sinqi_(&mr, &wy[1]);
    goto L124;
L123:
    cosqi_(&mr, &wy[1]);
L124:
    ifwrd = 1;
    is = 1;
L125:

/*     transform x */

    i__1 = mr;
    for (j = 1; j <= i__1; ++j) {
	i__2 = nr;
	for (k = 1; k <= i__2; ++k) {
	    i__3 = lr;
	    for (i__ = 1; i__ <= i__3; ++i__) {
		t[i__] = f[i__ + (j + k * f_dim2) * f_dim1];
/* L126: */
	    }
	    switch (*lp) {
		case 1:  goto L127;
		case 2:  goto L130;
		case 3:  goto L131;
		case 4:  goto L134;
		case 5:  goto L135;
	    }
L127:
	    switch (ifwrd) {
		case 1:  goto L128;
		case 2:  goto L129;
	    }
L128:
	    rfftf_(&lr, &t[1], &wx[1]);
	    goto L138;
L129:
	    rfftb_(&lr, &t[1], &wx[1]);
	    goto L138;
L130:
	    sint_(&lr, &t[1], &wx[1]);
	    goto L138;
L131:
	    switch (ifwrd) {
		case 1:  goto L132;
		case 2:  goto L133;
	    }
L132:
	    sinqf_(&lr, &t[1], &wx[1]);
	    goto L138;
L133:
	    sinqb_(&lr, &t[1], &wx[1]);
	    goto L138;
L134:
	    cost_(&lr, &t[1], &wx[1]);
	    goto L138;
L135:
	    switch (ifwrd) {
		case 1:  goto L136;
		case 2:  goto L137;
	    }
L136:
	    cosqf_(&lr, &t[1], &wx[1]);
	    goto L138;
L137:
	    cosqb_(&lr, &t[1], &wx[1]);
L138:
	    i__3 = lr;
	    for (i__ = 1; i__ <= i__3; ++i__) {
		f[i__ + (j + k * f_dim2) * f_dim1] = t[i__];
/* L139: */
	    }
/* L140: */
	}
/* L141: */
    }
    switch (ifwrd) {
	case 1:  goto L142;
	case 2:  goto L164;
    }

/*     transform y */

L142:
    i__1 = lr;
    for (i__ = 1; i__ <= i__1; ++i__) {
	i__2 = nr;
	for (k = 1; k <= i__2; ++k) {
	    i__3 = mr;
	    for (j = 1; j <= i__3; ++j) {
		t[j] = f[i__ + (j + k * f_dim2) * f_dim1];
/* L143: */
	    }
	    switch (*mp) {
		case 1:  goto L144;
		case 2:  goto L147;
		case 3:  goto L148;
		case 4:  goto L151;
		case 5:  goto L152;
	    }
L144:
	    switch (ifwrd) {
		case 1:  goto L145;
		case 2:  goto L146;
	    }
L145:
	    rfftf_(&mr, &t[1], &wy[1]);
	    goto L155;
L146:
	    rfftb_(&mr, &t[1], &wy[1]);
	    goto L155;
L147:
	    sint_(&mr, &t[1], &wy[1]);
	    goto L155;
L148:
	    switch (ifwrd) {
		case 1:  goto L149;
		case 2:  goto L150;
	    }
L149:
	    sinqf_(&mr, &t[1], &wy[1]);
	    goto L155;
L150:
	    sinqb_(&mr, &t[1], &wy[1]);
	    goto L155;
L151:
	    cost_(&mr, &t[1], &wy[1]);
	    goto L155;
L152:
	    switch (ifwrd) {
		case 1:  goto L153;
		case 2:  goto L154;
	    }
L153:
	    cosqf_(&mr, &t[1], &wy[1]);
	    goto L155;
L154:
	    cosqb_(&mr, &t[1], &wy[1]);
L155:
	    i__3 = mr;
	    for (j = 1; j <= i__3; ++j) {
		f[i__ + (j + k * f_dim2) * f_dim1] = t[j];
/* L156: */
	    }
/* L157: */
	}
/* L158: */
    }
    switch (ifwrd) {
	case 1:  goto L159;
	case 2:  goto L125;
    }
L159:

/*     solve tridiagonal systems in z */

    i__1 = lr;
    for (i__ = 1; i__ <= i__1; ++i__) {
	i__2 = mr;
	for (j = 1; j <= i__2; ++j) {
	    i__3 = nr;
	    for (k = 1; k <= i__3; ++k) {
		bb[k] = b[k] + xrt[i__] + yrt[j];
		t[k] = f[i__ + (j + k * f_dim2) * f_dim1];
/* L160: */
	    }
	    trid_(&nr, &a[1], &bb[1], &c__[1], &t[1], &d__[1]);
	    i__3 = nr;
	    for (k = 1; k <= i__3; ++k) {
		f[i__ + (j + k * f_dim2) * f_dim1] = t[k];
/* L161: */
	    }
/* L162: */
	}
/* L163: */
    }
    ifwrd = 2;
    is = -1;
    goto L142;
L164:
    i__1 = lr;
    for (i__ = 1; i__ <= i__1; ++i__) {
	i__2 = mr;
	for (j = 1; j <= i__2; ++j) {
	    i__3 = nr;
	    for (k = 1; k <= i__3; ++k) {
		f[i__ + (j + k * f_dim2) * f_dim1] /= scalx * scaly;
/* L165: */
	    }
/* L166: */
	}
/* L167: */
    }
    return 0;
} /* pos3d1_ */

/* Subroutine */ int trid_(integer *mr, real *a, real *b, real *c__, real *y, 
	real *d__)
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    static integer i__, m;
    static real z__;
    static integer ip, mm1;

    /* Parameter adjustments */
    --d__;
    --y;
    --c__;
    --b;
    --a;

    /* Function Body */
    m = *mr;
    mm1 = m - 1;
    z__ = 1.f / b[1];
    d__[1] = c__[1] * z__;
    y[1] *= z__;
    i__1 = mm1;
    for (i__ = 2; i__ <= i__1; ++i__) {
	z__ = 1.f / (b[i__] - a[i__] * d__[i__ - 1]);
	d__[i__] = c__[i__] * z__;
	y[i__] = (y[i__] - a[i__] * y[i__ - 1]) * z__;
/* L101: */
    }
    z__ = b[m] - a[m] * d__[mm1];
    if (z__ != 0.f) {
	goto L102;
    }
    y[m] = 0.f;
    goto L103;
L102:
    y[m] = (y[m] - a[m] * y[mm1]) / z__;
L103:
    i__1 = mm1;
    for (ip = 1; ip <= i__1; ++ip) {
	i__ = m - ip;
	y[i__] -= d__[i__] * y[i__ + 1];
/* L104: */
    }
    return 0;
} /* trid_ */

/* Subroutine */ int cosqb_(integer *n, real *x, real *wsave)
{
    /* Initialized data */

    static real tsqrt2 = 2.82842712474619f;

    /* System generated locals */
    integer i__1;

    /* Local variables */
    static real x1;
    extern /* Subroutine */ int cosqb1_(integer *, real *, real *, real *);

    /* Parameter adjustments */
    --wsave;
    --x;

    /* Function Body */
    if ((i__1 = *n - 2) < 0) {
	goto L101;
    } else if (i__1 == 0) {
	goto L102;
    } else {
	goto L103;
    }
L101:
    x[1] *= 4.f;
    return 0;
L102:
    x1 = (x[1] + x[2]) * 4.f;
    x[2] = tsqrt2 * (x[1] - x[2]);
    x[1] = x1;
    return 0;
L103:
    cosqb1_(n, &x[1], &wsave[1], &wsave[*n + 1]);
    return 0;
} /* cosqb_ */

/* Subroutine */ int cosqb1_(integer *n, real *x, real *w, real *xh)
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    static integer modn, i__, k;
    extern /* Subroutine */ int rfftb_(integer *, real *, real *);
    static integer kc, np2, ns2;
    static real xim1;

    /* Parameter adjustments */
    --xh;
    --w;
    --x;

    /* Function Body */
    ns2 = (*n + 1) / 2;
    np2 = *n + 2;
    i__1 = *n;
    for (i__ = 3; i__ <= i__1; i__ += 2) {
	xim1 = x[i__ - 1] + x[i__];
	x[i__] -= x[i__ - 1];
	x[i__ - 1] = xim1;
/* L101: */
    }
    x[1] += x[1];
    modn = *n % 2;
    if (modn == 0) {
	x[*n] += x[*n];
    }
    rfftb_(n, &x[1], &xh[1]);
    i__1 = ns2;
    for (k = 2; k <= i__1; ++k) {
	kc = np2 - k;
	xh[k] = w[k - 1] * x[kc] + w[kc - 1] * x[k];
	xh[kc] = w[k - 1] * x[k] - w[kc - 1] * x[kc];
/* L102: */
    }
    if (modn == 0) {
	x[ns2 + 1] = w[ns2] * (x[ns2 + 1] + x[ns2 + 1]);
    }
    i__1 = ns2;
    for (k = 2; k <= i__1; ++k) {
	kc = np2 - k;
	x[k] = xh[k] + xh[kc];
	x[kc] = xh[k] - xh[kc];
/* L103: */
    }
    x[1] += x[1];
    return 0;
} /* cosqb1_ */

/* Subroutine */ int cosqf_(integer *n, real *x, real *wsave)
{
    /* Initialized data */

    static real sqrt2 = 1.4142135623731f;

    /* System generated locals */
    integer i__1;

    /* Local variables */
    static real tsqx;
    extern /* Subroutine */ int cosqf1_(integer *, real *, real *, real *);

    /* Parameter adjustments */
    --wsave;
    --x;

    /* Function Body */
    if ((i__1 = *n - 2) < 0) {
	goto L102;
    } else if (i__1 == 0) {
	goto L101;
    } else {
	goto L103;
    }
L101:
    tsqx = sqrt2 * x[2];
    x[2] = x[1] - tsqx;
    x[1] += tsqx;
L102:
    return 0;
L103:
    cosqf1_(n, &x[1], &wsave[1], &wsave[*n + 1]);
    return 0;
} /* cosqf_ */

/* Subroutine */ int cosqf1_(integer *n, real *x, real *w, real *xh)
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    static integer modn, i__, k;
    extern /* Subroutine */ int rfftf_(integer *, real *, real *);
    static integer kc, np2, ns2;
    static real xim1;

    /* Parameter adjustments */
    --xh;
    --w;
    --x;

    /* Function Body */
    ns2 = (*n + 1) / 2;
    np2 = *n + 2;
    i__1 = ns2;
    for (k = 2; k <= i__1; ++k) {
	kc = np2 - k;
	xh[k] = x[k] + x[kc];
	xh[kc] = x[k] - x[kc];
/* L101: */
    }
    modn = *n % 2;
    if (modn == 0) {
	xh[ns2 + 1] = x[ns2 + 1] + x[ns2 + 1];
    }
    i__1 = ns2;
    for (k = 2; k <= i__1; ++k) {
	kc = np2 - k;
	x[k] = w[k - 1] * xh[kc] + w[kc - 1] * xh[k];
	x[kc] = w[k - 1] * xh[k] - w[kc - 1] * xh[kc];
/* L102: */
    }
    if (modn == 0) {
	x[ns2 + 1] = w[ns2] * xh[ns2 + 1];
    }
    rfftf_(n, &x[1], &xh[1]);
    i__1 = *n;
    for (i__ = 3; i__ <= i__1; i__ += 2) {
	xim1 = x[i__ - 1] - x[i__];
	x[i__] = x[i__ - 1] + x[i__];
	x[i__ - 1] = xim1;
/* L103: */
    }
    return 0;
} /* cosqf1_ */

/* Subroutine */ int cosqi_(integer *n, real *wsave)
{
    /* Initialized data */

    static real pih = 1.57079632679491f;

    /* System generated locals */
    integer i__1;

    /* Builtin functions */
    double cos(doublereal);

    /* Local variables */
    static integer k;
    extern /* Subroutine */ int rffti_(integer *, real *);
    static real fk, dt;

    /* Parameter adjustments */
    --wsave;

    /* Function Body */
    dt = pih / (real) (*n);
    fk = 0.f;
    i__1 = *n;
    for (k = 1; k <= i__1; ++k) {
	fk += 1.f;
	wsave[k] = cos(fk * dt);
/* L101: */
    }
    rffti_(n, &wsave[*n + 1]);
    return 0;
} /* cosqi_ */

/* Subroutine */ int cost_(integer *n, real *x, real *wsave)
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    static integer modn, i__, k;
    extern /* Subroutine */ int rfftf_(integer *, real *, real *);
    static real c1, t1, t2;
    static integer kc;
    static real xi;
    static integer nm1, np1;
    static real x1h;
    static integer ns2;
    static real tx2, x1p3, xim2;

    /* Parameter adjustments */
    --wsave;
    --x;

    /* Function Body */
    nm1 = *n - 1;
    np1 = *n + 1;
    ns2 = *n / 2;
    if ((i__1 = *n - 2) < 0) {
	goto L106;
    } else if (i__1 == 0) {
	goto L101;
    } else {
	goto L102;
    }
L101:
    x1h = x[1] + x[2];
    x[2] = x[1] - x[2];
    x[1] = x1h;
    return 0;
L102:
    if (*n > 3) {
	goto L103;
    }
    x1p3 = x[1] + x[3];
    tx2 = x[2] + x[2];
    x[2] = x[1] - x[3];
    x[1] = x1p3 + tx2;
    x[3] = x1p3 - tx2;
    return 0;
L103:
    c1 = x[1] - x[*n];
    x[1] += x[*n];
    i__1 = ns2;
    for (k = 2; k <= i__1; ++k) {
	kc = np1 - k;
	t1 = x[k] + x[kc];
	t2 = x[k] - x[kc];
	c1 += wsave[kc] * t2;
	t2 = wsave[k] * t2;
	x[k] = t1 - t2;
	x[kc] = t1 + t2;
/* L104: */
    }
    modn = *n % 2;
    if (modn != 0) {
	x[ns2 + 1] += x[ns2 + 1];
    }
    rfftf_(&nm1, &x[1], &wsave[*n + 1]);
    xim2 = x[2];
    x[2] = c1;
    i__1 = *n;
    for (i__ = 4; i__ <= i__1; i__ += 2) {
	xi = x[i__];
	x[i__] = x[i__ - 2] - x[i__ - 1];
	x[i__ - 1] = xim2;
	xim2 = xi;
/* L105: */
    }
    if (modn != 0) {
	x[*n] = xim2;
    }
L106:
    return 0;
} /* cost_ */

/* Subroutine */ int costi_(integer *n, real *wsave)
{
    /* Initialized data */

    static real pi = 3.14159265358979f;

    /* System generated locals */
    integer i__1;

    /* Builtin functions */
    double sin(doublereal), cos(doublereal);

    /* Local variables */
    static integer k;
    extern /* Subroutine */ int rffti_(integer *, real *);
    static integer kc;
    static real fk, dt;
    static integer nm1, np1, ns2;

    /* Parameter adjustments */
    --wsave;

    /* Function Body */
    if (*n <= 3) {
	return 0;
    }
    nm1 = *n - 1;
    np1 = *n + 1;
    ns2 = *n / 2;
    dt = pi / (real) nm1;
    fk = 0.f;
    i__1 = ns2;
    for (k = 2; k <= i__1; ++k) {
	kc = np1 - k;
	fk += 1.f;
	wsave[k] = sin(fk * dt) * 2.f;
	wsave[kc] = cos(fk * dt) * 2.f;
/* L101: */
    }
    rffti_(&nm1, &wsave[*n + 1]);
    return 0;
} /* costi_ */


/* Subroutine */ int radb2_(integer *ido, integer *l1, real *cc, real *ch, 
	real *wa1)
{
    /* System generated locals */
    integer cc_dim1, cc_offset, ch_dim1, ch_dim2, ch_offset, i__1, i__2;

    /* Local variables */
    static integer i__, k, ic;
    static real ti2, tr2;
    static integer idp2;

    /* Parameter adjustments */
    ch_dim1 = *ido;
    ch_dim2 = *l1;
    ch_offset = 1 + ch_dim1 * (1 + ch_dim2 * 1);
    ch -= ch_offset;
    cc_dim1 = *ido;
    cc_offset = 1 + cc_dim1 * 3;
    cc -= cc_offset;
    --wa1;

    /* Function Body */
    i__1 = *l1;
    for (k = 1; k <= i__1; ++k) {
	ch[(k + ch_dim2) * ch_dim1 + 1] = cc[((k << 1) + 1) * cc_dim1 + 1] + 
		cc[*ido + ((k << 1) + 2) * cc_dim1];
	ch[(k + (ch_dim2 << 1)) * ch_dim1 + 1] = cc[((k << 1) + 1) * cc_dim1 
		+ 1] - cc[*ido + ((k << 1) + 2) * cc_dim1];
/* L101: */
    }
    if ((i__1 = *ido - 2) < 0) {
	goto L107;
    } else if (i__1 == 0) {
	goto L105;
    } else {
	goto L102;
    }
L102:
    idp2 = *ido + 2;
    i__1 = *l1;
    for (k = 1; k <= i__1; ++k) {
	i__2 = *ido;
	for (i__ = 3; i__ <= i__2; i__ += 2) {
	    ic = idp2 - i__;
	    ch[i__ - 1 + (k + ch_dim2) * ch_dim1] = cc[i__ - 1 + ((k << 1) + 
		    1) * cc_dim1] + cc[ic - 1 + ((k << 1) + 2) * cc_dim1];
	    tr2 = cc[i__ - 1 + ((k << 1) + 1) * cc_dim1] - cc[ic - 1 + ((k << 
		    1) + 2) * cc_dim1];
	    ch[i__ + (k + ch_dim2) * ch_dim1] = cc[i__ + ((k << 1) + 1) * 
		    cc_dim1] - cc[ic + ((k << 1) + 2) * cc_dim1];
	    ti2 = cc[i__ + ((k << 1) + 1) * cc_dim1] + cc[ic + ((k << 1) + 2) 
		    * cc_dim1];
	    ch[i__ - 1 + (k + (ch_dim2 << 1)) * ch_dim1] = wa1[i__ - 2] * tr2 
		    - wa1[i__ - 1] * ti2;
	    ch[i__ + (k + (ch_dim2 << 1)) * ch_dim1] = wa1[i__ - 2] * ti2 + 
		    wa1[i__ - 1] * tr2;
/* L103: */
	}
/* L104: */
    }
    if (*ido % 2 == 1) {
	return 0;
    }
L105:
    i__1 = *l1;
    for (k = 1; k <= i__1; ++k) {
	ch[*ido + (k + ch_dim2) * ch_dim1] = cc[*ido + ((k << 1) + 1) * 
		cc_dim1] + cc[*ido + ((k << 1) + 1) * cc_dim1];
	ch[*ido + (k + (ch_dim2 << 1)) * ch_dim1] = -(cc[((k << 1) + 2) * 
		cc_dim1 + 1] + cc[((k << 1) + 2) * cc_dim1 + 1]);
/* L106: */
    }
L107:
    return 0;
} /* radb2_ */

/* Subroutine */ int radb3_(integer *ido, integer *l1, real *cc, real *ch, 
	real *wa1, real *wa2)
{
    /* Initialized data */

    static real taur = -.5f;
    static real taui = .866025403784439f;

    /* System generated locals */
    integer cc_dim1, cc_offset, ch_dim1, ch_dim2, ch_offset, i__1, i__2;

    /* Local variables */
    static integer i__, k, ic;
    static real ci2, ci3, di2, di3, cr2, cr3, dr2, dr3, ti2, tr2;
    static integer idp2;

    /* Parameter adjustments */
    ch_dim1 = *ido;
    ch_dim2 = *l1;
    ch_offset = 1 + ch_dim1 * (1 + ch_dim2 * 1);
    ch -= ch_offset;
    cc_dim1 = *ido;
    cc_offset = 1 + (cc_dim1 << 2);
    cc -= cc_offset;
    --wa1;
    --wa2;

    /* Function Body */
    i__1 = *l1;
    for (k = 1; k <= i__1; ++k) {
	tr2 = cc[*ido + (k * 3 + 2) * cc_dim1] + cc[*ido + (k * 3 + 2) * 
		cc_dim1];
	cr2 = cc[(k * 3 + 1) * cc_dim1 + 1] + taur * tr2;
	ch[(k + ch_dim2) * ch_dim1 + 1] = cc[(k * 3 + 1) * cc_dim1 + 1] + tr2;
	ci3 = taui * (cc[(k * 3 + 3) * cc_dim1 + 1] + cc[(k * 3 + 3) * 
		cc_dim1 + 1]);
	ch[(k + (ch_dim2 << 1)) * ch_dim1 + 1] = cr2 - ci3;
	ch[(k + ch_dim2 * 3) * ch_dim1 + 1] = cr2 + ci3;
/* L101: */
    }
    if (*ido == 1) {
	return 0;
    }
    idp2 = *ido + 2;
    i__1 = *l1;
    for (k = 1; k <= i__1; ++k) {
	i__2 = *ido;
	for (i__ = 3; i__ <= i__2; i__ += 2) {
	    ic = idp2 - i__;
	    tr2 = cc[i__ - 1 + (k * 3 + 3) * cc_dim1] + cc[ic - 1 + (k * 3 + 
		    2) * cc_dim1];
	    cr2 = cc[i__ - 1 + (k * 3 + 1) * cc_dim1] + taur * tr2;
	    ch[i__ - 1 + (k + ch_dim2) * ch_dim1] = cc[i__ - 1 + (k * 3 + 1) *
		     cc_dim1] + tr2;
	    ti2 = cc[i__ + (k * 3 + 3) * cc_dim1] - cc[ic + (k * 3 + 2) * 
		    cc_dim1];
	    ci2 = cc[i__ + (k * 3 + 1) * cc_dim1] + taur * ti2;
	    ch[i__ + (k + ch_dim2) * ch_dim1] = cc[i__ + (k * 3 + 1) * 
		    cc_dim1] + ti2;
	    cr3 = taui * (cc[i__ - 1 + (k * 3 + 3) * cc_dim1] - cc[ic - 1 + (
		    k * 3 + 2) * cc_dim1]);
	    ci3 = taui * (cc[i__ + (k * 3 + 3) * cc_dim1] + cc[ic + (k * 3 + 
		    2) * cc_dim1]);
	    dr2 = cr2 - ci3;
	    dr3 = cr2 + ci3;
	    di2 = ci2 + cr3;
	    di3 = ci2 - cr3;
	    ch[i__ - 1 + (k + (ch_dim2 << 1)) * ch_dim1] = wa1[i__ - 2] * dr2 
		    - wa1[i__ - 1] * di2;
	    ch[i__ + (k + (ch_dim2 << 1)) * ch_dim1] = wa1[i__ - 2] * di2 + 
		    wa1[i__ - 1] * dr2;
	    ch[i__ - 1 + (k + ch_dim2 * 3) * ch_dim1] = wa2[i__ - 2] * dr3 - 
		    wa2[i__ - 1] * di3;
	    ch[i__ + (k + ch_dim2 * 3) * ch_dim1] = wa2[i__ - 2] * di3 + wa2[
		    i__ - 1] * dr3;
/* L102: */
	}
/* L103: */
    }
    return 0;
} /* radb3_ */

/* Subroutine */ int radb4_(integer *ido, integer *l1, real *cc, real *ch, 
	real *wa1, real *wa2, real *wa3)
{
    /* Initialized data */

    static real sqrt2 = 1.414213562373095f;

    /* System generated locals */
    integer cc_dim1, cc_offset, ch_dim1, ch_dim2, ch_offset, i__1, i__2;

    /* Local variables */
    static integer i__, k, ic;
    static real ci2, ci3, ci4, cr2, cr3, cr4, ti1, ti2, ti3, ti4, tr1, tr2, 
	    tr3, tr4;
    static integer idp2;

    /* Parameter adjustments */
    ch_dim1 = *ido;
    ch_dim2 = *l1;
    ch_offset = 1 + ch_dim1 * (1 + ch_dim2 * 1);
    ch -= ch_offset;
    cc_dim1 = *ido;
    cc_offset = 1 + cc_dim1 * 5;
    cc -= cc_offset;
    --wa1;
    --wa2;
    --wa3;

    /* Function Body */
    i__1 = *l1;
    for (k = 1; k <= i__1; ++k) {
	tr1 = cc[((k << 2) + 1) * cc_dim1 + 1] - cc[*ido + ((k << 2) + 4) * 
		cc_dim1];
	tr2 = cc[((k << 2) + 1) * cc_dim1 + 1] + cc[*ido + ((k << 2) + 4) * 
		cc_dim1];
	tr3 = cc[*ido + ((k << 2) + 2) * cc_dim1] + cc[*ido + ((k << 2) + 2) *
		 cc_dim1];
	tr4 = cc[((k << 2) + 3) * cc_dim1 + 1] + cc[((k << 2) + 3) * cc_dim1 
		+ 1];
	ch[(k + ch_dim2) * ch_dim1 + 1] = tr2 + tr3;
	ch[(k + (ch_dim2 << 1)) * ch_dim1 + 1] = tr1 - tr4;
	ch[(k + ch_dim2 * 3) * ch_dim1 + 1] = tr2 - tr3;
	ch[(k + (ch_dim2 << 2)) * ch_dim1 + 1] = tr1 + tr4;
/* L101: */
    }
    if ((i__1 = *ido - 2) < 0) {
	goto L107;
    } else if (i__1 == 0) {
	goto L105;
    } else {
	goto L102;
    }
L102:
    idp2 = *ido + 2;
    i__1 = *l1;
    for (k = 1; k <= i__1; ++k) {
	i__2 = *ido;
	for (i__ = 3; i__ <= i__2; i__ += 2) {
	    ic = idp2 - i__;
	    ti1 = cc[i__ + ((k << 2) + 1) * cc_dim1] + cc[ic + ((k << 2) + 4) 
		    * cc_dim1];
	    ti2 = cc[i__ + ((k << 2) + 1) * cc_dim1] - cc[ic + ((k << 2) + 4) 
		    * cc_dim1];
	    ti3 = cc[i__ + ((k << 2) + 3) * cc_dim1] - cc[ic + ((k << 2) + 2) 
		    * cc_dim1];
	    tr4 = cc[i__ + ((k << 2) + 3) * cc_dim1] + cc[ic + ((k << 2) + 2) 
		    * cc_dim1];
	    tr1 = cc[i__ - 1 + ((k << 2) + 1) * cc_dim1] - cc[ic - 1 + ((k << 
		    2) + 4) * cc_dim1];
	    tr2 = cc[i__ - 1 + ((k << 2) + 1) * cc_dim1] + cc[ic - 1 + ((k << 
		    2) + 4) * cc_dim1];
	    ti4 = cc[i__ - 1 + ((k << 2) + 3) * cc_dim1] - cc[ic - 1 + ((k << 
		    2) + 2) * cc_dim1];
	    tr3 = cc[i__ - 1 + ((k << 2) + 3) * cc_dim1] + cc[ic - 1 + ((k << 
		    2) + 2) * cc_dim1];
	    ch[i__ - 1 + (k + ch_dim2) * ch_dim1] = tr2 + tr3;
	    cr3 = tr2 - tr3;
	    ch[i__ + (k + ch_dim2) * ch_dim1] = ti2 + ti3;
	    ci3 = ti2 - ti3;
	    cr2 = tr1 - tr4;
	    cr4 = tr1 + tr4;
	    ci2 = ti1 + ti4;
	    ci4 = ti1 - ti4;
	    ch[i__ - 1 + (k + (ch_dim2 << 1)) * ch_dim1] = wa1[i__ - 2] * cr2 
		    - wa1[i__ - 1] * ci2;
	    ch[i__ + (k + (ch_dim2 << 1)) * ch_dim1] = wa1[i__ - 2] * ci2 + 
		    wa1[i__ - 1] * cr2;
	    ch[i__ - 1 + (k + ch_dim2 * 3) * ch_dim1] = wa2[i__ - 2] * cr3 - 
		    wa2[i__ - 1] * ci3;
	    ch[i__ + (k + ch_dim2 * 3) * ch_dim1] = wa2[i__ - 2] * ci3 + wa2[
		    i__ - 1] * cr3;
	    ch[i__ - 1 + (k + (ch_dim2 << 2)) * ch_dim1] = wa3[i__ - 2] * cr4 
		    - wa3[i__ - 1] * ci4;
	    ch[i__ + (k + (ch_dim2 << 2)) * ch_dim1] = wa3[i__ - 2] * ci4 + 
		    wa3[i__ - 1] * cr4;
/* L103: */
	}
/* L104: */
    }
    if (*ido % 2 == 1) {
	return 0;
    }
L105:
    i__1 = *l1;
    for (k = 1; k <= i__1; ++k) {
	ti1 = cc[((k << 2) + 2) * cc_dim1 + 1] + cc[((k << 2) + 4) * cc_dim1 
		+ 1];
	ti2 = cc[((k << 2) + 4) * cc_dim1 + 1] - cc[((k << 2) + 2) * cc_dim1 
		+ 1];
	tr1 = cc[*ido + ((k << 2) + 1) * cc_dim1] - cc[*ido + ((k << 2) + 3) *
		 cc_dim1];
	tr2 = cc[*ido + ((k << 2) + 1) * cc_dim1] + cc[*ido + ((k << 2) + 3) *
		 cc_dim1];
	ch[*ido + (k + ch_dim2) * ch_dim1] = tr2 + tr2;
	ch[*ido + (k + (ch_dim2 << 1)) * ch_dim1] = sqrt2 * (tr1 - ti1);
	ch[*ido + (k + ch_dim2 * 3) * ch_dim1] = ti2 + ti2;
	ch[*ido + (k + (ch_dim2 << 2)) * ch_dim1] = -sqrt2 * (tr1 + ti1);
/* L106: */
    }
L107:
    return 0;
} /* radb4_ */

/* Subroutine */ int radb5_(integer *ido, integer *l1, real *cc, real *ch, 
	real *wa1, real *wa2, real *wa3, real *wa4)
{
    /* Initialized data */

    static real tr11 = .309016994374947f;
    static real ti11 = .951056516295154f;
    static real tr12 = -.809016994374947f;
    static real ti12 = .587785252292473f;

    /* System generated locals */
    integer cc_dim1, cc_offset, ch_dim1, ch_dim2, ch_offset, i__1, i__2;

    /* Local variables */
    static integer i__, k, ic;
    static real ci2, ci3, ci4, ci5, di3, di4, di5, di2, cr2, cr3, cr5, cr4, 
	    ti2, ti3, ti4, ti5, dr3, dr4, dr5, dr2, tr2, tr3, tr4, tr5;
    static integer idp2;

    /* Parameter adjustments */
    ch_dim1 = *ido;
    ch_dim2 = *l1;
    ch_offset = 1 + ch_dim1 * (1 + ch_dim2 * 1);
    ch -= ch_offset;
    cc_dim1 = *ido;
    cc_offset = 1 + cc_dim1 * 6;
    cc -= cc_offset;
    --wa1;
    --wa2;
    --wa3;
    --wa4;

    /* Function Body */
    i__1 = *l1;
    for (k = 1; k <= i__1; ++k) {
	ti5 = cc[(k * 5 + 3) * cc_dim1 + 1] + cc[(k * 5 + 3) * cc_dim1 + 1];
	ti4 = cc[(k * 5 + 5) * cc_dim1 + 1] + cc[(k * 5 + 5) * cc_dim1 + 1];
	tr2 = cc[*ido + (k * 5 + 2) * cc_dim1] + cc[*ido + (k * 5 + 2) * 
		cc_dim1];
	tr3 = cc[*ido + (k * 5 + 4) * cc_dim1] + cc[*ido + (k * 5 + 4) * 
		cc_dim1];
	ch[(k + ch_dim2) * ch_dim1 + 1] = cc[(k * 5 + 1) * cc_dim1 + 1] + tr2 
		+ tr3;
	cr2 = cc[(k * 5 + 1) * cc_dim1 + 1] + tr11 * tr2 + tr12 * tr3;
	cr3 = cc[(k * 5 + 1) * cc_dim1 + 1] + tr12 * tr2 + tr11 * tr3;
	ci5 = ti11 * ti5 + ti12 * ti4;
	ci4 = ti12 * ti5 - ti11 * ti4;
	ch[(k + (ch_dim2 << 1)) * ch_dim1 + 1] = cr2 - ci5;
	ch[(k + ch_dim2 * 3) * ch_dim1 + 1] = cr3 - ci4;
	ch[(k + (ch_dim2 << 2)) * ch_dim1 + 1] = cr3 + ci4;
	ch[(k + ch_dim2 * 5) * ch_dim1 + 1] = cr2 + ci5;
/* L101: */
    }
    if (*ido == 1) {
	return 0;
    }
    idp2 = *ido + 2;
    i__1 = *l1;
    for (k = 1; k <= i__1; ++k) {
	i__2 = *ido;
	for (i__ = 3; i__ <= i__2; i__ += 2) {
	    ic = idp2 - i__;
	    ti5 = cc[i__ + (k * 5 + 3) * cc_dim1] + cc[ic + (k * 5 + 2) * 
		    cc_dim1];
	    ti2 = cc[i__ + (k * 5 + 3) * cc_dim1] - cc[ic + (k * 5 + 2) * 
		    cc_dim1];
	    ti4 = cc[i__ + (k * 5 + 5) * cc_dim1] + cc[ic + (k * 5 + 4) * 
		    cc_dim1];
	    ti3 = cc[i__ + (k * 5 + 5) * cc_dim1] - cc[ic + (k * 5 + 4) * 
		    cc_dim1];
	    tr5 = cc[i__ - 1 + (k * 5 + 3) * cc_dim1] - cc[ic - 1 + (k * 5 + 
		    2) * cc_dim1];
	    tr2 = cc[i__ - 1 + (k * 5 + 3) * cc_dim1] + cc[ic - 1 + (k * 5 + 
		    2) * cc_dim1];
	    tr4 = cc[i__ - 1 + (k * 5 + 5) * cc_dim1] - cc[ic - 1 + (k * 5 + 
		    4) * cc_dim1];
	    tr3 = cc[i__ - 1 + (k * 5 + 5) * cc_dim1] + cc[ic - 1 + (k * 5 + 
		    4) * cc_dim1];
	    ch[i__ - 1 + (k + ch_dim2) * ch_dim1] = cc[i__ - 1 + (k * 5 + 1) *
		     cc_dim1] + tr2 + tr3;
	    ch[i__ + (k + ch_dim2) * ch_dim1] = cc[i__ + (k * 5 + 1) * 
		    cc_dim1] + ti2 + ti3;
	    cr2 = cc[i__ - 1 + (k * 5 + 1) * cc_dim1] + tr11 * tr2 + tr12 * 
		    tr3;
	    ci2 = cc[i__ + (k * 5 + 1) * cc_dim1] + tr11 * ti2 + tr12 * ti3;
	    cr3 = cc[i__ - 1 + (k * 5 + 1) * cc_dim1] + tr12 * tr2 + tr11 * 
		    tr3;
	    ci3 = cc[i__ + (k * 5 + 1) * cc_dim1] + tr12 * ti2 + tr11 * ti3;
	    cr5 = ti11 * tr5 + ti12 * tr4;
	    ci5 = ti11 * ti5 + ti12 * ti4;
	    cr4 = ti12 * tr5 - ti11 * tr4;
	    ci4 = ti12 * ti5 - ti11 * ti4;
	    dr3 = cr3 - ci4;
	    dr4 = cr3 + ci4;
	    di3 = ci3 + cr4;
	    di4 = ci3 - cr4;
	    dr5 = cr2 + ci5;
	    dr2 = cr2 - ci5;
	    di5 = ci2 - cr5;
	    di2 = ci2 + cr5;
	    ch[i__ - 1 + (k + (ch_dim2 << 1)) * ch_dim1] = wa1[i__ - 2] * dr2 
		    - wa1[i__ - 1] * di2;
	    ch[i__ + (k + (ch_dim2 << 1)) * ch_dim1] = wa1[i__ - 2] * di2 + 
		    wa1[i__ - 1] * dr2;
	    ch[i__ - 1 + (k + ch_dim2 * 3) * ch_dim1] = wa2[i__ - 2] * dr3 - 
		    wa2[i__ - 1] * di3;
	    ch[i__ + (k + ch_dim2 * 3) * ch_dim1] = wa2[i__ - 2] * di3 + wa2[
		    i__ - 1] * dr3;
	    ch[i__ - 1 + (k + (ch_dim2 << 2)) * ch_dim1] = wa3[i__ - 2] * dr4 
		    - wa3[i__ - 1] * di4;
	    ch[i__ + (k + (ch_dim2 << 2)) * ch_dim1] = wa3[i__ - 2] * di4 + 
		    wa3[i__ - 1] * dr4;
	    ch[i__ - 1 + (k + ch_dim2 * 5) * ch_dim1] = wa4[i__ - 2] * dr5 - 
		    wa4[i__ - 1] * di5;
	    ch[i__ + (k + ch_dim2 * 5) * ch_dim1] = wa4[i__ - 2] * di5 + wa4[
		    i__ - 1] * dr5;
/* L102: */
	}
/* L103: */
    }
    return 0;
} /* radb5_ */

/* Subroutine */ int radbg_(integer *ido, integer *ip, integer *l1, integer *
	idl1, real *cc, real *c1, real *c2, real *ch, real *ch2, real *wa)
{
    /* Initialized data */

    static real tpi = 6.28318530717959f;

    /* System generated locals */
    integer ch_dim1, ch_dim2, ch_offset, cc_dim1, cc_dim2, cc_offset, c1_dim1,
	     c1_dim2, c1_offset, c2_dim1, c2_offset, ch2_dim1, ch2_offset, 
	    i__1, i__2, i__3;

    /* Builtin functions */
    double cos(doublereal), sin(doublereal);

    /* Local variables */
    static integer idij, ipph, i__, j, k, l, j2, ic, jc, lc, ik, is;
    static real dc2, ai1, ai2, ar1, ar2, ds2;
    static integer nbd;
    static real dcp, arg, dsp, ar1h, ar2h;
    static integer idp2, ipp2;

    /* Parameter adjustments */
    ch_dim1 = *ido;
    ch_dim2 = *l1;
    ch_offset = 1 + ch_dim1 * (1 + ch_dim2 * 1);
    ch -= ch_offset;
    c1_dim1 = *ido;
    c1_dim2 = *l1;
    c1_offset = 1 + c1_dim1 * (1 + c1_dim2 * 1);
    c1 -= c1_offset;
    cc_dim1 = *ido;
    cc_dim2 = *ip;
    cc_offset = 1 + cc_dim1 * (1 + cc_dim2 * 1);
    cc -= cc_offset;
    ch2_dim1 = *idl1;
    ch2_offset = 1 + ch2_dim1 * 1;
    ch2 -= ch2_offset;
    c2_dim1 = *idl1;
    c2_offset = 1 + c2_dim1 * 1;
    c2 -= c2_offset;
    --wa;

    /* Function Body */
    arg = tpi / (real) (*ip);
    dcp = cos(arg);
    dsp = sin(arg);
    idp2 = *ido + 2;
    nbd = (*ido - 1) / 2;
    ipp2 = *ip + 2;
    ipph = (*ip + 1) / 2;
    if (*ido < *l1) {
	goto L103;
    }
    i__1 = *l1;
    for (k = 1; k <= i__1; ++k) {
	i__2 = *ido;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    ch[i__ + (k + ch_dim2) * ch_dim1] = cc[i__ + (k * cc_dim2 + 1) * 
		    cc_dim1];
/* L101: */
	}
/* L102: */
    }
    goto L106;
L103:
    i__1 = *ido;
    for (i__ = 1; i__ <= i__1; ++i__) {
	i__2 = *l1;
	for (k = 1; k <= i__2; ++k) {
	    ch[i__ + (k + ch_dim2) * ch_dim1] = cc[i__ + (k * cc_dim2 + 1) * 
		    cc_dim1];
/* L104: */
	}
/* L105: */
    }
L106:
    i__1 = ipph;
    for (j = 2; j <= i__1; ++j) {
	jc = ipp2 - j;
	j2 = j + j;
	i__2 = *l1;
	for (k = 1; k <= i__2; ++k) {
	    ch[(k + j * ch_dim2) * ch_dim1 + 1] = cc[*ido + (j2 - 2 + k * 
		    cc_dim2) * cc_dim1] + cc[*ido + (j2 - 2 + k * cc_dim2) * 
		    cc_dim1];
	    ch[(k + jc * ch_dim2) * ch_dim1 + 1] = cc[(j2 - 1 + k * cc_dim2) *
		     cc_dim1 + 1] + cc[(j2 - 1 + k * cc_dim2) * cc_dim1 + 1];
/* L107: */
	}
/* L108: */
    }
    if (*ido == 1) {
	goto L116;
    }
    if (nbd < *l1) {
	goto L112;
    }
    i__1 = ipph;
    for (j = 2; j <= i__1; ++j) {
	jc = ipp2 - j;
	i__2 = *l1;
	for (k = 1; k <= i__2; ++k) {
	    i__3 = *ido;
	    for (i__ = 3; i__ <= i__3; i__ += 2) {
		ic = idp2 - i__;
		ch[i__ - 1 + (k + j * ch_dim2) * ch_dim1] = cc[i__ - 1 + ((j 
			<< 1) - 1 + k * cc_dim2) * cc_dim1] + cc[ic - 1 + ((j 
			<< 1) - 2 + k * cc_dim2) * cc_dim1];
		ch[i__ - 1 + (k + jc * ch_dim2) * ch_dim1] = cc[i__ - 1 + ((j 
			<< 1) - 1 + k * cc_dim2) * cc_dim1] - cc[ic - 1 + ((j 
			<< 1) - 2 + k * cc_dim2) * cc_dim1];
		ch[i__ + (k + j * ch_dim2) * ch_dim1] = cc[i__ + ((j << 1) - 
			1 + k * cc_dim2) * cc_dim1] - cc[ic + ((j << 1) - 2 + 
			k * cc_dim2) * cc_dim1];
		ch[i__ + (k + jc * ch_dim2) * ch_dim1] = cc[i__ + ((j << 1) - 
			1 + k * cc_dim2) * cc_dim1] + cc[ic + ((j << 1) - 2 + 
			k * cc_dim2) * cc_dim1];
/* L109: */
	    }
/* L110: */
	}
/* L111: */
    }
    goto L116;
L112:
    i__1 = ipph;
    for (j = 2; j <= i__1; ++j) {
	jc = ipp2 - j;
	i__2 = *ido;
	for (i__ = 3; i__ <= i__2; i__ += 2) {
	    ic = idp2 - i__;
	    i__3 = *l1;
	    for (k = 1; k <= i__3; ++k) {
		ch[i__ - 1 + (k + j * ch_dim2) * ch_dim1] = cc[i__ - 1 + ((j 
			<< 1) - 1 + k * cc_dim2) * cc_dim1] + cc[ic - 1 + ((j 
			<< 1) - 2 + k * cc_dim2) * cc_dim1];
		ch[i__ - 1 + (k + jc * ch_dim2) * ch_dim1] = cc[i__ - 1 + ((j 
			<< 1) - 1 + k * cc_dim2) * cc_dim1] - cc[ic - 1 + ((j 
			<< 1) - 2 + k * cc_dim2) * cc_dim1];
		ch[i__ + (k + j * ch_dim2) * ch_dim1] = cc[i__ + ((j << 1) - 
			1 + k * cc_dim2) * cc_dim1] - cc[ic + ((j << 1) - 2 + 
			k * cc_dim2) * cc_dim1];
		ch[i__ + (k + jc * ch_dim2) * ch_dim1] = cc[i__ + ((j << 1) - 
			1 + k * cc_dim2) * cc_dim1] + cc[ic + ((j << 1) - 2 + 
			k * cc_dim2) * cc_dim1];
/* L113: */
	    }
/* L114: */
	}
/* L115: */
    }
L116:
    ar1 = 1.f;
    ai1 = 0.f;
    i__1 = ipph;
    for (l = 2; l <= i__1; ++l) {
	lc = ipp2 - l;
	ar1h = dcp * ar1 - dsp * ai1;
	ai1 = dcp * ai1 + dsp * ar1;
	ar1 = ar1h;
	i__2 = *idl1;
	for (ik = 1; ik <= i__2; ++ik) {
	    c2[ik + l * c2_dim1] = ch2[ik + ch2_dim1] + ar1 * ch2[ik + (
		    ch2_dim1 << 1)];
	    c2[ik + lc * c2_dim1] = ai1 * ch2[ik + *ip * ch2_dim1];
/* L117: */
	}
	dc2 = ar1;
	ds2 = ai1;
	ar2 = ar1;
	ai2 = ai1;
	i__2 = ipph;
	for (j = 3; j <= i__2; ++j) {
	    jc = ipp2 - j;
	    ar2h = dc2 * ar2 - ds2 * ai2;
	    ai2 = dc2 * ai2 + ds2 * ar2;
	    ar2 = ar2h;
	    i__3 = *idl1;
	    for (ik = 1; ik <= i__3; ++ik) {
		c2[ik + l * c2_dim1] += ar2 * ch2[ik + j * ch2_dim1];
		c2[ik + lc * c2_dim1] += ai2 * ch2[ik + jc * ch2_dim1];
/* L118: */
	    }
/* L119: */
	}
/* L120: */
    }
    i__1 = ipph;
    for (j = 2; j <= i__1; ++j) {
	i__2 = *idl1;
	for (ik = 1; ik <= i__2; ++ik) {
	    ch2[ik + ch2_dim1] += ch2[ik + j * ch2_dim1];
/* L121: */
	}
/* L122: */
    }
    i__1 = ipph;
    for (j = 2; j <= i__1; ++j) {
	jc = ipp2 - j;
	i__2 = *l1;
	for (k = 1; k <= i__2; ++k) {
	    ch[(k + j * ch_dim2) * ch_dim1 + 1] = c1[(k + j * c1_dim2) * 
		    c1_dim1 + 1] - c1[(k + jc * c1_dim2) * c1_dim1 + 1];
	    ch[(k + jc * ch_dim2) * ch_dim1 + 1] = c1[(k + j * c1_dim2) * 
		    c1_dim1 + 1] + c1[(k + jc * c1_dim2) * c1_dim1 + 1];
/* L123: */
	}
/* L124: */
    }
    if (*ido == 1) {
	goto L132;
    }
    if (nbd < *l1) {
	goto L128;
    }
    i__1 = ipph;
    for (j = 2; j <= i__1; ++j) {
	jc = ipp2 - j;
	i__2 = *l1;
	for (k = 1; k <= i__2; ++k) {
	    i__3 = *ido;
	    for (i__ = 3; i__ <= i__3; i__ += 2) {
		ch[i__ - 1 + (k + j * ch_dim2) * ch_dim1] = c1[i__ - 1 + (k + 
			j * c1_dim2) * c1_dim1] - c1[i__ + (k + jc * c1_dim2) 
			* c1_dim1];
		ch[i__ - 1 + (k + jc * ch_dim2) * ch_dim1] = c1[i__ - 1 + (k 
			+ j * c1_dim2) * c1_dim1] + c1[i__ + (k + jc * 
			c1_dim2) * c1_dim1];
		ch[i__ + (k + j * ch_dim2) * ch_dim1] = c1[i__ + (k + j * 
			c1_dim2) * c1_dim1] + c1[i__ - 1 + (k + jc * c1_dim2) 
			* c1_dim1];
		ch[i__ + (k + jc * ch_dim2) * ch_dim1] = c1[i__ + (k + j * 
			c1_dim2) * c1_dim1] - c1[i__ - 1 + (k + jc * c1_dim2) 
			* c1_dim1];
/* L125: */
	    }
/* L126: */
	}
/* L127: */
    }
    goto L132;
L128:
    i__1 = ipph;
    for (j = 2; j <= i__1; ++j) {
	jc = ipp2 - j;
	i__2 = *ido;
	for (i__ = 3; i__ <= i__2; i__ += 2) {
	    i__3 = *l1;
	    for (k = 1; k <= i__3; ++k) {
		ch[i__ - 1 + (k + j * ch_dim2) * ch_dim1] = c1[i__ - 1 + (k + 
			j * c1_dim2) * c1_dim1] - c1[i__ + (k + jc * c1_dim2) 
			* c1_dim1];
		ch[i__ - 1 + (k + jc * ch_dim2) * ch_dim1] = c1[i__ - 1 + (k 
			+ j * c1_dim2) * c1_dim1] + c1[i__ + (k + jc * 
			c1_dim2) * c1_dim1];
		ch[i__ + (k + j * ch_dim2) * ch_dim1] = c1[i__ + (k + j * 
			c1_dim2) * c1_dim1] + c1[i__ - 1 + (k + jc * c1_dim2) 
			* c1_dim1];
		ch[i__ + (k + jc * ch_dim2) * ch_dim1] = c1[i__ + (k + j * 
			c1_dim2) * c1_dim1] - c1[i__ - 1 + (k + jc * c1_dim2) 
			* c1_dim1];
/* L129: */
	    }
/* L130: */
	}
/* L131: */
    }
L132:
    if (*ido == 1) {
	return 0;
    }
    i__1 = *idl1;
    for (ik = 1; ik <= i__1; ++ik) {
	c2[ik + c2_dim1] = ch2[ik + ch2_dim1];
/* L133: */
    }
    i__1 = *ip;
    for (j = 2; j <= i__1; ++j) {
	i__2 = *l1;
	for (k = 1; k <= i__2; ++k) {
	    c1[(k + j * c1_dim2) * c1_dim1 + 1] = ch[(k + j * ch_dim2) * 
		    ch_dim1 + 1];
/* L134: */
	}
/* L135: */
    }
    if (nbd > *l1) {
	goto L139;
    }
    is = -(*ido);
    i__1 = *ip;
    for (j = 2; j <= i__1; ++j) {
	is += *ido;
	idij = is;
	i__2 = *ido;
	for (i__ = 3; i__ <= i__2; i__ += 2) {
	    idij += 2;
	    i__3 = *l1;
	    for (k = 1; k <= i__3; ++k) {
		c1[i__ - 1 + (k + j * c1_dim2) * c1_dim1] = wa[idij - 1] * ch[
			i__ - 1 + (k + j * ch_dim2) * ch_dim1] - wa[idij] * 
			ch[i__ + (k + j * ch_dim2) * ch_dim1];
		c1[i__ + (k + j * c1_dim2) * c1_dim1] = wa[idij - 1] * ch[i__ 
			+ (k + j * ch_dim2) * ch_dim1] + wa[idij] * ch[i__ - 
			1 + (k + j * ch_dim2) * ch_dim1];
/* L136: */
	    }
/* L137: */
	}
/* L138: */
    }
    goto L143;
L139:
    is = -(*ido);
    i__1 = *ip;
    for (j = 2; j <= i__1; ++j) {
	is += *ido;
	i__2 = *l1;
	for (k = 1; k <= i__2; ++k) {
	    idij = is;
	    i__3 = *ido;
	    for (i__ = 3; i__ <= i__3; i__ += 2) {
		idij += 2;
		c1[i__ - 1 + (k + j * c1_dim2) * c1_dim1] = wa[idij - 1] * ch[
			i__ - 1 + (k + j * ch_dim2) * ch_dim1] - wa[idij] * 
			ch[i__ + (k + j * ch_dim2) * ch_dim1];
		c1[i__ + (k + j * c1_dim2) * c1_dim1] = wa[idij - 1] * ch[i__ 
			+ (k + j * ch_dim2) * ch_dim1] + wa[idij] * ch[i__ - 
			1 + (k + j * ch_dim2) * ch_dim1];
/* L140: */
	    }
/* L141: */
	}
/* L142: */
    }
L143:
    return 0;
} /* radbg_ */

/* Subroutine */ int radf2_(integer *ido, integer *l1, real *cc, real *ch, 
	real *wa1)
{
    /* System generated locals */
    integer ch_dim1, ch_offset, cc_dim1, cc_dim2, cc_offset, i__1, i__2;

    /* Local variables */
    static integer i__, k, ic;
    static real ti2, tr2;
    static integer idp2;

    /* Parameter adjustments */
    ch_dim1 = *ido;
    ch_offset = 1 + ch_dim1 * 3;
    ch -= ch_offset;
    cc_dim1 = *ido;
    cc_dim2 = *l1;
    cc_offset = 1 + cc_dim1 * (1 + cc_dim2 * 1);
    cc -= cc_offset;
    --wa1;

    /* Function Body */
    i__1 = *l1;
    for (k = 1; k <= i__1; ++k) {
	ch[((k << 1) + 1) * ch_dim1 + 1] = cc[(k + cc_dim2) * cc_dim1 + 1] + 
		cc[(k + (cc_dim2 << 1)) * cc_dim1 + 1];
	ch[*ido + ((k << 1) + 2) * ch_dim1] = cc[(k + cc_dim2) * cc_dim1 + 1] 
		- cc[(k + (cc_dim2 << 1)) * cc_dim1 + 1];
/* L101: */
    }
    if ((i__1 = *ido - 2) < 0) {
	goto L107;
    } else if (i__1 == 0) {
	goto L105;
    } else {
	goto L102;
    }
L102:
    idp2 = *ido + 2;
    i__1 = *l1;
    for (k = 1; k <= i__1; ++k) {
	i__2 = *ido;
	for (i__ = 3; i__ <= i__2; i__ += 2) {
	    ic = idp2 - i__;
	    tr2 = wa1[i__ - 2] * cc[i__ - 1 + (k + (cc_dim2 << 1)) * cc_dim1] 
		    + wa1[i__ - 1] * cc[i__ + (k + (cc_dim2 << 1)) * cc_dim1];
	    ti2 = wa1[i__ - 2] * cc[i__ + (k + (cc_dim2 << 1)) * cc_dim1] - 
		    wa1[i__ - 1] * cc[i__ - 1 + (k + (cc_dim2 << 1)) * 
		    cc_dim1];
	    ch[i__ + ((k << 1) + 1) * ch_dim1] = cc[i__ + (k + cc_dim2) * 
		    cc_dim1] + ti2;
	    ch[ic + ((k << 1) + 2) * ch_dim1] = ti2 - cc[i__ + (k + cc_dim2) *
		     cc_dim1];
	    ch[i__ - 1 + ((k << 1) + 1) * ch_dim1] = cc[i__ - 1 + (k + 
		    cc_dim2) * cc_dim1] + tr2;
	    ch[ic - 1 + ((k << 1) + 2) * ch_dim1] = cc[i__ - 1 + (k + cc_dim2)
		     * cc_dim1] - tr2;
/* L103: */
	}
/* L104: */
    }
    if (*ido % 2 == 1) {
	return 0;
    }
L105:
    i__1 = *l1;
    for (k = 1; k <= i__1; ++k) {
	ch[((k << 1) + 2) * ch_dim1 + 1] = -cc[*ido + (k + (cc_dim2 << 1)) * 
		cc_dim1];
	ch[*ido + ((k << 1) + 1) * ch_dim1] = cc[*ido + (k + cc_dim2) * 
		cc_dim1];
/* L106: */
    }
L107:
    return 0;
} /* radf2_ */

/* Subroutine */ int radf3_(integer *ido, integer *l1, real *cc, real *ch, 
	real *wa1, real *wa2)
{
    /* Initialized data */

    static real taur = -.5f;
    static real taui = .866025403784439f;

    /* System generated locals */
    integer ch_dim1, ch_offset, cc_dim1, cc_dim2, cc_offset, i__1, i__2;

    /* Local variables */
    static integer i__, k, ic;
    static real ci2, di2, di3, cr2, dr2, dr3, ti2, ti3, tr2, tr3;
    static integer idp2;

    /* Parameter adjustments */
    ch_dim1 = *ido;
    ch_offset = 1 + (ch_dim1 << 2);
    ch -= ch_offset;
    cc_dim1 = *ido;
    cc_dim2 = *l1;
    cc_offset = 1 + cc_dim1 * (1 + cc_dim2 * 1);
    cc -= cc_offset;
    --wa1;
    --wa2;

    /* Function Body */
    i__1 = *l1;
    for (k = 1; k <= i__1; ++k) {
	cr2 = cc[(k + (cc_dim2 << 1)) * cc_dim1 + 1] + cc[(k + cc_dim2 * 3) * 
		cc_dim1 + 1];
	ch[(k * 3 + 1) * ch_dim1 + 1] = cc[(k + cc_dim2) * cc_dim1 + 1] + cr2;
	ch[(k * 3 + 3) * ch_dim1 + 1] = taui * (cc[(k + cc_dim2 * 3) * 
		cc_dim1 + 1] - cc[(k + (cc_dim2 << 1)) * cc_dim1 + 1]);
	ch[*ido + (k * 3 + 2) * ch_dim1] = cc[(k + cc_dim2) * cc_dim1 + 1] + 
		taur * cr2;
/* L101: */
    }
    if (*ido == 1) {
	return 0;
    }
    idp2 = *ido + 2;
    i__1 = *l1;
    for (k = 1; k <= i__1; ++k) {
	i__2 = *ido;
	for (i__ = 3; i__ <= i__2; i__ += 2) {
	    ic = idp2 - i__;
	    dr2 = wa1[i__ - 2] * cc[i__ - 1 + (k + (cc_dim2 << 1)) * cc_dim1] 
		    + wa1[i__ - 1] * cc[i__ + (k + (cc_dim2 << 1)) * cc_dim1];
	    di2 = wa1[i__ - 2] * cc[i__ + (k + (cc_dim2 << 1)) * cc_dim1] - 
		    wa1[i__ - 1] * cc[i__ - 1 + (k + (cc_dim2 << 1)) * 
		    cc_dim1];
	    dr3 = wa2[i__ - 2] * cc[i__ - 1 + (k + cc_dim2 * 3) * cc_dim1] + 
		    wa2[i__ - 1] * cc[i__ + (k + cc_dim2 * 3) * cc_dim1];
	    di3 = wa2[i__ - 2] * cc[i__ + (k + cc_dim2 * 3) * cc_dim1] - wa2[
		    i__ - 1] * cc[i__ - 1 + (k + cc_dim2 * 3) * cc_dim1];
	    cr2 = dr2 + dr3;
	    ci2 = di2 + di3;
	    ch[i__ - 1 + (k * 3 + 1) * ch_dim1] = cc[i__ - 1 + (k + cc_dim2) *
		     cc_dim1] + cr2;
	    ch[i__ + (k * 3 + 1) * ch_dim1] = cc[i__ + (k + cc_dim2) * 
		    cc_dim1] + ci2;
	    tr2 = cc[i__ - 1 + (k + cc_dim2) * cc_dim1] + taur * cr2;
	    ti2 = cc[i__ + (k + cc_dim2) * cc_dim1] + taur * ci2;
	    tr3 = taui * (di2 - di3);
	    ti3 = taui * (dr3 - dr2);
	    ch[i__ - 1 + (k * 3 + 3) * ch_dim1] = tr2 + tr3;
	    ch[ic - 1 + (k * 3 + 2) * ch_dim1] = tr2 - tr3;
	    ch[i__ + (k * 3 + 3) * ch_dim1] = ti2 + ti3;
	    ch[ic + (k * 3 + 2) * ch_dim1] = ti3 - ti2;
/* L102: */
	}
/* L103: */
    }
    return 0;
} /* radf3_ */

/* Subroutine */ int radf4_(integer *ido, integer *l1, real *cc, real *ch, 
	real *wa1, real *wa2, real *wa3)
{
    /* Initialized data */

    static real hsqt2 = .7071067811865475f;

    /* System generated locals */
    integer cc_dim1, cc_dim2, cc_offset, ch_dim1, ch_offset, i__1, i__2;

    /* Local variables */
    static integer i__, k, ic;
    static real ci2, ci3, ci4, cr2, cr3, cr4, ti1, ti2, ti3, ti4, tr1, tr2, 
	    tr3, tr4;
    static integer idp2;

    /* Parameter adjustments */
    ch_dim1 = *ido;
    ch_offset = 1 + ch_dim1 * 5;
    ch -= ch_offset;
    cc_dim1 = *ido;
    cc_dim2 = *l1;
    cc_offset = 1 + cc_dim1 * (1 + cc_dim2 * 1);
    cc -= cc_offset;
    --wa1;
    --wa2;
    --wa3;

    /* Function Body */
    i__1 = *l1;
    for (k = 1; k <= i__1; ++k) {
	tr1 = cc[(k + (cc_dim2 << 1)) * cc_dim1 + 1] + cc[(k + (cc_dim2 << 2))
		 * cc_dim1 + 1];
	tr2 = cc[(k + cc_dim2) * cc_dim1 + 1] + cc[(k + cc_dim2 * 3) * 
		cc_dim1 + 1];
	ch[((k << 2) + 1) * ch_dim1 + 1] = tr1 + tr2;
	ch[*ido + ((k << 2) + 4) * ch_dim1] = tr2 - tr1;
	ch[*ido + ((k << 2) + 2) * ch_dim1] = cc[(k + cc_dim2) * cc_dim1 + 1] 
		- cc[(k + cc_dim2 * 3) * cc_dim1 + 1];
	ch[((k << 2) + 3) * ch_dim1 + 1] = cc[(k + (cc_dim2 << 2)) * cc_dim1 
		+ 1] - cc[(k + (cc_dim2 << 1)) * cc_dim1 + 1];
/* L101: */
    }
    if ((i__1 = *ido - 2) < 0) {
	goto L107;
    } else if (i__1 == 0) {
	goto L105;
    } else {
	goto L102;
    }
L102:
    idp2 = *ido + 2;
    i__1 = *l1;
    for (k = 1; k <= i__1; ++k) {
	i__2 = *ido;
	for (i__ = 3; i__ <= i__2; i__ += 2) {
	    ic = idp2 - i__;
	    cr2 = wa1[i__ - 2] * cc[i__ - 1 + (k + (cc_dim2 << 1)) * cc_dim1] 
		    + wa1[i__ - 1] * cc[i__ + (k + (cc_dim2 << 1)) * cc_dim1];
	    ci2 = wa1[i__ - 2] * cc[i__ + (k + (cc_dim2 << 1)) * cc_dim1] - 
		    wa1[i__ - 1] * cc[i__ - 1 + (k + (cc_dim2 << 1)) * 
		    cc_dim1];
	    cr3 = wa2[i__ - 2] * cc[i__ - 1 + (k + cc_dim2 * 3) * cc_dim1] + 
		    wa2[i__ - 1] * cc[i__ + (k + cc_dim2 * 3) * cc_dim1];
	    ci3 = wa2[i__ - 2] * cc[i__ + (k + cc_dim2 * 3) * cc_dim1] - wa2[
		    i__ - 1] * cc[i__ - 1 + (k + cc_dim2 * 3) * cc_dim1];
	    cr4 = wa3[i__ - 2] * cc[i__ - 1 + (k + (cc_dim2 << 2)) * cc_dim1] 
		    + wa3[i__ - 1] * cc[i__ + (k + (cc_dim2 << 2)) * cc_dim1];
	    ci4 = wa3[i__ - 2] * cc[i__ + (k + (cc_dim2 << 2)) * cc_dim1] - 
		    wa3[i__ - 1] * cc[i__ - 1 + (k + (cc_dim2 << 2)) * 
		    cc_dim1];
	    tr1 = cr2 + cr4;
	    tr4 = cr4 - cr2;
	    ti1 = ci2 + ci4;
	    ti4 = ci2 - ci4;
	    ti2 = cc[i__ + (k + cc_dim2) * cc_dim1] + ci3;
	    ti3 = cc[i__ + (k + cc_dim2) * cc_dim1] - ci3;
	    tr2 = cc[i__ - 1 + (k + cc_dim2) * cc_dim1] + cr3;
	    tr3 = cc[i__ - 1 + (k + cc_dim2) * cc_dim1] - cr3;
	    ch[i__ - 1 + ((k << 2) + 1) * ch_dim1] = tr1 + tr2;
	    ch[ic - 1 + ((k << 2) + 4) * ch_dim1] = tr2 - tr1;
	    ch[i__ + ((k << 2) + 1) * ch_dim1] = ti1 + ti2;
	    ch[ic + ((k << 2) + 4) * ch_dim1] = ti1 - ti2;
	    ch[i__ - 1 + ((k << 2) + 3) * ch_dim1] = ti4 + tr3;
	    ch[ic - 1 + ((k << 2) + 2) * ch_dim1] = tr3 - ti4;
	    ch[i__ + ((k << 2) + 3) * ch_dim1] = tr4 + ti3;
	    ch[ic + ((k << 2) + 2) * ch_dim1] = tr4 - ti3;
/* L103: */
	}
/* L104: */
    }
    if (*ido % 2 == 1) {
	return 0;
    }
L105:
    i__1 = *l1;
    for (k = 1; k <= i__1; ++k) {
	ti1 = -hsqt2 * (cc[*ido + (k + (cc_dim2 << 1)) * cc_dim1] + cc[*ido + 
		(k + (cc_dim2 << 2)) * cc_dim1]);
	tr1 = hsqt2 * (cc[*ido + (k + (cc_dim2 << 1)) * cc_dim1] - cc[*ido + (
		k + (cc_dim2 << 2)) * cc_dim1]);
	ch[*ido + ((k << 2) + 1) * ch_dim1] = tr1 + cc[*ido + (k + cc_dim2) * 
		cc_dim1];
	ch[*ido + ((k << 2) + 3) * ch_dim1] = cc[*ido + (k + cc_dim2) * 
		cc_dim1] - tr1;
	ch[((k << 2) + 2) * ch_dim1 + 1] = ti1 - cc[*ido + (k + cc_dim2 * 3) *
		 cc_dim1];
	ch[((k << 2) + 4) * ch_dim1 + 1] = ti1 + cc[*ido + (k + cc_dim2 * 3) *
		 cc_dim1];
/* L106: */
    }
L107:
    return 0;
} /* radf4_ */

/* Subroutine */ int radf5_(integer *ido, integer *l1, real *cc, real *ch, 
	real *wa1, real *wa2, real *wa3, real *wa4)
{
    /* Initialized data */

    static real tr11 = .309016994374947f;
    static real ti11 = .951056516295154f;
    static real tr12 = -.809016994374947f;
    static real ti12 = .587785252292473f;

    /* System generated locals */
    integer cc_dim1, cc_dim2, cc_offset, ch_dim1, ch_offset, i__1, i__2;

    /* Local variables */
    static integer i__, k, ic;
    static real ci2, di2, ci4, ci5, di3, di4, di5, ci3, cr2, cr3, dr2, dr3, 
	    dr4, dr5, cr5, cr4, ti2, ti3, ti5, ti4, tr2, tr3, tr4, tr5;
    static integer idp2;

    /* Parameter adjustments */
    ch_dim1 = *ido;
    ch_offset = 1 + ch_dim1 * 6;
    ch -= ch_offset;
    cc_dim1 = *ido;
    cc_dim2 = *l1;
    cc_offset = 1 + cc_dim1 * (1 + cc_dim2 * 1);
    cc -= cc_offset;
    --wa1;
    --wa2;
    --wa3;
    --wa4;

    /* Function Body */
    i__1 = *l1;
    for (k = 1; k <= i__1; ++k) {
	cr2 = cc[(k + cc_dim2 * 5) * cc_dim1 + 1] + cc[(k + (cc_dim2 << 1)) * 
		cc_dim1 + 1];
	ci5 = cc[(k + cc_dim2 * 5) * cc_dim1 + 1] - cc[(k + (cc_dim2 << 1)) * 
		cc_dim1 + 1];
	cr3 = cc[(k + (cc_dim2 << 2)) * cc_dim1 + 1] + cc[(k + cc_dim2 * 3) * 
		cc_dim1 + 1];
	ci4 = cc[(k + (cc_dim2 << 2)) * cc_dim1 + 1] - cc[(k + cc_dim2 * 3) * 
		cc_dim1 + 1];
	ch[(k * 5 + 1) * ch_dim1 + 1] = cc[(k + cc_dim2) * cc_dim1 + 1] + cr2 
		+ cr3;
	ch[*ido + (k * 5 + 2) * ch_dim1] = cc[(k + cc_dim2) * cc_dim1 + 1] + 
		tr11 * cr2 + tr12 * cr3;
	ch[(k * 5 + 3) * ch_dim1 + 1] = ti11 * ci5 + ti12 * ci4;
	ch[*ido + (k * 5 + 4) * ch_dim1] = cc[(k + cc_dim2) * cc_dim1 + 1] + 
		tr12 * cr2 + tr11 * cr3;
	ch[(k * 5 + 5) * ch_dim1 + 1] = ti12 * ci5 - ti11 * ci4;
/* L101: */
    }
    if (*ido == 1) {
	return 0;
    }
    idp2 = *ido + 2;
    i__1 = *l1;
    for (k = 1; k <= i__1; ++k) {
	i__2 = *ido;
	for (i__ = 3; i__ <= i__2; i__ += 2) {
	    ic = idp2 - i__;
	    dr2 = wa1[i__ - 2] * cc[i__ - 1 + (k + (cc_dim2 << 1)) * cc_dim1] 
		    + wa1[i__ - 1] * cc[i__ + (k + (cc_dim2 << 1)) * cc_dim1];
	    di2 = wa1[i__ - 2] * cc[i__ + (k + (cc_dim2 << 1)) * cc_dim1] - 
		    wa1[i__ - 1] * cc[i__ - 1 + (k + (cc_dim2 << 1)) * 
		    cc_dim1];
	    dr3 = wa2[i__ - 2] * cc[i__ - 1 + (k + cc_dim2 * 3) * cc_dim1] + 
		    wa2[i__ - 1] * cc[i__ + (k + cc_dim2 * 3) * cc_dim1];
	    di3 = wa2[i__ - 2] * cc[i__ + (k + cc_dim2 * 3) * cc_dim1] - wa2[
		    i__ - 1] * cc[i__ - 1 + (k + cc_dim2 * 3) * cc_dim1];
	    dr4 = wa3[i__ - 2] * cc[i__ - 1 + (k + (cc_dim2 << 2)) * cc_dim1] 
		    + wa3[i__ - 1] * cc[i__ + (k + (cc_dim2 << 2)) * cc_dim1];
	    di4 = wa3[i__ - 2] * cc[i__ + (k + (cc_dim2 << 2)) * cc_dim1] - 
		    wa3[i__ - 1] * cc[i__ - 1 + (k + (cc_dim2 << 2)) * 
		    cc_dim1];
	    dr5 = wa4[i__ - 2] * cc[i__ - 1 + (k + cc_dim2 * 5) * cc_dim1] + 
		    wa4[i__ - 1] * cc[i__ + (k + cc_dim2 * 5) * cc_dim1];
	    di5 = wa4[i__ - 2] * cc[i__ + (k + cc_dim2 * 5) * cc_dim1] - wa4[
		    i__ - 1] * cc[i__ - 1 + (k + cc_dim2 * 5) * cc_dim1];
	    cr2 = dr2 + dr5;
	    ci5 = dr5 - dr2;
	    cr5 = di2 - di5;
	    ci2 = di2 + di5;
	    cr3 = dr3 + dr4;
	    ci4 = dr4 - dr3;
	    cr4 = di3 - di4;
	    ci3 = di3 + di4;
	    ch[i__ - 1 + (k * 5 + 1) * ch_dim1] = cc[i__ - 1 + (k + cc_dim2) *
		     cc_dim1] + cr2 + cr3;
	    ch[i__ + (k * 5 + 1) * ch_dim1] = cc[i__ + (k + cc_dim2) * 
		    cc_dim1] + ci2 + ci3;
	    tr2 = cc[i__ - 1 + (k + cc_dim2) * cc_dim1] + tr11 * cr2 + tr12 * 
		    cr3;
	    ti2 = cc[i__ + (k + cc_dim2) * cc_dim1] + tr11 * ci2 + tr12 * ci3;
	    tr3 = cc[i__ - 1 + (k + cc_dim2) * cc_dim1] + tr12 * cr2 + tr11 * 
		    cr3;
	    ti3 = cc[i__ + (k + cc_dim2) * cc_dim1] + tr12 * ci2 + tr11 * ci3;
	    tr5 = ti11 * cr5 + ti12 * cr4;
	    ti5 = ti11 * ci5 + ti12 * ci4;
	    tr4 = ti12 * cr5 - ti11 * cr4;
	    ti4 = ti12 * ci5 - ti11 * ci4;
	    ch[i__ - 1 + (k * 5 + 3) * ch_dim1] = tr2 + tr5;
	    ch[ic - 1 + (k * 5 + 2) * ch_dim1] = tr2 - tr5;
	    ch[i__ + (k * 5 + 3) * ch_dim1] = ti2 + ti5;
	    ch[ic + (k * 5 + 2) * ch_dim1] = ti5 - ti2;
	    ch[i__ - 1 + (k * 5 + 5) * ch_dim1] = tr3 + tr4;
	    ch[ic - 1 + (k * 5 + 4) * ch_dim1] = tr3 - tr4;
	    ch[i__ + (k * 5 + 5) * ch_dim1] = ti3 + ti4;
	    ch[ic + (k * 5 + 4) * ch_dim1] = ti4 - ti3;
/* L102: */
	}
/* L103: */
    }
    return 0;
} /* radf5_ */

/* Subroutine */ int radfg_(integer *ido, integer *ip, integer *l1, integer *
	idl1, real *cc, real *c1, real *c2, real *ch, real *ch2, real *wa)
{
    /* Initialized data */

    static real tpi = 6.28318530717959f;

    /* System generated locals */
    integer ch_dim1, ch_dim2, ch_offset, cc_dim1, cc_dim2, cc_offset, c1_dim1,
	     c1_dim2, c1_offset, c2_dim1, c2_offset, ch2_dim1, ch2_offset, 
	    i__1, i__2, i__3;

    /* Builtin functions */
    double cos(doublereal), sin(doublereal);

    /* Local variables */
    static integer idij, ipph, i__, j, k, l, j2, ic, jc, lc, ik, is;
    static real dc2, ai1, ai2, ar1, ar2, ds2;
    static integer nbd;
    static real dcp, arg, dsp, ar1h, ar2h;
    static integer idp2, ipp2;

    /* Parameter adjustments */
    ch_dim1 = *ido;
    ch_dim2 = *l1;
    ch_offset = 1 + ch_dim1 * (1 + ch_dim2 * 1);
    ch -= ch_offset;
    c1_dim1 = *ido;
    c1_dim2 = *l1;
    c1_offset = 1 + c1_dim1 * (1 + c1_dim2 * 1);
    c1 -= c1_offset;
    cc_dim1 = *ido;
    cc_dim2 = *ip;
    cc_offset = 1 + cc_dim1 * (1 + cc_dim2 * 1);
    cc -= cc_offset;
    ch2_dim1 = *idl1;
    ch2_offset = 1 + ch2_dim1 * 1;
    ch2 -= ch2_offset;
    c2_dim1 = *idl1;
    c2_offset = 1 + c2_dim1 * 1;
    c2 -= c2_offset;
    --wa;

    /* Function Body */
    arg = tpi / (real) (*ip);
    dcp = cos(arg);
    dsp = sin(arg);
    ipph = (*ip + 1) / 2;
    ipp2 = *ip + 2;
    idp2 = *ido + 2;
    nbd = (*ido - 1) / 2;
    if (*ido == 1) {
	goto L119;
    }
    i__1 = *idl1;
    for (ik = 1; ik <= i__1; ++ik) {
	ch2[ik + ch2_dim1] = c2[ik + c2_dim1];
/* L101: */
    }
    i__1 = *ip;
    for (j = 2; j <= i__1; ++j) {
	i__2 = *l1;
	for (k = 1; k <= i__2; ++k) {
	    ch[(k + j * ch_dim2) * ch_dim1 + 1] = c1[(k + j * c1_dim2) * 
		    c1_dim1 + 1];
/* L102: */
	}
/* L103: */
    }
    if (nbd > *l1) {
	goto L107;
    }
    is = -(*ido);
    i__1 = *ip;
    for (j = 2; j <= i__1; ++j) {
	is += *ido;
	idij = is;
	i__2 = *ido;
	for (i__ = 3; i__ <= i__2; i__ += 2) {
	    idij += 2;
	    i__3 = *l1;
	    for (k = 1; k <= i__3; ++k) {
		ch[i__ - 1 + (k + j * ch_dim2) * ch_dim1] = wa[idij - 1] * c1[
			i__ - 1 + (k + j * c1_dim2) * c1_dim1] + wa[idij] * 
			c1[i__ + (k + j * c1_dim2) * c1_dim1];
		ch[i__ + (k + j * ch_dim2) * ch_dim1] = wa[idij - 1] * c1[i__ 
			+ (k + j * c1_dim2) * c1_dim1] - wa[idij] * c1[i__ - 
			1 + (k + j * c1_dim2) * c1_dim1];
/* L104: */
	    }
/* L105: */
	}
/* L106: */
    }
    goto L111;
L107:
    is = -(*ido);
    i__1 = *ip;
    for (j = 2; j <= i__1; ++j) {
	is += *ido;
	i__2 = *l1;
	for (k = 1; k <= i__2; ++k) {
	    idij = is;
	    i__3 = *ido;
	    for (i__ = 3; i__ <= i__3; i__ += 2) {
		idij += 2;
		ch[i__ - 1 + (k + j * ch_dim2) * ch_dim1] = wa[idij - 1] * c1[
			i__ - 1 + (k + j * c1_dim2) * c1_dim1] + wa[idij] * 
			c1[i__ + (k + j * c1_dim2) * c1_dim1];
		ch[i__ + (k + j * ch_dim2) * ch_dim1] = wa[idij - 1] * c1[i__ 
			+ (k + j * c1_dim2) * c1_dim1] - wa[idij] * c1[i__ - 
			1 + (k + j * c1_dim2) * c1_dim1];
/* L108: */
	    }
/* L109: */
	}
/* L110: */
    }
L111:
    if (nbd < *l1) {
	goto L115;
    }
    i__1 = ipph;
    for (j = 2; j <= i__1; ++j) {
	jc = ipp2 - j;
	i__2 = *l1;
	for (k = 1; k <= i__2; ++k) {
	    i__3 = *ido;
	    for (i__ = 3; i__ <= i__3; i__ += 2) {
		c1[i__ - 1 + (k + j * c1_dim2) * c1_dim1] = ch[i__ - 1 + (k + 
			j * ch_dim2) * ch_dim1] + ch[i__ - 1 + (k + jc * 
			ch_dim2) * ch_dim1];
		c1[i__ - 1 + (k + jc * c1_dim2) * c1_dim1] = ch[i__ + (k + j *
			 ch_dim2) * ch_dim1] - ch[i__ + (k + jc * ch_dim2) * 
			ch_dim1];
		c1[i__ + (k + j * c1_dim2) * c1_dim1] = ch[i__ + (k + j * 
			ch_dim2) * ch_dim1] + ch[i__ + (k + jc * ch_dim2) * 
			ch_dim1];
		c1[i__ + (k + jc * c1_dim2) * c1_dim1] = ch[i__ - 1 + (k + jc 
			* ch_dim2) * ch_dim1] - ch[i__ - 1 + (k + j * ch_dim2)
			 * ch_dim1];
/* L112: */
	    }
/* L113: */
	}
/* L114: */
    }
    goto L121;
L115:
    i__1 = ipph;
    for (j = 2; j <= i__1; ++j) {
	jc = ipp2 - j;
	i__2 = *ido;
	for (i__ = 3; i__ <= i__2; i__ += 2) {
	    i__3 = *l1;
	    for (k = 1; k <= i__3; ++k) {
		c1[i__ - 1 + (k + j * c1_dim2) * c1_dim1] = ch[i__ - 1 + (k + 
			j * ch_dim2) * ch_dim1] + ch[i__ - 1 + (k + jc * 
			ch_dim2) * ch_dim1];
		c1[i__ - 1 + (k + jc * c1_dim2) * c1_dim1] = ch[i__ + (k + j *
			 ch_dim2) * ch_dim1] - ch[i__ + (k + jc * ch_dim2) * 
			ch_dim1];
		c1[i__ + (k + j * c1_dim2) * c1_dim1] = ch[i__ + (k + j * 
			ch_dim2) * ch_dim1] + ch[i__ + (k + jc * ch_dim2) * 
			ch_dim1];
		c1[i__ + (k + jc * c1_dim2) * c1_dim1] = ch[i__ - 1 + (k + jc 
			* ch_dim2) * ch_dim1] - ch[i__ - 1 + (k + j * ch_dim2)
			 * ch_dim1];
/* L116: */
	    }
/* L117: */
	}
/* L118: */
    }
    goto L121;
L119:
    i__1 = *idl1;
    for (ik = 1; ik <= i__1; ++ik) {
	c2[ik + c2_dim1] = ch2[ik + ch2_dim1];
/* L120: */
    }
L121:
    i__1 = ipph;
    for (j = 2; j <= i__1; ++j) {
	jc = ipp2 - j;
	i__2 = *l1;
	for (k = 1; k <= i__2; ++k) {
	    c1[(k + j * c1_dim2) * c1_dim1 + 1] = ch[(k + j * ch_dim2) * 
		    ch_dim1 + 1] + ch[(k + jc * ch_dim2) * ch_dim1 + 1];
	    c1[(k + jc * c1_dim2) * c1_dim1 + 1] = ch[(k + jc * ch_dim2) * 
		    ch_dim1 + 1] - ch[(k + j * ch_dim2) * ch_dim1 + 1];
/* L122: */
	}
/* L123: */
    }

    ar1 = 1.f;
    ai1 = 0.f;
    i__1 = ipph;
    for (l = 2; l <= i__1; ++l) {
	lc = ipp2 - l;
	ar1h = dcp * ar1 - dsp * ai1;
	ai1 = dcp * ai1 + dsp * ar1;
	ar1 = ar1h;
	i__2 = *idl1;
	for (ik = 1; ik <= i__2; ++ik) {
	    ch2[ik + l * ch2_dim1] = c2[ik + c2_dim1] + ar1 * c2[ik + (
		    c2_dim1 << 1)];
	    ch2[ik + lc * ch2_dim1] = ai1 * c2[ik + *ip * c2_dim1];
/* L124: */
	}
	dc2 = ar1;
	ds2 = ai1;
	ar2 = ar1;
	ai2 = ai1;
	i__2 = ipph;
	for (j = 3; j <= i__2; ++j) {
	    jc = ipp2 - j;
	    ar2h = dc2 * ar2 - ds2 * ai2;
	    ai2 = dc2 * ai2 + ds2 * ar2;
	    ar2 = ar2h;
	    i__3 = *idl1;
	    for (ik = 1; ik <= i__3; ++ik) {
		ch2[ik + l * ch2_dim1] += ar2 * c2[ik + j * c2_dim1];
		ch2[ik + lc * ch2_dim1] += ai2 * c2[ik + jc * c2_dim1];
/* L125: */
	    }
/* L126: */
	}
/* L127: */
    }
    i__1 = ipph;
    for (j = 2; j <= i__1; ++j) {
	i__2 = *idl1;
	for (ik = 1; ik <= i__2; ++ik) {
	    ch2[ik + ch2_dim1] += c2[ik + j * c2_dim1];
/* L128: */
	}
/* L129: */
    }

    if (*ido < *l1) {
	goto L132;
    }
    i__1 = *l1;
    for (k = 1; k <= i__1; ++k) {
	i__2 = *ido;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    cc[i__ + (k * cc_dim2 + 1) * cc_dim1] = ch[i__ + (k + ch_dim2) * 
		    ch_dim1];
/* L130: */
	}
/* L131: */
    }
    goto L135;
L132:
    i__1 = *ido;
    for (i__ = 1; i__ <= i__1; ++i__) {
	i__2 = *l1;
	for (k = 1; k <= i__2; ++k) {
	    cc[i__ + (k * cc_dim2 + 1) * cc_dim1] = ch[i__ + (k + ch_dim2) * 
		    ch_dim1];
/* L133: */
	}
/* L134: */
    }
L135:
    i__1 = ipph;
    for (j = 2; j <= i__1; ++j) {
	jc = ipp2 - j;
	j2 = j + j;
	i__2 = *l1;
	for (k = 1; k <= i__2; ++k) {
	    cc[*ido + (j2 - 2 + k * cc_dim2) * cc_dim1] = ch[(k + j * ch_dim2)
		     * ch_dim1 + 1];
	    cc[(j2 - 1 + k * cc_dim2) * cc_dim1 + 1] = ch[(k + jc * ch_dim2) *
		     ch_dim1 + 1];
/* L136: */
	}
/* L137: */
    }
    if (*ido == 1) {
	return 0;
    }
    if (nbd < *l1) {
	goto L141;
    }
    i__1 = ipph;
    for (j = 2; j <= i__1; ++j) {
	jc = ipp2 - j;
	j2 = j + j;
	i__2 = *l1;
	for (k = 1; k <= i__2; ++k) {
	    i__3 = *ido;
	    for (i__ = 3; i__ <= i__3; i__ += 2) {
		ic = idp2 - i__;
		cc[i__ - 1 + (j2 - 1 + k * cc_dim2) * cc_dim1] = ch[i__ - 1 + 
			(k + j * ch_dim2) * ch_dim1] + ch[i__ - 1 + (k + jc * 
			ch_dim2) * ch_dim1];
		cc[ic - 1 + (j2 - 2 + k * cc_dim2) * cc_dim1] = ch[i__ - 1 + (
			k + j * ch_dim2) * ch_dim1] - ch[i__ - 1 + (k + jc * 
			ch_dim2) * ch_dim1];
		cc[i__ + (j2 - 1 + k * cc_dim2) * cc_dim1] = ch[i__ + (k + j *
			 ch_dim2) * ch_dim1] + ch[i__ + (k + jc * ch_dim2) * 
			ch_dim1];
		cc[ic + (j2 - 2 + k * cc_dim2) * cc_dim1] = ch[i__ + (k + jc *
			 ch_dim2) * ch_dim1] - ch[i__ + (k + j * ch_dim2) * 
			ch_dim1];
/* L138: */
	    }
/* L139: */
	}
/* L140: */
    }
    return 0;
L141:
    i__1 = ipph;
    for (j = 2; j <= i__1; ++j) {
	jc = ipp2 - j;
	j2 = j + j;
	i__2 = *ido;
	for (i__ = 3; i__ <= i__2; i__ += 2) {
	    ic = idp2 - i__;
	    i__3 = *l1;
	    for (k = 1; k <= i__3; ++k) {
		cc[i__ - 1 + (j2 - 1 + k * cc_dim2) * cc_dim1] = ch[i__ - 1 + 
			(k + j * ch_dim2) * ch_dim1] + ch[i__ - 1 + (k + jc * 
			ch_dim2) * ch_dim1];
		cc[ic - 1 + (j2 - 2 + k * cc_dim2) * cc_dim1] = ch[i__ - 1 + (
			k + j * ch_dim2) * ch_dim1] - ch[i__ - 1 + (k + jc * 
			ch_dim2) * ch_dim1];
		cc[i__ + (j2 - 1 + k * cc_dim2) * cc_dim1] = ch[i__ + (k + j *
			 ch_dim2) * ch_dim1] + ch[i__ + (k + jc * ch_dim2) * 
			ch_dim1];
		cc[ic + (j2 - 2 + k * cc_dim2) * cc_dim1] = ch[i__ + (k + jc *
			 ch_dim2) * ch_dim1] - ch[i__ + (k + j * ch_dim2) * 
			ch_dim1];
/* L142: */
	    }
/* L143: */
	}
/* L144: */
    }
    return 0;
} /* radfg_ */

/* Subroutine */ int rfftb_(integer *n, real *r__, real *wsave)
{
    extern /* Subroutine */ int rfftb1_(integer *, real *, real *, real *, 
	    integer *);

    /* Parameter adjustments */
    --wsave;
    --r__;

    /* Function Body */
    if (*n == 1) {
	return 0;
    }
    rfftb1_(n, &r__[1], &wsave[1], &wsave[*n + 1], (integer*)&wsave[(*n << 1) 
	    + 1]);
/*     call rfftb1 (n,r,wsave,wsave(n+1),int(wsave(2*n+1))) */
    return 0;
} /* rfftb_ */

/* Subroutine */ int rfftb1_(integer *n, real *c__, real *ch, real *wa, 
	integer *ifac)
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    extern /* Subroutine */ int radb2_(integer *, integer *, real *, real *, 
	    real *), radb3_(integer *, integer *, real *, real *, real *, 
	    real *), radb4_(integer *, integer *, real *, real *, real *, 
	    real *, real *), radb5_(integer *, integer *, real *, real *, 
	    real *, real *, real *, real *);
    static integer i__;
    extern /* Subroutine */ int radbg_(integer *, integer *, integer *, 
	    integer *, real *, real *, real *, real *, real *, real *);
    static integer k1, l1, l2, na, nf, ip, iw, ix2, ix3, ix4, ido, idl1;

    /* Parameter adjustments */
    --ifac;
    --wa;
    --ch;
    --c__;

    /* Function Body */
    nf = ifac[2];
    na = 0;
    l1 = 1;
    iw = 1;
    i__1 = nf;
    for (k1 = 1; k1 <= i__1; ++k1) {
	ip = ifac[k1 + 2];
	l2 = ip * l1;
	ido = *n / l2;
	idl1 = ido * l1;
	if (ip != 4) {
	    goto L103;
	}
	ix2 = iw + ido;
	ix3 = ix2 + ido;
	if (na != 0) {
	    goto L101;
	}
	radb4_(&ido, &l1, &c__[1], &ch[1], &wa[iw], &wa[ix2], &wa[ix3]);
	goto L102;
L101:
	radb4_(&ido, &l1, &ch[1], &c__[1], &wa[iw], &wa[ix2], &wa[ix3]);
L102:
	na = 1 - na;
	goto L115;
L103:
	if (ip != 2) {
	    goto L106;
	}
	if (na != 0) {
	    goto L104;
	}
	radb2_(&ido, &l1, &c__[1], &ch[1], &wa[iw]);
	goto L105;
L104:
	radb2_(&ido, &l1, &ch[1], &c__[1], &wa[iw]);
L105:
	na = 1 - na;
	goto L115;
L106:
	if (ip != 3) {
	    goto L109;
	}
	ix2 = iw + ido;
	if (na != 0) {
	    goto L107;
	}
	radb3_(&ido, &l1, &c__[1], &ch[1], &wa[iw], &wa[ix2]);
	goto L108;
L107:
	radb3_(&ido, &l1, &ch[1], &c__[1], &wa[iw], &wa[ix2]);
L108:
	na = 1 - na;
	goto L115;
L109:
	if (ip != 5) {
	    goto L112;
	}
	ix2 = iw + ido;
	ix3 = ix2 + ido;
	ix4 = ix3 + ido;
	if (na != 0) {
	    goto L110;
	}
	radb5_(&ido, &l1, &c__[1], &ch[1], &wa[iw], &wa[ix2], &wa[ix3], &wa[
		ix4]);
	goto L111;
L110:
	radb5_(&ido, &l1, &ch[1], &c__[1], &wa[iw], &wa[ix2], &wa[ix3], &wa[
		ix4]);
L111:
	na = 1 - na;
	goto L115;
L112:
	if (na != 0) {
	    goto L113;
	}
	radbg_(&ido, &ip, &l1, &idl1, &c__[1], &c__[1], &c__[1], &ch[1], &ch[
		1], &wa[iw]);
	goto L114;
L113:
	radbg_(&ido, &ip, &l1, &idl1, &ch[1], &ch[1], &ch[1], &c__[1], &c__[1]
		, &wa[iw]);
L114:
	if (ido == 1) {
	    na = 1 - na;
	}
L115:
	l1 = l2;
	iw += (ip - 1) * ido;
/* L116: */
    }
    if (na == 0) {
	return 0;
    }
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	c__[i__] = ch[i__];
/* L117: */
    }
    return 0;
} /* rfftb1_ */

/* Subroutine */ int rfftf_(integer *n, real *r__, real *wsave)
{
    extern /* Subroutine */ int rfftf1_(integer *, real *, real *, real *, 
	    integer *);

    /* Parameter adjustments */
    --wsave;
    --r__;

    /* Function Body */
    if (*n == 1) {
	return 0;
    }
    rfftf1_(n, &r__[1], &wsave[1], &wsave[*n + 1], (integer*)&wsave[(*n << 1) 
	    + 1]);
/*     call rfftf1 (n,r,wsave,wsave(n+1),int(wsave(2*n+1))) */
    return 0;
} /* rfftf_ */

/* Subroutine */ int rfftf1_(integer *n, real *c__, real *ch, real *wa, 
	integer *ifac)
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    extern /* Subroutine */ int radf2_(integer *, integer *, real *, real *, 
	    real *), radf3_(integer *, integer *, real *, real *, real *, 
	    real *), radf4_(integer *, integer *, real *, real *, real *, 
	    real *, real *), radf5_(integer *, integer *, real *, real *, 
	    real *, real *, real *, real *);
    static integer i__;
    extern /* Subroutine */ int radfg_(integer *, integer *, integer *, 
	    integer *, real *, real *, real *, real *, real *, real *);
    static integer k1, l1, l2, na, kh, nf, ip, iw, ix2, ix3, ix4, ido, idl1;

    /* Parameter adjustments */
    --ifac;
    --wa;
    --ch;
    --c__;

    /* Function Body */
    nf = ifac[2];
    na = 1;
    l2 = *n;
    iw = *n;
    i__1 = nf;
    for (k1 = 1; k1 <= i__1; ++k1) {
	kh = nf - k1;
	ip = ifac[kh + 3];
	l1 = l2 / ip;
	ido = *n / l2;
	idl1 = ido * l1;
	iw -= (ip - 1) * ido;
	na = 1 - na;
	if (ip != 4) {
	    goto L102;
	}
	ix2 = iw + ido;
	ix3 = ix2 + ido;
	if (na != 0) {
	    goto L101;
	}
	radf4_(&ido, &l1, &c__[1], &ch[1], &wa[iw], &wa[ix2], &wa[ix3]);
	goto L110;
L101:
	radf4_(&ido, &l1, &ch[1], &c__[1], &wa[iw], &wa[ix2], &wa[ix3]);
	goto L110;
L102:
	if (ip != 2) {
	    goto L104;
	}
	if (na != 0) {
	    goto L103;
	}
	radf2_(&ido, &l1, &c__[1], &ch[1], &wa[iw]);
	goto L110;
L103:
	radf2_(&ido, &l1, &ch[1], &c__[1], &wa[iw]);
	goto L110;
L104:
	if (ip != 3) {
	    goto L106;
	}
	ix2 = iw + ido;
	if (na != 0) {
	    goto L105;
	}
	radf3_(&ido, &l1, &c__[1], &ch[1], &wa[iw], &wa[ix2]);
	goto L110;
L105:
	radf3_(&ido, &l1, &ch[1], &c__[1], &wa[iw], &wa[ix2]);
	goto L110;
L106:
	if (ip != 5) {
	    goto L108;
	}
	ix2 = iw + ido;
	ix3 = ix2 + ido;
	ix4 = ix3 + ido;
	if (na != 0) {
	    goto L107;
	}
	radf5_(&ido, &l1, &c__[1], &ch[1], &wa[iw], &wa[ix2], &wa[ix3], &wa[
		ix4]);
	goto L110;
L107:
	radf5_(&ido, &l1, &ch[1], &c__[1], &wa[iw], &wa[ix2], &wa[ix3], &wa[
		ix4]);
	goto L110;
L108:
	if (ido == 1) {
	    na = 1 - na;
	}
	if (na != 0) {
	    goto L109;
	}
	radfg_(&ido, &ip, &l1, &idl1, &c__[1], &c__[1], &c__[1], &ch[1], &ch[
		1], &wa[iw]);
	na = 1;
	goto L110;
L109:
	radfg_(&ido, &ip, &l1, &idl1, &ch[1], &ch[1], &ch[1], &c__[1], &c__[1]
		, &wa[iw]);
	na = 0;
L110:
	l2 = l1;
/* L111: */
    }
    if (na == 1) {
	return 0;
    }
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	c__[i__] = ch[i__];
/* L112: */
    }
    return 0;
} /* rfftf1_ */

/* Subroutine */ int rffti_(integer *n, real *wsave)
{
    extern /* Subroutine */ int rffti1_(integer *, real *, integer *);

    /* Parameter adjustments */
    --wsave;

    /* Function Body */
    if (*n == 1) {
	return 0;
    }
    rffti1_(n, &wsave[*n + 1], (integer*)&wsave[(*n << 1) + 1]);
/*     call rffti1 (n,wsave(n+1),int(wsave(2*n+1))) */
    return 0;
} /* rffti_ */

/* Subroutine */ int rffti1_(integer *n, real *wa, integer *ifac)
{
    /* Initialized data */

    static integer ntryh[4] = { 4,2,3,5 };

    /* System generated locals */
    integer i__1, i__2, i__3;

    /* Builtin functions */
    double cos(doublereal), sin(doublereal);

    /* Local variables */
    static real argh;
    static integer ntry, i__, j;
    static real argld;
    static integer k1, l1, l2, ib;
    static real fi;
    static integer ld, ii, nf, ip, nl, is, nq, nr;
    static real arg;
    static integer ido, ipm;
    static real tpi;
    static integer nfm1;

    /* Parameter adjustments */
    --ifac;
    --wa;

    /* Function Body */
    nl = *n;
    nf = 0;
    j = 0;
L101:
    ++j;
    if (j - 4 <= 0) {
	goto L102;
    } else {
	goto L103;
    }
L102:
    ntry = ntryh[j - 1];
    goto L104;
L103:
    ntry += 2;
L104:
    nq = nl / ntry;
    nr = nl - ntry * nq;
    if (nr != 0) {
	goto L101;
    } else {
	goto L105;
    }
L105:
    ++nf;
    ifac[nf + 2] = ntry;
    nl = nq;
    if (ntry != 2) {
	goto L107;
    }
    if (nf == 1) {
	goto L107;
    }
    i__1 = nf;
    for (i__ = 2; i__ <= i__1; ++i__) {
	ib = nf - i__ + 2;
	ifac[ib + 2] = ifac[ib + 1];
/* L106: */
    }
    ifac[3] = 2;
L107:
    if (nl != 1) {
	goto L104;
    }
    ifac[1] = *n;
    ifac[2] = nf;
    tpi = 6.28318530717959f;
    argh = tpi / (real) (*n);
    is = 0;
    nfm1 = nf - 1;
    l1 = 1;
    if (nfm1 == 0) {
	return 0;
    }
    i__1 = nfm1;
    for (k1 = 1; k1 <= i__1; ++k1) {
	ip = ifac[k1 + 2];
	ld = 0;
	l2 = l1 * ip;
	ido = *n / l2;
	ipm = ip - 1;
	i__2 = ipm;
	for (j = 1; j <= i__2; ++j) {
	    ld += l1;
	    i__ = is;
	    argld = (real) ld * argh;
	    fi = 0.f;
	    i__3 = ido;
	    for (ii = 3; ii <= i__3; ii += 2) {
		i__ += 2;
		fi += 1.f;
		arg = fi * argld;
		wa[i__ - 1] = cos(arg);
		wa[i__] = sin(arg);
/* L108: */
	    }
	    is += ido;
/* L109: */
	}
	l1 = l2;
/* L110: */
    }
    return 0;
} /* rffti1_ */

/* Subroutine */ int sinqb_(integer *n, real *x, real *wsave)
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    static integer k;
    extern /* Subroutine */ int cosqb_(integer *, real *, real *);
    static real xhold;
    static integer kc, ns2;

    /* Parameter adjustments */
    --wsave;
    --x;

    /* Function Body */
    if (*n > 1) {
	goto L101;
    }
    x[1] *= 4.f;
    return 0;
L101:
    ns2 = *n / 2;
    i__1 = *n;
    for (k = 2; k <= i__1; k += 2) {
	x[k] = -x[k];
/* L102: */
    }
    cosqb_(n, &x[1], &wsave[1]);
    i__1 = ns2;
    for (k = 1; k <= i__1; ++k) {
	kc = *n - k;
	xhold = x[k];
	x[k] = x[kc + 1];
	x[kc + 1] = xhold;
/* L103: */
    }
    return 0;
} /* sinqb_ */

/* Subroutine */ int sinqf_(integer *n, real *x, real *wsave)
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    static integer k;
    extern /* Subroutine */ int cosqf_(integer *, real *, real *);
    static real xhold;
    static integer kc, ns2;

    /* Parameter adjustments */
    --wsave;
    --x;

    /* Function Body */
    if (*n == 1) {
	return 0;
    }
    ns2 = *n / 2;
    i__1 = ns2;
    for (k = 1; k <= i__1; ++k) {
	kc = *n - k;
	xhold = x[k];
	x[k] = x[kc + 1];
	x[kc + 1] = xhold;
/* L101: */
    }
    cosqf_(n, &x[1], &wsave[1]);
    i__1 = *n;
    for (k = 2; k <= i__1; k += 2) {
	x[k] = -x[k];
/* L102: */
    }
    return 0;
} /* sinqf_ */

/* Subroutine */ int sinqi_(integer *n, real *wsave)
{
    extern /* Subroutine */ int cosqi_(integer *, real *);

    /* Parameter adjustments */
    --wsave;

    /* Function Body */
    cosqi_(n, &wsave[1]);
    return 0;
} /* sinqi_ */

/* Subroutine */ int sint_(integer *n, real *x, real *wsave)
{
    extern /* Subroutine */ int sint1_(integer *, real *, real *, real *, 
	    real *, integer *);
    static integer np1, iw1, iw2, iw3;

    /* Parameter adjustments */
    --wsave;
    --x;

    /* Function Body */
    np1 = *n + 1;
    iw1 = *n / 2 + 1;
    iw2 = iw1 + np1;
    iw3 = iw2 + np1;
    sint1_(n, &x[1], &wsave[1], &wsave[iw1], &wsave[iw2], (integer*)&wsave[
	    iw3]);
/*     call sint1(n,x,wsave,wsave(iw1),wsave(iw2),int(wsave(iw3))) */
    return 0;
} /* sint_ */

/* Subroutine */ int sint1_(integer *n, real *war, real *was, real *xh, real *
	x, integer *ifac)
{
    /* Initialized data */

    static real sqrt3 = 1.73205080756888f;

    /* System generated locals */
    integer i__1;

    /* Local variables */
    static integer modn, i__, k;
    static real xhold, t1, t2;
    extern /* Subroutine */ int rfftf1_(integer *, real *, real *, real *, 
	    integer *);
    static integer kc, np1, ns2;

    /* Parameter adjustments */
    --ifac;
    --x;
    --xh;
    --was;
    --war;

    /* Function Body */
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	xh[i__] = war[i__];
	war[i__] = x[i__];
/* L100: */
    }
    if ((i__1 = *n - 2) < 0) {
	goto L101;
    } else if (i__1 == 0) {
	goto L102;
    } else {
	goto L103;
    }
L101:
    xh[1] += xh[1];
    goto L106;
L102:
    xhold = sqrt3 * (xh[1] + xh[2]);
    xh[2] = sqrt3 * (xh[1] - xh[2]);
    xh[1] = xhold;
    goto L106;
L103:
    np1 = *n + 1;
    ns2 = *n / 2;
    x[1] = 0.f;
    i__1 = ns2;
    for (k = 1; k <= i__1; ++k) {
	kc = np1 - k;
	t1 = xh[k] - xh[kc];
	t2 = was[k] * (xh[k] + xh[kc]);
	x[k + 1] = t1 + t2;
	x[kc + 1] = t2 - t1;
/* L104: */
    }
    modn = *n % 2;
    if (modn != 0) {
	x[ns2 + 2] = xh[ns2 + 1] * 4.f;
    }
    rfftf1_(&np1, &x[1], &xh[1], &war[1], &ifac[1]);
    xh[1] = x[1] * .5f;
    i__1 = *n;
    for (i__ = 3; i__ <= i__1; i__ += 2) {
	xh[i__ - 1] = -x[i__];
	xh[i__] = xh[i__ - 2] + x[i__ - 1];
/* L105: */
    }
    if (modn != 0) {
	goto L106;
    }
    xh[*n] = -x[*n + 1];
L106:
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	x[i__] = war[i__];
	war[i__] = xh[i__];
/* L107: */
    }
    return 0;
} /* sint1_ */

/* Subroutine */ int sinti_(integer *n, real *wsave)
{
    /* Initialized data */

    static real pi = 3.14159265358979f;

    /* System generated locals */
    integer i__1;

    /* Builtin functions */
    double sin(doublereal);

    /* Local variables */
    static integer k;
    extern /* Subroutine */ int rffti_(integer *, real *);
    static real dt;
    static integer np1, ns2;

    /* Parameter adjustments */
    --wsave;

    /* Function Body */
    if (*n <= 1) {
	return 0;
    }
    ns2 = *n / 2;
    np1 = *n + 1;
    dt = pi / (real) np1;
    i__1 = ns2;
    for (k = 1; k <= i__1; ++k) {
	wsave[k] = sin(k * dt) * 2.f;
/* L101: */
    }
    rffti_(&np1, &wsave[ns2 + 1]);
    return 0;
} /* sinti_ */

