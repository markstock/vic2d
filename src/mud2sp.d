c
c     file mud2sp.d
c
c     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c     *                                                               *
c     *                  copyright (c) 2008 by UCAR                   *
c     *                                                               *
c     *       University Corporation for Atmospheric Research         *
c     *                                                               *
c     *                      all rights reserved                      *
c     *                                                               *
c     *                     MUDPACK  version 5.0.1                    *
c     *                                                               *
c     *                 A Fortran Package of Multigrid                *
c     *                                                               *
c     *                Subroutines and Example Programs               *
c     *                                                               *
c     *      for Solving Elliptic Partial Differential Equations      *
c     *                                                               *
c     *                             by                                *
c     *                                                               *
c     *                         John Adams                            *
c     *                                                               *
c     *                             of                                *
c     *                                                               *
c     *         the National Center for Atmospheric Research          *
c     *                                                               *
c     *                Boulder, Colorado  (80307)  U.S.A.             *
c     *                                                               *
c     *                   which is sponsored by                       *
c     *                                                               *
c     *              the National Science Foundation                  *
c     *                                                               *
c     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c ... file mud2sp.d
c
c     contains documentation for:
c     subroutine mud2sp(iparm,fparm,work,cofx,cofy,bndyc,rhs,phi,mgopt,ierror)
c     A sample fortran driver is file "tmud2sp.f".
c
c ... required MUDPACK files
c
c     mudcom.f
c
c ... purpose
c
c     subroutine mud2sp automatically discretizes and attempts to compute
c     the second-order difference approximation to the two-dimensional
c     linear separable elliptic partial differential equation on a
c     rectangle.  the approximation is generated on a uniform grid covering
c     the rectangle (see mesh description below).  boundary conditions
c     may be specified (dirchlet), periodic, or mixed derivative in any
c     combination.  the form of the pde solved is:
c
c
c          cxx(x)*pxx + cx(x)*px + cex(x)*p(x,y) +
c
c          cyy(y)*pyy + cy(y)*py + cey(y)*p(x,y) = r(x,y)
c
c     pxx,pyy,px,py are second and first partial derivatives of the
c     unknown real solution function p(x,y) with respect to the
c     independent variables x,y.  cxx,cx,cex,cyy,cy,cey are the known
c     real coefficients of the elliptic pde and r(x,y) is the known
c     real right hand side of the equation.  cxx and cyy should be
c     positive for all x,y in the solution region.  If some of the
c     coefficients depend on both x and y then the PDE is nonseparable.
c     In this case subroutine muh2 or mud2 must be used instead of mud2sp
c     (see the files muh2.d or mud2.d)
c                                                                               
c
c ... mesh description . . .
c
c     the approximation is generated on a uniform nx by ny grid.  the grid
c     is superimposed on the rectangular solution region
c
c          [xa,xb] x [yc,yd].
c
c     let
c
c          dlx = (xb-xa)/(nx-1), dly = (yd-yc)/(ny-1)
c
c     be the uniform grid increments in the x,y directions. then
c
c          xi=xa+(i-1)*dlx,  yj=yc+(j-1)*dly
c
c     for i=1,...,nx and j=1,...,ny  denote the x,y uniform mesh points
c
c
c ... language                                                                  
c
c     fortran90/fortran77
c
c
c ... portability                                                               
c
c     mudpack5.0.1 software has been compiled and tested with Fortran77
c     and Fortran90 on a variety of platforms.
c                                                                               
c ... methods
c
c     details of the methods employeed by the solvers in mudpack are given
c     in [1,9].  [1,2,9] contain performance measurements on a variety of
c     elliptic pdes (see "references" in the file "readme").  in summary:
c
c *** discretization and solution (second-order solvers) (see [1])
c
c     the pde and boundary conditions are automatically discretized at all
c     grid levels using second-order finite difference formula.  diagonal
c     dominance at coarser grid levels is maintained in the presence of
c     nonzero first-order terms by adjusting the second-order coefficient
c     when necessary.  the resulting block tri-diagonal linear system is
c     approximated using multigrid iteration [10,11,13,15,16,18].  version
c     5.0.1 of mudpack uses only fully weighted residual restriction.  defaults
c     include cubic prolongation and w(2,1) cycles.  these can be overridden
c     with  selected multigrid options (see "mgopt").  error control based on
c     maximum relative differences is available. full multigrid cycling (fmg)
c     or cycling beginning or restarting at the finest grid level can be
c     selected. a menu of relaxation methods including gauss-seidel point,
c     line relaxation(s) (in any combination of directions) and planar
c     relaxation (for three-dimensional anisotropic problems) are provided.
c     all methods use ordering based on alternating points (red/black),
c     lines, or planes for cray vectorization and improved convergence
c     rates [14].
c
c *** higher order solution (fourth-order solvers) (see [9,19,21])
c
c     if the multigrid cycling results in a second-order estimate (i.e.,
c     discretization level error is reached) then this can be improved to a
c     fourth-order estimate using the technique of "deferred corrections."
c     the values in the solution array are used to generate a fourth-order
c     approximation to the truncation error.  second-order finite difference
c     formula are used to approximate third and fourth partial derivatives
c     of the solution function [3].  the truncation error estimate is
c     transferred down to all grid levels using weighted averaging where
c     it serves as a new right hand side.  the default multigrid options
c     are used to compute the fourth-order correction term which is added
c     to the original solution array.
c
c
c ... references (partial)
c
c
c     [1] J. Adams, "MUDPACK: Multigrid Fortran Software for the Efficient
c     Solution of Linear Elliptic Partial Differential Equations,"
c     Applied Math. and Comput. vol.34, Nov 1989, pp.113-146.
c
c     [2] J. Adams,"FMG Results with the Multigrid Software Package MUDPACK,"
c     proceedings of the fourth Copper Mountain Conference on Multigrid, SIAM,
c     1989, pp.1-12.
c     .
c     .
c     .
c     [7] J. Adams, R. Garcia, B. Gross, J. Hack, D. Haidvogel, and V. Pizzo,
c     "Applications of Multigrid Software in the Atmospheric Sciences,"
c      Mon. Wea. Rev.,vol. 120 # 7, July 1992, pp. 1447-1458.
c     .
c     .
c     .
c     [9] J. Adams, "Recent Enhancements in MUDPACK, a Multigrid Software
c     package for Elliptic Partial Differential Equations," Applied Math.
c     and Comp., 1991, vol. 43, May 1991, pp. 79-94.
c
c     [10]J. Adams, "MUDPACK-2: Multigrid Software for Approximating
c     Elliptic Partial Differential Equations on Uniform Grids with
c     any Resolution," Applied Math. and Comp., 1993, vol. 53, February
c     1993, pp. 235-249
c     .
c     .
c     .
c
c ... argument description
c                                                                               
c
c **********************************************************************
c *** input arguments *************************************************
c **********************************************************************
c
c
c ... iparm                                                                     
c
c          an integer vector of length 17 used to pass integer
c          arguments.  iparm is set internally and defined as
c          follows:
c
c
c ... intl=iparm(1)
c
c          an initialization argument.  intl=0  must be input
c          on an initial call. in this case input arguments will
c          be checked for errors and the elliptic partial differential
c          equation and boundary conditions will be discretized using
c          second order finite difference formula.
c
c ***      An approximation is NOT generated after an intl=0 call!
c          mud2sp should be called with intl=1 to approximate the elliptic
c          PDE discretized by the intl=0 call.  intl=1 should also
c          be input if mud2sp has been called earlier and only the
c          values in in rhs (see below) or gbdy (see bndyc below)
c          or phi (see below) have changed.  This will bypass
c          redundant pde discretization and argument checking
c          and save computational time.  Some examples of when
c          intl=1 calls should be used are:
c
c          (0) after a intl=0 argument checking and discretization call
c
c          (1) mud2sp is being recalled for additional accuracy.  In
c              this case iguess=iparm(12)=1 should also be used.
c
c          (2) mud2sp is being called every time step in a time dependent
c              problem (see discussion below) where the elliptic operator
c              does not depend on time.
c
c          (3) mud2sp is being used to solve the same elliptic equation
c              for several different right hand sides (iguess=0 should
c              probably be used for each new righthand side).
c
c          intl = 0 must be input before calling with intl = 1 when any
c          of the following conditions hold:
c
c          (a) the initial call to mud2sp

c          (b) any of the integer arguments other than iguess=iparm(12)
c              or maxcy=iparm(13) or mgopt have changed since the previous
c              call.
c
c          (c) any of the floating point arguments other than tolmax=
c              fparm(5) have changed since the previous call
c
c          (d) any of the coefficients input by cofx,cofy (see below) have
c              changed since the previous call
c
c          (e) any of the constant "alfa" coefficients input by bndyc
c              (see below) have changed since the previous call.
c
c          If any of (a) through (e) are true then the elliptic PDE
c          must be discretized or rediscretized.  If none of (a)
c          through (e) holds, calls can be made with intl=1.
c          Incorrect calls with intl=1 will produce erroneous results.
c  ***     The values set in the saved work space "work" (see below) with
c          an intl=0 call must be preserved with subsequent intl=1 calls.
c
c          MUDPACK software performance should be monitored for intl=1
c          calls.  The intl=0 discretization call performance depends
c          primarily on the efficiency or lack of efficiency of the
c          user provided subroutines for pde coefficients and
c          boundary conditions.
c
c ... nxa=iparm(2)
c
c          flags boundary conditions on the edge x=xa
c
c        = 0 if p(x,y) is periodic in x on [xa,xb]
c            (i.e., p(x+xb-xa,y) = p(x,y) for all x,y
c            (if nxa=0 then nxb=0 is required, see ierror = 2)
c
c        = 1 if p(xa,y) is specified (this must be input thru phi(1,j))
c
c        = 2 if there are mixed derivative boundary conditions at x=xa
c            (see bndyc)
c
c
c ... nxb=iparm(3)
c
c          flags boundary conditions on the edge x=xb
c
c        = 0 if p(x,y) is periodic in x on [xa,xb]
c            (i.e., p(x+xb-xa,y) = p(x,y) for all x,y)
c            (if nxb=0 then nxa=0 is required, see ierror = 2)
c
c        = 1 if p(xb,y) is specified (this must be input thru phi(nx,j))
c
c        = 2 if there are mixed derivative boundary conditions at x=xb
c            (see bndyc)
c
c
c ... nyc=iparm(4)
c
c          flags boundary conditions on the edge y=yc
c
c        = 0 if p(x,y) is periodic in y on [yc,yd]
c            (i.e., p(x,y+yd-yc) = p(x,y) for all x,y
c            (if nyc=0 then nyd=0 is required, see ierror = 2)
c
c        = 1 if p(x,yc) is specified (this must be input thru phi(i,1))
c
c        = 2 if there are mixed derivative boundary conditions at y=yc
c            (see bndyc)
c
c
c ... nyd=iparm(5)
c
c          flags boundary conditions on the edge y=yd
c
c        = 0 if p(x,y) is periodic in y on [yc,yd]
c            (i.e., p(x,y+yd-yc) = p(x,y) for all x,y
c             (if nyd=0 then nyc=0 is required, see ierror = 2)
c
c        = 1 if p(x,yd) is specified (this must be input thru phi(i,ny))
c
c        = 2 if there are mixed derivative boundary conditions at y=yd
c            (see bndyc)
c
c
c *** grid size arguments
c
c
c ... ixp = iparm(6)
c
c          an integer greater than one which is used in defining the number
c          of grid points in the x direction (see nx = iparm(10)).  "ixp+1"
c          is the number of points on the coarsest x grid visited during
c          multigrid cycling.  ixp should be chosen as small as possible.
c          recommended values are the small primes 2 or 3.
c          larger values can reduce multigrid convergence rates considerably,
c          especially if line relaxation in the x direction is not used.
c          if ixp > 2 then it should be 2 or a small odd value since a power
c          of 2 factor of ixp can be removed by increasing iex = iparm(8)
c          without changing nx = iparm(10).
c
c
c ... jyq = iparm(7)
c
c          an integer greater than one which is used in defining the number
c          of grid points in the y direction (see ny = iparm(11)).  "jyq+1"
c          is the number of points on the coarsest y grid visited during
c          multigrid cycling.  jyq should be chosen as small as possible.
c          recommended values are the small primes 2 or 3.
c          larger values can reduce multigrid convergence rates considerably,
c          especially if line relaxation in the y direction is not used.
c          if jyq > 2 then it should be 2 or a small odd value since a power
c          of 2 factor of jyq can be removed by increasing jey = iparm(9)
c          without changing ny = iparm(11).
c
c
c ... iex = iparm(8)
c
c          a positive integer exponent of 2 used in defining the number
c          of grid points in the x direction (see nx = iparm(10)).
c          iex .le. 50 is required.  for efficient multigrid cycling,
c          iex should be chosen as large as possible and ixp=iparm(8)
c          as small as possible within grid size constraints when
c          defining nx.
c
c
c ... jey = iparm(9)
c
c          a positive integer exponent of 2 used in defining the number
c          of grid points in the y direction (see ny = iparm(11)).
c          jey .le. 50 is required.  for efficient multigrid cycling,
c          jey should be chosen as large as possible and jyq=iparm(7)
c          as small as possible within grid size constraints when
c          defining ny.
c
c
c
c ... nx = iparm(10)
c
c          the number of equally spaced grid points in the interval [xa,xb]
c          (including the boundaries).  nx must have the form
c
c               nx = ixp*(2**(iex-1)) + 1
c
c          where ixp = iparm(6), iex = iparm(8).
c
c
c ... ny = iparm(11)
c
c          the number of equally spaced grid points in the interval [yc,yd]
c          (including the boundaries).  ny must have the form:
c
c               ny = jyq*(2**(jey-1)) + 1
c
c          where jyq = iparm(7), jey = iparm(9).
c
c
c *** example
c
c         suppose a solution is wanted on a 33 by 97 grid.  then
c         ixp=2, jyq=6 and iex=jey=5 could be used.  a better
c         choice would be ixp=2, jyq=3, and iex=5, jey=6.
c
c
c *** note
c
c     let G be the nx by ny fine grid on which the approximation is
c     generated and let n = max0(iex,jey).  in mudpack, multigrid
c     cycling is implemented on the ascending chain of grids
c
c         G(1) < ... < G(k) < ... < G(n) = G.
c
c     each G(k) (k=1,...,n) has mx(k) by my(k) grid points
c     given by:
c
c         mx(k) = ixp*[2**(max0(iex+k-n,1)-1)] + 1
c
c         my(k) = jyq*[2**(max0(jey+k-n,1)-1)] + 1
c
c
c
c ... iguess=iparm(12)
c
c          = 0 if no initial guess to the pde is provided
c
c          = 1 if an initial guess to the pde is at the finest grid
c              level is provided in phi (see below)
c
c     comments on iguess = 0 or 1 . . .
c
c     even if iguess = 0, phi must be initialized at all grid points (this
c     is not checked).  phi can be set to 0.0 at non-dirchlet grid points
c     if nothing better is available.  the values set in phi when iguess = 0
c     are passed down and serve as an initial guess to the pde at the coarsest
c     grid level where cycling commences.  in this sense, values input in
c     phi always serve as an initial guess.  setting iguess = 0 forces full
c     multigrid cycling beginning at the coarsest and finishing at the finest
c     grid level.
c
c     if iguess = 1 then the values input in phi are an initial guess to the
c     pde at the finest grid level where cycling begins.  this option should
c     be used only if a "very good" initial guess is available (as, for
c     example, when restarting from a previous iguess=0 call).
c
c     time dependent problems . . .
c
c *** assume we are solving an elliptic pde every time step in a
c     marching problem of the form:
c
c          l(p(t)) = r(t)
c
c     where the differential operator "l" has no time dependence,
c     "p(t)" is the solution and "r(t)" is the right hand side at
c     current time "t". let "dt" be the increment between time steps.
c     then p(t) can be used as an initial guess to p(t+dt) with
c     intl = 1 when solving
c
c          l(p(t+dt)) = r(t+dt).
c
c     after the first two time steps, rather than continue, it would
c     be better to define the "correction" term:
c
c          e(t,dt) = p(t+dt) - p(t)
c
c     this clearly satisfies the equation
c
c          l(e(t,dt)) = r(t+dt) - r(t).
c
c     this should be solved with iguess = 0 and intl = 1. boundary
c     conditions for e(t,dt) are obtained from the boundary conditions
c     for p(t) by subtracting given values at t from given values at
c     t+dt. for example if
c
c          d(p(t))/dx = f(t), d(p(t+dt))/dx = f(t+dt)
c
c     at some x boundary then e(t,dt) satisfies the derivative
c     boundary condition
c
c          d(e(t,dt))/dx = f(t+dt) - f(t).
c
c     e(t,dt) can be preset to 0.0 (at nondirchlet points) or (if p(t-dt)
c     is saved) to p(t)-p(t-dt).  with iguess = 0, these values will serve
c     as an initial guess to e(t,dt) at the coarsest grid level.  this
c     approach has the advantage that a full sequence of multigrid cycles,
c     beginning at the coarsest grid level, is invoked every time step in
c     solving for e(t,dt).  a few digits of accuracy in e(t,dt), which is
c     ordinarily much smaller than p(t), will yield several more digits of
c     accuracy in the final approximation:
c
c          p(t+dt) = p(t) + e(t,dt).
c
c     using this approach to integrate in time will give more accuracy
c     then using p(t) as an initial guess to p(t+dt) for all time steps.
c     it does require additional storage.
c
c     if the differential operator "l" has time dependence (either thru
c     the coefficients in the pde or the coefficients in the derivative
c     boundary conditions) then use p(t) as an initial guess to p(t+dt)
c     when solving
c
c          l(t+dt)(p(t+dt)) = r(t+dt)
c
c     with intl = 0 for all time steps (the discretization must be repeated
c     for each new "t"). either iguess = 0 (p(t) will then be an initial
c     guess at the coarsest grid level where cycles will commence) or
c     iguess = 1 (p(t) will then be an initial guess at the finest grid
c     level where cycles will remain fixed) can be tried.
c
c
c ... maxcy = iparm(13)
c
c          the exact number of cycles executed between the finest (nx by
c          ny) and the coarsest ((ixp+1) by (jyq+1)) grid levels when
c          tolmax=fparm(5)=0.0 (no error control).  when tolmax > 0.0
c          is input (error control) then maxcy is a limit on the number
c          of cycles between the finest and coarsest grid levels.  in
c          any case, at most maxcy*(iprer+ipost) relaxation sweeps are
c          are performed at the finest grid level (see iprer=mgopt(2),
c          ipost=mgopt(3) below).  when multigrid iteration is working
c          "correctly" only a few are required for convergence.  large
c          values for maxcy should not be necessary.
c
c
c ... method = iparm(14) determines the method of relaxation
c              (gauss-seidel based on alternating points or lines)
c
c          = 0 for  point relaxation
c
c          = 1 for line relaxation in the x direction
c
c          = 2 for line relaxation in the y direction
c
c          = 3 for line relaxation in both the x and y direction
c
c
c *** choice of method. . .
c
c     let fx represent the quantity cxx(x,y)/dlx**2 over the solution region.
c
c     let fy represent the quantity cyy(x,y)/dly**2 over the solution region
c
c     if fx,fy are roughly the same size and do not vary too much over
c     the solution region choose method = 0.  if this fails try method=3.
c
c     if fx is much greater than fy choose method = 1.
c
c     if fy is much greater than fx choose method = 2
c
c     if neither fx or fy dominates over the solution region and they
c     both vary considerably choose method = 3.
c
c
c ... length = iparm(15)
c
c          the length of the work space provided in vector work (see below).
c          let isx = 0 if method = 0 or method = 2
c          let isx = 3 if method = 1 or method = 3 and nxa.ne.0
c          let isx = 5 if method = 1 or method = 3 and nxa.eq.0
c          let jsy = 0 if method = 0 or method = 1
c          let jsy = 3 if method = 2 or method = 3 and nyc.ne.0
c          let jsy = 5 if method = 2 or method = 3 and nyc.eq.0
c          then . . .
c
c               length = nx*ny*(5+3*(isx+jsy)/2)+ 10*(nx+ny)
c
c          will suffice in all cases but very small nx and ny.
c          the exact minimal work space length required for the
c          current set of input arugments is output in iparm(16).
c          (even if iparm(15) is too small).  this will be usually
c          be less then the value given by the simplified formula
c          above.  * Notice that mud2sp requires considerably less
c          work space than the nonseparable solvers muh2,mud2 if
c          and only if method=0 is chosen.
c
c ... fparm                                                                     
c
c          a floating point vector of length 6 used to efficiently
c          pass floating point arguments.  fparm is set internally
c          in mud2sp and defined as follows . . .
c
c
c ... xa=fparm(1), xb=fparm(2)
c
c          the range of the x independent variable. xa must                     
c          be less than xb                                                      
c
c
c ... yc=fparm(3), yd=fparm(4)
c
c          the range of the y independent variable.  yc must                    
c          be less than yd.                                                     
c
c
c ... tolmax = fparm(5)
c
c          when input positive, tolmax is a maximum relative error tolerance
c          used to terminate the relaxation iterations. assume phi1(i,j)
c          and phi2(i,j) are the last two computed approximations at all
c          grid points of the finest grid level. if we define
c
c              phdif = max(abs(phi2(i,j)-phi1(i,j))) for all i,j
c
c          and
c
c              phmax = max(abs(phi2(i,j))) for all i,j
c
c          then "convergence" is considered to have occurred if and only if
c
c              phdif/phmax < tolmax.
c
c
c          if tolmax=fparm(5)=0.0 is input then there is no error control
c          and maxcy cycles from the finest grid level are executed. maxcy
c          is a limit which cannot be exceeded even with error control.
c     ***  calls with tolmax=0.0, when appropriate because of known
c          convergence behavior, are more efficient than calls with tolmax
c          positive (i.e., if possible DO NOT use error control!).
c
c ... work                                                                      
c
c          a one dimensional real saved work space (see iparm(15) for
c          length) which must be preserved from the previous call when
c          calling with intl=iparm(1)=1.
c
c ... bndyc                                                                     
c
c          a subroutine with  arguments (kbdy,xory,alfa,gbdy) which
c          are used to input mixed boundary conditions to mud2sp. bndyc
c          must be declared "external" in the program calling mud2sp.
c          the boundaries are numbered one thru four and the mixed
c          derivative boundary conditions are described below (see the
c          sample driver code "tmud2sp.f" for an example of how bndyc is
c          can beset up).
c                                                                               
c          * * * * * * * * * * * *  y=yd
c          *       kbdy=4        *
c          *                     *
c          *                     *
c          *                     *
c          * kbdy=1       kbdy=2 *
c          *                     *
c          *                     *
c          *                     *
c          *       kbdy=3        *
c          * * * * * * * * * * * *  y=yc
c
c          x=xa                  x=xb
c
c
c          (1) the kbdy=1 boundary
c
c          this is the edge x=xa where nxa=iparm(2)=2 flags
c          a mixed boundary condition of the form
c
c             dp/dx + alfxa*p(xa,y) = gbdxa(y)
c
c          in this case kbdy=1,xory=y will be input to bndyc and
c          alfa,gbdy corresponding to alfxa,gbdxa(y) must be returned
c
c
c          (2) the kbdy=2 boundary
c
c          this is the edge x=xb where nxb=iparm(3)=2 flags
c          a mixed boundary condition of the form
c
c             dp/dx + alfxb*p(xb,y) = gbdxb(y)
c
c          in this case kbdy=2,xory=y, will be input to bndyc and
c          alfa,gbdy corresponding to alfxb,gbdxb(y) must be returned.
c
c
c          (3) the kbdy=3 boundary
c
c          this is the edge y=yc where nyc=iparm(4)=2 flags
c          a mixed boundary condition of the form
c
c             dp/dy + alfyc*p(x,yc) = gbdyc(x)
c
c          in this case kbdy=3,xory=x will be input to bndyc and
c          alfa,gbdy corresponding to alfyc,gbdyc(x) must be returned.
c
c
c          (4) the kbdy=4 boundary
c
c          this is the edge y=yd where nyd=iparm(5)=2 flags
c          a mixed boundary condition of the form
c
c             dp/dy + alfyd*p(x,yd) = gbdyd(x)
c
c          in this case kbdy=4,xory=x will be input to bndyc and
c          alfa,gbdy corresponding to alfyd,gbdyd(x) must be returned.
c
c                                                                               
c ***      alfxa,alfxb,alfyc,alfyd must be constants for mud2sp.
c          Use muh2 or mud2 if any of these depend on x or y.
c          bndyc must provide the mixed boundary condition values
c          in correspondence with those flagged in iparm(2) thru
c          iparm(5).  if all boundaries are specified or periodic
c          mud2sp will never call bndyc.  even then it must be entered
c          as a dummy subroutine. bndyc must be declared "external"
c          in the routine calling mud2sp the actual name chosen may
c          be different.
c
c
c ... cofx
c
c         a subroutine with arguments (x,cxx,cx,cex) which provides
c         the known real x dependent coefficients for the separable
c         elliptic pde at any x grid point.  the name chosen in the calling
c         routine may be different where the coefficient routine must be declared
c         "external."
c
c ... cofy
c
c         a subroutine with arguments (y,cyy,cy,cey) which provides
c         the known real y dependent coefficients for the separable
c         elliptic pde at any y grid point.  the name chosen in the calling
c         routine may be different where the coefficient routine must be declared
c         "external."
c
c ... rhs                                                                       
c
c          an array dimensioned nx by ny which contains the given
c          right hand side values on the uniform 2-d mesh.
c
c              rhs(i,j) = r(xi,yj) for i=1,...,nx and j=1,...,ny
c
c ... phi                                                                       
c
c          an array dimensioned nx by ny.  on input phi must contain
c          specified boundary values. for example, if nyd=iparm(5)=1
c          then phi(i,ny) must be set equal to p(xi,yd) for i=1,...nx
c          prior to calling mud2sp.  these values are preserved by mud2sp.
c          if an initial guess is provided (iguess=iparm(11)=1) it must
c          be input thru phi.
c
c
c ***      if no initial guess is given (iguess=0) then phi must still
c          be initialized at all grid points (this is not checked).  these
c          values will serve as an initial guess to the pde at the coarsest
c          grid level after a transfer from the fine solution grid.  set phi
c          equal to to 0.0 at all internal and non-specified boundaries
c          grid points if nothing better is available.
c
c
c ... mgopt
c
c           an integer vector of length 4 which allows the user to select
c           among various multigrid options.  if mgopt(1)=0 is input then
c           a default set of multigrid arguments (chosen for robustness)
c           will be internally selected and the remaining values in mgopt
c           will be ignored.  if mgopt(1) is nonzero then the arguments
c           in mgopt are set internally and defined as follows:  (see the
c           basic coarse grid correction algorithm below)
c
c
c     kcycle = mgopt(1)
c
c            = 0 if default multigrid options are to be used
c
c            = 1 if v cycling is to be used (the least expensive per cycle)
c
c            = 2 if w cycling is to be used (the default)
c
c            > 2 if more general k cycling is to be used
c             *** warning--values larger than 2 increase
c                 the execution time per cycle considerably and
c                 result in the nonfatal error ierror = -5
c                 which indicates inefficient multigrid cycling.
c
c    iprer = mgopt(2)
c
c           the number of "pre-relaxation" sweeps executed before the
c           residual is restricted and cycling is invoked at the next
c           coarser grid level (default value is 2 whenever mgopt(1)=0)
c
c    ipost = mgopt(3)
c
c           the number of "post relaxation" sweeps executed after cycling
c           has been invoked at the next coarser grid level and the residual
c           correction has been transferred back (default value is 1
c           whenever mgopt(1)=0).
c
c *** if iprer, ipost, or (especially) kcycle is greater than 2
c     than inefficient multigrid cycling has probably been chosen and
c     the nonfatal error (see below) ierror = -5 will be set.  note
c     this warning may be overridden by any other nonzero value
c     for ierror.
c
c   intpol = mgopt(4)
c
c          = 1 if multilinear prolongation (interpolation) is used to
c              transfer residual corrections and the pde approximation
c              from coarse to fine grids within full multigrid cycling.
c
c          = 3 if multicubic prolongation (interpolation) is used to
c              transfer residual corrections and the pde approximation
c              from coarse to fine grids within full multigrid cycling.
c              (this is the default value whenever mgopt(1)=0).
c
c *** the default values (2,2,1,3) in the vector mgopt were chosen for
c     robustness.  in some cases v(2,1) cycles with linear prolongation will
c     give good results with less computation (especially in two-dimensions).
c     this  was the default and only choice in an earlier version of mudpack
c     (see [1]) and can be set with the integer vector (1,2,1,1) in mgopt.
c
c *** the schedules for one full multigrid cycle (iguess=0) using v(2,1)
c     cycles and w(2,1) cycles are depicted for a four level grid below.
c     the number of relaxation sweeps when each grid is visited are indicated.
c     the "*" stands for prolongation of the full approximation and the "."
c     stands for transfer of residuals and residual corrections within the
c     coarse grid correction algorithm (see below).  all version 5.0.1
c     mudpack solvers use only fully weighted residual restriction
c
c     one fmg with v(2,1) cycles:
c
c
c     ------------------------------2-----------------1------     level 4
c                                  * .               .
c                                 *   .             .
c     ---------------2-----------1-----2-----------1---------     level 3
c                   * .         .       .         .
c                  *   .       .         .       .
c     ------2-----1-----2-----1-----------2-----1------------     level 2
c          * .   .       .   .             .   .
c         *   . .         . .               . .
c     ---3-----3-----------3-----------------3---------------     level 1
c
c
c     one fmg with w(2,1) cycles:
c
c     ------------------------2---------------------------1--     level 4
c                            * .                         .
c     ----------2-----------1---2-----------3-----------1----     level 3
c              * .         .     .         . .         .
c     ----2---1---2---3---1-------2---3---1---2---3---1------     level 2
c        * . .     . . . .         . . . .     . . . .
c     --6---6-------6---6-----------6---6-------6---6--------     level 1
c
c
c     the form of the "recursive" coarse grid correction cycling used
c     when kcycle.ge.0 is input is described below in pseudo-algorithmic
c     language.  it is implemented non-recursively in fortran in mudpack.
c
c     algorithm cgc(k,l(k),u(k),r(k),kcycle,iprer,ipost,iresw,intpol)
c
c *** approximately solve l(k)*u(k) = r(k) using multigrid iteration
c *** k is the current grid level
c *** l(k) is the discretized pde operator at level k
c *** u(k) is the initial guess at level k
c *** r(k) is the right hand side at level k
c *** i(k,k-1) is the restriction operator from level k to level k-1
c *** (the form of i(k,k-1) depends on iresw)
c *** i(k-1,k) is the prolongation operator from level k-1 to level k
c *** (the form of i(k-1,k) depends on intpol)
c
c     begin algorithm cgc
c
c ***   pre-relax at level k
c
c     . do (i=1,iprer)
c
c     . .  relax(l(k),u(k),r(k))
c
c     . end do
c
c     . if (k > 1) then
c
c ***     restrict the residual from level k to level k-1
c
c     . . r(k-1) = i(k,k-1)(r(k)-l(k)*u(k))
c
c     . . kount = 0
c
c     . . repeat
c
c ***     solve for the residual correction at level k-1 in u(k-1)
c ***     using algorithm cgc "kcycle" times (this is the recursion)
c
c     . . . kount = kount+1
c
c     . . . invoke cgc(k-1,l(k-1),u(k-1),r(k-1),kcycle,iprer,ipost,iresw)
c
c
c     . . until (kount.eq.kcycle)
c
c ***     transfer residual correction in u(k-1) to level k
c ***     with the prolongation operator and add to u(k)
c
c     . . u(k) = u(k) + i(k-1,k)(u(k-1))
c
c     . end if
c
c ***   post relax at level k
c
c     . do (i=1,ipost)
c
c     . . relax(l(k),u(k),r(k))
c
c     . end do
c
c     . return
c
c     end algorithm cgc
c
c
c **********************************************************************
c *** output arguments ************************************************
c **********************************************************************
c
c
c ... iparm(16)  *** set for intl=0 calls only
c
c          on output iparm(16) contains the actual work space length
c          required.  this will usually be less than that given by the
c          simplified formula for length=iparm(15) (see as input argument)
c
c
c ... iparm(17)  *** set for intl=1 calls only
c
c          on output iparm(17) contains the actual number of multigrid cycles
c          between the finest and coarsest grid levels used to obtain the
c          approximation when error control (tolmax > 0.0) is set.
c
c
c ... fparm(6)   *** set for intl=1 calls with fparm(5) > 0. only
c
c          on output fparm(6) contains the final computed maximum relative
c          difference between the last two iterates at the finest grid level.
c          fparm(6) is computed only if there is error control (tolmax > 0.0)
c          assume phi1(i,j,k) and phi2(i,j,k) are the last two computed
c          values for phi(i,j,k) at all points of the finest grid level.
c          if we define
c
c               phdif = max(abs(phi2(i,j)-phi1(i,j))) over all i,j
c
c          and
c
c               phmax = max(abs(phi2(i,j)) over all i,j
c
c          then
c
c               fparm(6) = phdif/phmax
c
c          is returned whenever phmax > 0.0. in the degenerate case
c          phmax = 0.0, fparm(6) = phdif is returned.
c
c
c ... work
c
c          on output work contains intermediate values that must not
c          be destroyed if mud2sp is to be called again with intl=1
c
c
c ... phi   *** for intl=1 calls only
c
c          on output phi(i,j) contains the approximation to p(xi,yj)
c          for all mesh points i = 1,...,nx and j=1,...,ny.  the last
c          computed iterate in phi is returned even if convergence is
c          not obtained
c
c
c ... ierror
c
c          For intl=iparm(1)=0 initialization calls, ierror is an
c          error flag that indicates invalid input arguments when
c          returned positive and nonfatal warnings when returned
c          negative.  Argument checking and discretization
c          is bypassed for intl=1 calls which can only return
c          ierror = -1 or 0 or 1.
c
c
c     non-fatal warnings * * *
c
c
c     =-5 if kcycle=mgopt(1) is greater than 2. values larger than 2 results
c         in an algorithm which probably does far more computation than
c         necessary.  kcycle = 1 (v cycles)  or kcycle=2 (w cycles) should
c         suffice for most problems.  ierror = -5 is also set if either
c         iprer = mgopt(2) or ipost=mgopt(3) is greater than 2.  the
c         ierror=-5 flag is overridden by any other fatal or non-fatal
c         error.
c
c     =-4 if there are dominant nonzero first order terms in the pde which
c         make it "hyperbolic" at the finest grid level. numerically, this
c         happens if:
c
c              abs(cx)*dlx > 2.*abs(cxx)   (dlx = (xb-xa)/(nx-1))
c
c                         (or)
c
c              abs(cy)*dly > 2.*abs(cyy)   (dly = (yd-yc)/(ny-1))
c
c
c         at some fine grid point (xi,yj).  if an adjustment is not made the
c         condition can lead to a matrix coming from the discretization
c         which is not diagonally dominant and divergence is possible. since
c         the condition is "likely" at coarser grid levels for pde's with
c         nonzero first order terms, the adjustments (actually first order
c         approximations to the pde)
c
c
c             cxx = amax1(cxx,0.5*abs(cx)*dx)
c
c                          (and)
c
c             cyy = amax1(cyy,0.5*abs(cy)*dy)
c
c
c         (here dx,dy are the x,y mesh sizes of the subgrid)
c
c         are made to preserve convergence of multigrid iteration. if made
c         at the finest grid level, it can lead to convergence to an
c         erroneous solution (flagged by ierror = -4).  a possible remedy
c         is to increase resolution. the ierror = -4 flag overrides the
c         nonfatal ierror = -5 flag.
c
c
c     =-3  if the continuous elliptic pde is singular.  this means the
c          boundary conditions are periodic or pure derivative at all
c          boundaries and ce(x,y) = 0.0 for all x,y.  a solution is still
c          attempted but convergence may not occur due to ill-conditioning
c          of the linear system coming from the discretization.  the
c          ierror = -3 flag overrides the ierror=-4 and ierror=-5 nonfatal
c          flags.
c
c
c     =-2  if the pde is not elliptic (i.e., cxx*cyy.le.0.0 for some (xi,yj))
c          in this case a solution is still attempted although convergence
c          may not occur due to ill-conditioning of the linear system.
c          the ierror = -2 flag overrides the ierror=-5,-4,-3 nonfatal
c          flags.
c
c
c     =-1  if convergence to the tolerance specified in tolmax=fparm(5)>0.
c          is not obtained in maxcy=iparm(13) multigrid cycles between the
c          coarsest (ixp+1,jyq+1) and finest (nx,ny) grid levels.
c          in this case the last computed iterate is still returned.
c          the ierror = -1 flag overrides all other nonfatal flags
c
c
c     no errors * * *
c
c     = 0
c
c     fatal argument errors * * *
c
c     = 1 if intl=iparm(1) is not 0 on the first call or not 0 or 1
c         on subsequent calls
c
c     = 2 if any of the boundary condition flags nxa,nxb,nyc,nyd
c         in iparm(2),iparm(3),iparm(4),iparm(5) are not 0,1 or 2
c         or if nxa,nxb or nyc,nyd are not pairwise zero.
c
c     = 3 if mino(ixp,jyq) < 2 (ixp = iparm(6), jyq = iparm(7))
c
c     = 4 if min0(iex,jey) < 1 (iex = iparm(8), jey = iparm(9)) or
c         if max0(iex,jey) > 50
c
c     = 5 if nx.ne.ixp*2**(iex-1)+1 or ny.ne.jyq*2**(jey-1)+1
c         (nx = iparm(10), ny = iparm(11))
c
c     = 6 if iguess = iparm(12) is not equal to 0 or 1
c
c     = 7 if maxcy = iparm(13) < 1
c
c     = 8 if method = iparm(14) is not 0,1,2, or 3
c
c     = 9 if length = iparm(15) is too small (see iparm(16) on output
c         for minimum required work space length)
c
c     =10 if xa >= xb or yc >= yd
c         (xa=fparm(1),xb=fparm(2),yc=fparm(3),yd=fparm(4))
c
c     =11 if tolmax = fparm(5) < 0.0
c
c     errors in setting multigrid options * * * (see also ierror=-5)
c
c     =12 if kcycle = mgopt(1) < 0 or
c         if iprer = mgopt(1) < 1 or
c         if ipost = mgopt(3) < 1 or
c         if intpol = mgopt(4) is not 1 or 3
c
c *********************************************************
c *********************************************************
c
c     end of mud2sp documentation
c
c **********************************************************
c **********************************************************
c
c
