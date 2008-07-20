c
c     file mud2sp.f
c
c  . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
c  .                                                             .
c  .                  copyright (c) 1999 by UCAR                 .
c  .                                                             .
c  .       UNIVERSITY CORPORATION for ATMOSPHERIC RESEARCH       .
c  .                                                             .
c  .                      all rights reserved                    .
c  .                                                             .
c  .                                                             .
c  .                      MUDPACK version 5.0                    .
c  .                                                             .
c  . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
c
c
c
c ... author and specialist
c
c          John C. Adams (National Center for Atmospheric Research)
c          email: johnad@ucar.edu, phone: 303-497-1213

c ... For MUDPACK 5.0 information, visit the website:
c     (http://www.scd.ucar.edu/css/software/mudpack)
c
c ... purpose (see mud2sp.d for details)
c
c     mud2sp attempts to produce a second order finite difference
c     approximation to the two dimensional separable elliptic
c     partial differential equation of the form:
c
c       cxx(x)*pxx + cx(x)*px + cex(x)*p(x,y) +
c
c       cyy(y)*pyy + cy(y)*py + cey(y)*p(x,y) = r(x,y)
c
c
c ... documentation and test files
c
c     see the documentation file "mud2sp.d" for a complete discussion
c     of how to use subroutine mud2sp.  file "tmud2sp.f" is a test/driver
c     sample program illustrating use of mud2sp
c
c ... required MUDPACK files
c
c     mudcom.f
c
c
      subroutine mud2sp(iparm,fparm,work,cfx,cfy,bndyc,rhs,phi,
     +                  mgopt,ierror)
      implicit none
      integer iparm,mgopt,ierror
      integer intl,nxa,nxb,nyc,nyd,ixp,jyq,iex,jey,nfx,nfy,iguess,
     +             maxcy,method,nwork,lwork,itero,ngrid,klevel,kcur,
     +             kcycle,iprer,ipost,intpol,kps
      real fparm,xa,xb,yc,yd,tolmax,relmax
      integer kpbgn,krbgn,kcxbgn,kcybgn,ktxbgn,ktybgn,nxk,nyk,isx,jsy
      integer int,iw,k,kb,nx,ny,itx,ity,icx,icy
      dimension iparm(17),fparm(6),mgopt(4)
      real work(*),phi(*),rhs(*)
      common/imud2/intl,nxa,nxb,nyc,nyd,ixp,jyq,iex,jey,nfx,nfy,iguess,
     +             maxcy,method,nwork,lwork,itero,ngrid,klevel,kcur,
     +             kcycle,iprer,ipost,intpol,kps
      common/fmud2/xa,xb,yc,yd,tolmax,relmax
      common/mud2spc/kpbgn(50),krbgn(50),kcxbgn(50),kcybgn(50),
     +ktxbgn(50),ktybgn(50),nxk(50),nyk(50),isx,jsy
      external cfx,cfy,bndyc
      data int / 0 /
      save int

c     first executable statement
c     write(6,*) "Inside mud2sp"
c     call flush(6)

      ierror = 1
      intl = iparm(1)    ! set and check intl on all calls
      if (intl*(intl-1).ne.0) return
      if (int.eq.0) then
	int = 1
	if (intl.ne.0) return  ! very first call is not intl=0
      end if
      ierror = 0
c
c     set  arguments internally
c     these will not be rechecked if intl=1!
c
      nxa = iparm(2)
      nxb = iparm(3)
      nyc = iparm(4)
      nyd = iparm(5)
      ixp = iparm(6)
      jyq = iparm(7)
      iex = iparm(8)
      jey = iparm(9)
      ngrid = max0(iex,jey)
      nfx = iparm(10)
      nfy = iparm(11)
      iguess = iparm(12)
      maxcy = iparm(13)
      method = iparm(14)
      nwork = iparm(15)
      kcycle = mgopt(1)
      if (kcycle .eq. 0) then
c       set defaults
	kcycle = 2
	iprer = 2
	ipost = 1
	intpol = 3
      else
	iprer = mgopt(2)
	ipost = mgopt(3)
	intpol = mgopt(4)
      end if
      xa = fparm(1)
      xb = fparm(2)
      yc = fparm(3)
      yd = fparm(4)
      tolmax = fparm(5)
      if (intl .eq. 0) then  ! intialization call
c
c     check input arguments
c
	ierror = 2   ! check boundary condition flags
	if (max0(nxa,nxb,nyc,nyd).gt.2) return
	if (min0(nxa,nxb,nyc,nyd).lt.0) return
	if (nxa.eq.0.and.nxb.ne.0) return
	if (nxa.ne.0.and.nxb.eq.0) return
	if (nyc.eq.0.and.nyd.ne.0) return
	if (nyc.ne.0.and.nyd.eq.0) return
	ierror = 3   ! check grid sizes
	if (ixp.lt.2) return
	if (jyq.lt.2) return
	ierror = 4
	ngrid = max0(iex,jey)
	if (iex.lt.1) return
	if (jey.lt.1) return
	if (ngrid.gt.50) return

	ierror = 5
        ! do pow (**) manually because macosx gfortran balks on it
c       kb = 1
c       do k=1,iex-1
c         kb = kb*2
c       enddo
	if (nfx.ne.ixp*2**(iex-1)+1) return
c       if (nfx.ne.ixp*kb+1) return

c       kb = 1
c       do k=1,jey-1
c         kb = kb*2
c       enddo
        if (nfy.ne.jyq*2**(jey-1)+1) return
c       if (nfy.ne.jyq*jey+1) return

	ierror = 6
	if (iguess*(iguess-1).ne.0) return
	ierror = 7
	if (maxcy.lt.1) return
	ierror = 8
	if (method.lt.0 .or. method.gt.3) return
	ierror = 9
c       compute and test minimum work space
	isx = 0
	if (method.eq.1 .or. method.eq.3) then
	  if (nxa.ne.0) isx = 3
	  if (nxa.eq.0) isx = 5
	end if
	jsy = 0
	if (method.eq.2 .or. method.eq.3) then
	  if (nyc.ne.0) jsy = 3
	  if (nyc.eq.0) jsy = 5
	end if
	kps = 1
	do k=1,ngrid
c       set subgrid sizes
	  nxk(k) = ixp*2**(max0(k+iex-ngrid,1)-1)+1
	  nyk(k) = jyq*2**(max0(k+jey-ngrid,1)-1)+1
	  nx = nxk(k)
	  ny = nyk(k)
	  kps = kps+(nx+2)*(ny+2)+(1+isx+jsy)*nx*ny+3*(nx+ny)
	end do
	iparm(16) = kps+(nfx+2)*(nfy+2)   ! exact minimum work space
	lwork = iparm(16)
	if (lwork .gt. nwork) return
	ierror = 10   ! check solution region
	if (xb.le.xa .or. yd.le.yc) return
	ierror = 11
	if (tolmax .lt. 0.0) return
	ierror = 12   ! multigrid parameters
	if (kcycle.lt.0) return
	if (min0(iprer,ipost).lt.1) return
	if ((intpol-1)*(intpol-3).ne.0) return
	if (max0(kcycle,iprer,ipost).gt.2) then
	  ierror = -5   ! inefficient multigrid cycling
	end if
	if (ierror .gt. 0) ierror = 0   ! no fatal errors
c
c     set work space pointers and discretize pde at each grid level
c
	iw = 1
	do kb=1,ngrid
	  k = ngrid-kb+1
	  nx = nxk(k)
	  ny = nyk(k)
	  kpbgn(k) = iw
	  krbgn(k) = kpbgn(k)+(nx+2)*(ny+2)
	  kcxbgn(k) = krbgn(k)+nx*ny
	  kcybgn(k) = kcxbgn(k)+3*nx
	  ktxbgn(k) = kcybgn(k)+3*ny
	  ktybgn(k) = ktxbgn(k)+isx*nx*ny
	  iw = ktybgn(k)+jsy*nx*ny
	  icx = kcxbgn(k)
	  icy = kcybgn(k)
	  itx = ktxbgn(k)
	  ity = ktybgn(k)
	  klevel = k
	  call dismd2sp(nx,ny,work(icx),work(icy),work(itx),work(ity),
     +                  bndyc,cfx,cfy,work,ierror)
	  end do
	return
      end if   ! end of intl=0 initialization call block
      nx = nfx
      ny = nfy
      call mud2sp1(nx,ny,rhs,phi,cfx,cfy,bndyc,work)
      iparm(17) = itero
      if (tolmax.gt.0.0) then   ! check for convergence
	fparm(6) = relmax
	if (relmax.gt.tolmax) ierror = -1   ! flag convergenc failure
      end if
      return
      end

      subroutine mud2sp1(nx,ny,rhsf,phif,cfx,cfy,bndyc,wk)
      implicit none
      integer nx,ny
      real phif(nx,ny),rhsf(nx,ny),wk(*)
      integer intl,nxa,nxb,nyc,nyd,ixp,jyq,iex,jey,nfx,nfy,iguess,
     +             maxcy,method,nwork,lwork,itero,ngrid,klevel,kcur,
     +             kcycle,iprer,ipost,intpol,kps
      real xa,xb,yc,yd,tolmax,relmax,phmax
      integer kpbgn,krbgn,kcxbgn,kcybgn,ktxbgn,ktybgn,nxk,nyk,isx,jsy
      integer k,kb,ip,ipc,ir,irc
      integer ncx,ncy,jj,ij,i,j,iter
      common/imud2/intl,nxa,nxb,nyc,nyd,ixp,jyq,iex,jey,nfx,nfy,iguess,
     +             maxcy,method,nwork,lwork,itero,ngrid,klevel,kcur,
     +             kcycle,iprer,ipost,intpol,kps
      common/fmud2/xa,xb,yc,yd,tolmax,relmax
      common/mud2spc/kpbgn(50),krbgn(50),kcxbgn(50),kcybgn(50),
     +ktxbgn(50),ktybgn(50),nxk(50),nyk(50),isx,jsy
      external cfx,cfy,bndyc
      nx = nxk(ngrid)
      ny = nyk(ngrid)
      ip = kpbgn(ngrid)
      ir = krbgn(ngrid)
c
c     set phif,rhsf in wk and adjust right hand side
c
      call swk2(nx,ny,phif,rhsf,wk(ip),wk(ir))
      if (iguess.eq.0) then
c
c     no initial guess at finest grid level!
c
	do kb=2,ngrid
	  k = ngrid-kb+1
	  nx = nxk(k+1)
	  ny = nyk(k+1)
	  ip = kpbgn(k+1)
	  ir = krbgn(k+1)
	  ncx = nxk(k)
	  ncy = nyk(k)
	  ipc = kpbgn(k)
	  irc = krbgn(k)
c
c     transfer down to all grid levels
c
	  call trsfc2(nx,ny,wk(ip),wk(ir),ncx,ncy,
     +                wk(ipc),wk(irc))
	end do
c
c     adjust right hand side at all grid levels in case
c     rhs or specified b.c. in phi or gbdy changed
c
	do k=1,ngrid
	  nx = nxk(k)
	  ny = nyk(k)
	  ip = kpbgn(k)
	  ir = krbgn(k)
	  call adjmd2sp(nx,ny,wk(ip),wk(ir),bndyc,cfx,cfy)
	end do
c
c     execute one full multigrid cycle
c
	do k=1,ngrid-1
	  kcur = k
	  call kcymd2sp(wk)
	  nx = nxk(k+1)
	  ny = nyk(k+1)
	  ip = kpbgn(k+1)
	  ipc = kpbgn(k)
	  ncx = nxk(k)
	  ncy = nyk(k)
c
c     lift or prolong approximation from k to k+1
c
	  call prolon2(ncx,ncy,wk(ipc),nx,ny,wk(ip),nxa,nxb,
     +                 nyc,nyd,intpol)
	end do
      else
c
c     adjust rhs at finest grid level only
c
	nx = nxk(ngrid)
	ny = nyk(ngrid)
	ip = kpbgn(ngrid)
	ir = krbgn(ngrid)
	call adjmd2sp(nx,ny,wk(ip),wk(ir),bndyc,cfx,cfy)
      end if
c
c     execute maxcy more multigrid k cycles from finest level
c
      kcur = ngrid
      do iter=1,maxcy
	itero = iter
	call kcymd2sp(wk)
	if (tolmax.gt.0.0) then
c
c      error control
c
	  relmax = 0.0
	  phmax = 0.0
	  do j=1,nfy
	    jj = j*(nfx+2)
	    do i=1,nfx
	      ij = jj+i+1
	      phmax = amax1(phmax,abs(wk(ij)))
	      relmax = amax1(relmax,abs(wk(ij)-phif(i,j)))
	      phif(i,j) = wk(ij)
	    end do
	  end do
c
c     set maximum relative difference and check for convergence
c
	  if (phmax.gt.0.0) relmax = relmax/phmax
	  if (relmax.le.tolmax) return
	end if
      end do
c
c     set final interate after maxcy cycles in phif
c
      do j=1,nfy
	jj = j*(nfx+2)
	do i=1,nfx
	  ij = jj+i+1
	  phif(i,j) = wk(ij)
	end do
      end do
      return
      end

      subroutine kcymd2sp(wk)
c
c     execute multigrid k cycle from kcur grid level
c     kcycle=1 for v cycles, kcycle=2 for w cycles
c
      implicit none
      real wk(*)
      integer intl,nxa,nxb,nyc,nyd,ixp,jyq,iex,jey,nfx,nfy,iguess,
     +             maxcy,method,nwork,lwork,itero,ngrid,klevel,kcur,
     +             kcycle,iprer,ipost,intpol,kps
      integer nx,ny,ip,ipc,ir,irc,icx,icy,itx,ity,ncx,ncy,l,nrel
      real xa,xb,yc,yd,tolmax,relmax
      integer kpbgn,krbgn,kcxbgn,kcybgn,ktxbgn,ktybgn,nxk,nyk,isx,jsy
      common/imud2/intl,nxa,nxb,nyc,nyd,ixp,jyq,iex,jey,nfx,nfy,iguess,
     +             maxcy,method,nwork,lwork,itero,ngrid,klevel,kcur,
     +             kcycle,iprer,ipost,intpol,kps
      common/fmud2/xa,xb,yc,yd,tolmax,relmax
      common/mud2spc/kpbgn(50),krbgn(50),kcxbgn(50),kcybgn(50),
     +ktxbgn(50),ktybgn(50),nxk(50),nyk(50),isx,jsy
      integer kount(50)
      klevel = kcur
      nx = nxk(klevel)
      ny = nyk(klevel)
      ip = kpbgn(klevel)
      ir = krbgn(klevel)
      icx = kcxbgn(klevel)
      icy = kcybgn(klevel)
      itx = ktxbgn(klevel)
      ity = ktybgn(klevel)
c
c     prerelax at current finest grid level
c
      do l=1,iprer
	call relmd2sp(nx,ny,wk(ip),wk(ir),wk(icx),wk(icy),wk(itx),
     +                wk(ity),wk(kps))
      end do
      if (kcur .eq. 1) go to 5
c
c     restrict residual to kcur-1 level
c
      ipc = kpbgn(klevel-1)
      ncx = nxk(klevel-1)
      ncy = nyk(klevel-1)
      irc = krbgn(klevel-1)
      call resmd2sp(nx,ny,wk(ip),wk(ir),ncx,ncy,wk(ipc),wk(irc),
     +              wk(icx),wk(icy),wk(kps))
c
c    set counter for grid levels to zero
c
      do l = 1,kcur
	kount(l) = 0
      end do
c
c    set new grid level and continue k-cycling
c
      klevel = kcur-1
      nrel = iprer
c
c   kcycle control point
c
   10 continue
c
c      post relax when kcur revisited
c
      if (klevel .eq. kcur) go to 5
c
c   count hit at current level
c
      kount(klevel) = kount(klevel)+1
c
c   relax at current level
c
      nx = nxk(klevel)
      ny = nyk(klevel)
      ip = kpbgn(klevel)
      ir = krbgn(klevel)
      icx = kcxbgn(klevel)
      icy = kcybgn(klevel)
      itx = ktxbgn(klevel)
      ity = ktybgn(klevel)
      do l=1,nrel
	call relmd2sp(nx,ny,wk(ip),wk(ir),wk(icx),wk(icy),wk(itx),
     +                wk(ity),wk(kps))
      end do
      if (kount(klevel) .eq. kcycle+1) then
c
c     kcycle complete at klevel
c
	ipc = ip
	ncx = nx
	ncy = ny
	ip = kpbgn(klevel+1)
	ir = krbgn(klevel+1)
	icx = kcxbgn(klevel+1)
	icy = kcybgn(klevel+1)
	nx = nxk(klevel+1)
	ny = nyk(klevel+1)
c
c    inject correction to finer grid
c
	call cor2(nx,ny,wk(ip),ncx,ncy,wk(ipc),nxa,nxb,nyc,nyd,
     +            intpol,wk(kps))
c
c    reset counter to zero
c
	kount(klevel) = 0
c
c     ascend to next higher level and set to postrelax there
c
	klevel = klevel+1
	nrel = ipost
	go to 10
      else
	if (klevel .gt. 1) then
c
c    kcycle not complete so descend unless at coarsest grid
c
	  ipc = kpbgn(klevel-1)
	  ncx = nxk(klevel-1)
	  ncy = nyk(klevel-1)
	  irc = krbgn(klevel-1)
	  call resmd2sp(nx,ny,wk(ip),wk(ir),ncx,ncy,wk(ipc),wk(irc),
     +                  wk(icx),wk(icy),wk(kps))
c
c     prerelax at next coarser level
c
	  klevel = klevel-1
	  nrel = iprer
	  go to 10
	else
c
c    postrelax at coarsest level
c
	  do l=1,ipost
	    call relmd2sp(nx,ny,wk(ip),wk(ir),wk(icx),wk(icy),wk(itx),
     +                    wk(ity),wk(kps))
	  end do
	  ipc = ip
	  ncx = nx
	  ncy = ny
	  ip = kpbgn(2)
	  ir = krbgn(2)
	  icx = kcxbgn(2)
	  icy = kcybgn(2)
	  nx = nxk(2)
	  ny = nyk(2)
c
c    inject correction to level 2
c
	  call cor2(nx,ny,wk(ip),ncx,ncy,wk(ipc),nxa,nxb,nyc,nyd,
     +              intpol,wk(kps))
c
c     set to postrelax at level 2
c
	  nrel = ipost
	  klevel = 2
	  go to 10
	end if
      end if
    5 continue
c
c     post relax at current finest grid level
c
      nx = nxk(kcur)
      ny = nyk(kcur)
      ip = kpbgn(kcur)
      ir = krbgn(kcur)
      icx = kcxbgn(kcur)
      icy = kcybgn(kcur)
      itx = ktxbgn(kcur)
      ity = ktybgn(kcur)
      do l=1,ipost
	call relmd2sp(nx,ny,wk(ip),wk(ir),wk(icx),wk(icy),wk(itx),
     +                wk(ity),wk(kps))
      end do
      return
      end

      subroutine dismd2sp(nx,ny,cofx,cofy,tx,ty,bndyc,cfx,cfy,wk,ier)
c
c     discretize elliptic pde for mud2, set nonfatal errors
c
      implicit none
      integer intl,nxa,nxb,nyc,nyd,ixp,jyq,iex,jey,nfx,nfy,iguess,
     +             maxcy,method,nwork,lwork,itero,ngrid,klevel,kcur,
     +             kcycle,iprer,ipost,intpol,kps
      real xa,xb,yc,yd,tolmax,relmax
      integer nx,ny,ist,ifn,jst,jfn,i,j,kbdy,im1,jm1,ier
      real cofx(nx,3),cofy(ny,3),tx(nx,ny,*),ty(ny,nx,*)
      real wk(*),dlx,dlx2,dlxx,dly,dly2,dlyy,cmin,alfmax,cemax
      real x,y,cxx,cyy,cx,cy,cex,cey,c1,c2,c3,alfa,gbdy
      common/imud2/intl,nxa,nxb,nyc,nyd,ixp,jyq,iex,jey,nfx,nfy,iguess,
     +             maxcy,method,nwork,lwork,itero,ngrid,klevel,kcur,
     +             kcycle,iprer,ipost,intpol,kps
      common/fmud2/xa,xb,yc,yd,tolmax,relmax
      external bndyc,cfx,cfy
      dlx = (xb-xa)/(nx-1)
      dlx2 = dlx+dlx
      dlxx = dlx*dlx
      dly = (yd-yc)/(ny-1)
      dly2 = dly+dly
      dlyy = dly*dly
      cmin = 1.0
      alfmax = 0.0
      cemax = 0.0
c
c     set x,y subscript limits for calls to cofx,cofy
c     (avoid specified boundaries)
c
      ist = 1
      ifn = nx
      jst = 1
      jfn = ny
      if (nxa.eq.1) ist = 2
      if (nxb.eq.1) ifn = nx-1
      if (nyc.eq.1) jst = 2
      if (nyd.eq.1) jfn = ny-1
c
c     compute discretization coefficients on interior and
c     nonspecified boundaries
c
      do j=jst,jfn
	y = yc+(j-1)*dly
c	call cfy(y,cyy,cy,cey)
	cyy = 1.
	cy = 0.
	cey = 0.
	cmin = amin1(cmin,cyy)
	cemax = amax1(abs(cey),cemax)
c
c     flag hyperbolic pde if necessary
c
	if (klevel.eq.ngrid) then
	  if (abs(cy)*dly.gt.2*abs(cyy)) then
	    ier = -4
	  end if
	end if
	cyy = amax1(cyy,abs(cy)*dly*0.5)
	c1 = cyy/dlyy-cy/dly2
	c2 = cyy/dlyy+cy/dly2
	c3 = cey-(c1+c2)
	cofy(j,1) = c1
	cofy(j,2) = c2
	cofy(j,3) = c3
      end do
      do i=ist,ifn
	x = xa+(i-1)*dlx
c	call cfx(x,cxx,cx,cex)
	cxx = 1.
	cx = 0.
	cex = 0.
	cmin = amin1(cmin,cxx)
	cemax = amax1(abs(cex),cemax)
c
c      flag hyperbolic pde if necessary
c
	if (klevel.eq.ngrid) then
	  if (abs(cx)*dlx.gt.2*abs(cxx)) then
	    ier = -4
	  end if
	end if
	cxx = amax1(cxx,abs(cx)*dlx*0.5)
	c1 = cxx/dlxx-cx/dlx2
	c2 = cxx/dlxx+cx/dlx2
	c3 = cex-(c1+c2)
	cofx(i,1) = c1
	cofx(i,2) = c2
	cofx(i,3) = c3
      end do
c
c     adjust discretization for mixed derivative b.c.
c
      if (nxa.eq.2) then
	kbdy = 1
	i = 1
	c1 = cofx(i,1)
	cofx(i,1) = 0.0
	cofx(i,2) = cofx(i,2)+c1
	y = yc+dly
c
c     compute constant coefficient alfa
c
	call bndyc(kbdy,y,alfa,gbdy)
	alfmax = amax1(alfmax,abs(alfa))
	cofx(i,3) = cofx(i,3)+dlx2*alfa*c1
      end if
      if (nxb.eq.2) then
	kbdy = 2
	i = nx
	y = yc+dly
c
c     compute constant coefficient alfa
c
	call bndyc(kbdy,y,alfa,gbdy)
	c2 = cofx(i,2)
	cofx(i,1) = cofx(i,1)+c2
	cofx(i,2) = 0.0
	cofx(i,3) = cofx(i,3)-dlx2*alfa*c2
	alfmax = amax1(abs(alfa),alfmax)
      end if
      if (nyc.eq.2) then
	kbdy = 3
	j = 1
	x = xa+dlx
c
c     compute constant coefficient alfa
c
	call bndyc(kbdy,x,alfa,gbdy)
	c1 = cofy(j,1)
	cofy(j,1) = 0.0
	cofy(j,2) = cofy(j,2) + c1
	cofy(j,3) = cofy(j,3) + dly2*alfa*c1
	alfmax = amax1(abs(alfa),alfmax)
      end if
      if (nyd.eq.2) then
	kbdy = 4
	j = ny
	x = xa+dlx
c
c     compute constant coefficient alfa
c
	call bndyc(kbdy,x,alfa,gbdy)
	c2 = cofy(j,2)
	cofy(j,2) = 0.0
	cofy(j,1) = cofy(j,1) + c2
	cofy(j,3) = cofy(j,3) - dly2*c2*alfa
	alfmax = amax1(abs(alfa),alfmax)
      end if
c
c     if detected then flag singular pde
c
	if (cemax.eq.0.0.and.alfmax.eq.0.0) then
	  if (nxa.eq.0.or.(nxa.eq.2.and.nxb.eq.2)) then
	    if (nyc.eq.0.or.(nyc.eq.2.and.nyd.eq.2)) then
	      ier = -3
	   end if
	  end if
	end if
c
c     if detected then flag nonellipticity
c
      if (cmin.le.0.0) then
	ier = -2
      end if
c
c     set and factor tridiagonal matrices for line relaxation(s) if flagged
c
      if (method.eq.1.or.method.eq.3) then
	if (nxa.ne.0) then
c
c    nonperiodic x line relaxation
c
	  do i=1,nx
	    im1 = max0(i-1,1)
	    do j=1,ny
	      tx(im1,j,1) = cofx(i,1)
	      tx(i,j,2) = cofx(i,3)+cofy(j,3)
	      tx(i,j,3) = cofx(i,2)
	    end do
	  end do
	  if (nxa.eq.1) then
	    do j=1,ny
	      tx(1,j,2) = 1.0
	      tx(1,j,3) = 0.0
	    end do
	  end if
	  if (nxb.eq.1) then
	    do j=1,ny
	      tx(nx-1,j,1) = 0.0
	      tx(nx,j,2) = 1.0
	    end do
	  end if
	  call factri(ny,nx,tx(1,1,1),tx(1,1,2),tx(1,1,3))
	else
c
c     periodic x line relaxation
c
	  if (nx .gt. 3) then
c
c     set and factor iff nx > 3
c
	    do i=1,nx-1
	      do j=1,ny
		tx(i,j,1) = cofx(i,1)
		tx(i,j,2) = cofx(i,3)+cofy(j,3)
		tx(i,j,3) = cofx(i,2)
	      end do
	    end do
	    call factrp(ny,nx,tx,tx(1,1,2),tx(1,1,3),tx(1,1,4),
     +                  tx(1,1,5),wk(kps))
	  end if
	end if
      end if

      if (method.eq.2.or.method.eq.3) then
	if (nyc.ne.0) then
c
c     nonperiodic y line relaxation
c
	  do j=1,ny
	    jm1 = max0(j-1,1)
	    do i=1,nx
	      ty(jm1,i,1) = cofy(j,1)
	      ty(j,i,2) = cofy(j,3)+cofx(i,3)
	      ty(j,i,3) = cofy(j,2)
	    end do
	  end do
	  if (nyc.eq.1) then
	    do i=1,nx
	      ty(1,i,2) = 1.0
	      ty(1,i,3) = 0.0
	    end do
	  end if
	  if (nyd.eq.1) then
	    do i=1,nx
	      ty(ny-1,i,1) = 0.0
	      ty(ny,i,2) = 1.0
	    end do
	  end if
	  call factri(nx,ny,ty(1,1,1),ty(1,1,2),ty(1,1,3))
	else
c
c      periodic y line relaxation
c
	  if (ny .gt. 3) then
c
c     set and factor iff ny > 3
c
	    do j=1,ny-1
	      do i=1,nx
		ty(j,i,1) = cofy(j,1)
		ty(j,i,2) = cofy(j,3)+cofx(i,3)
		ty(j,i,3) = cofy(j,2)
	      end do
	    end do
	    call factrp(nx,ny,ty,ty(1,1,2),ty(1,1,3),ty(1,1,4),
     +                  ty(1,1,5),wk(kps))
	  end if
	end if
      end if
      return
      end

      subroutine adjmd2sp(nx,ny,phi,rhs,bndyc,cfx,cfy)
c
c     adjust righthand side for various boundary conditions
c
      implicit none
      integer intl,nxa,nxb,nyc,nyd,ixp,jyq,iex,jey,nfx,nfy,iguess,
     +             maxcy,method,nwork,lwork,itero,ngrid,klevel,kcur,
     +             kcycle,iprer,ipost,intpol,kps
      real xa,xb,yc,yd,tolmax,relmax
      integer nx,ny,ist,ifn,jst,jfn,i,j,kbdy
      real phi(0:nx+1,0:ny+1),rhs(nx,ny)
      real dlx,dlx2,dlxx,dly,dly2,dlyy
      real x,y,cxx,cyy,cx,cy,cex,cey,c1,c2,alfa,gbdy
      common/imud2/intl,nxa,nxb,nyc,nyd,ixp,jyq,iex,jey,nfx,nfy,iguess,
     +             maxcy,method,nwork,lwork,itero,ngrid,klevel,kcur,
     +             kcycle,iprer,ipost,intpol,kps
      common/fmud2/xa,xb,yc,yd,tolmax,relmax
      external bndyc,cfx,cfy
      dlx = (xb-xa)/(nx-1)
      dly = (yd-yc)/(ny-1)
      dlx2 = dlx+dlx
      dly2 = dly+dly
      dlxx = dlx*dlx
      dlyy = dly*dly
      ist = 1
      ifn = nx
      jst = 1
      jfn = ny
      if (nxa.eq.1) ist = 2
      if (nxb.eq.1) ifn = nx-1
      if (nyc.eq.1) jst = 2
      if (nyd.eq.1) jfn = ny-1
c
c     adjust right hand side at derivative boundaries
c
      if (nxa.eq.2) then
	kbdy = 1
	x = xa
	i = 1
c	call cfx(x,cxx,cx,cex)
	cxx = 1.
	cx = 0.
	cex = 0.
	cxx = amax1(cxx,abs(cx)*dlx*0.5)
	c1 = cxx/dlxx-cx/dlx2
	do j=jst,jfn
	  y = yc+(j-1)*dly
	  call bndyc(kbdy,y,alfa,gbdy)
	  rhs(i,j) = rhs(i,j)+dlx2*c1*gbdy
	end do
      end if
      if (nxb.eq.2) then
	kbdy = 2
	x = xb
	i = nx
c	call cfx(x,cxx,cx,cex)
	cxx = 1.
	cx = 0.
	cex = 0.
	cxx = amax1(cxx,abs(cx)*dlx*0.5)
	c2 = cxx/dlxx+cx/dlx2
	do j=jst,jfn
	  y = yc+(j-1)*dly
	  call bndyc(kbdy,y,alfa,gbdy)
	  rhs(i,j) = rhs(i,j)-dlx2*c2*gbdy
	end do
      end if
      if (nyc.eq.2) then
	kbdy = 3
	y = yc
	j = 1
c	call cfy(y,cyy,cy,cey)
	cyy = 1.
	cy = 0.
	cey = 0.
	cyy = amax1(cyy,abs(cy)*dly*0.5)
	c1 = cyy/dlyy-cy/dly2
	do i=ist,ifn
	  x = xa+(i-1)*dlx
	  call bndyc(kbdy,x,alfa,gbdy)
	  rhs(i,j) = rhs(i,j)+dly2*c1*gbdy
	end do
      end if
      if (nyd.eq.2) then
	kbdy = 4
	y = yd
	j = ny
c	call cfy(y,cyy,cy,cey)
	cyy = 1.
	cy = 0.
	cey = 0.
	cyy = amax1(cyy,abs(cy)*dly*0.5)
	c2 = cyy/dlyy+cy/dly2
	do i=ist,ifn
	  x = xa+(i-1)*dlx
	  call bndyc(kbdy,x,alfa,gbdy)
	  rhs(i,j) = rhs(i,j)-dly2*c2*gbdy
	end do
      end if
c
c     set specified boundaries in rhs from phi
c
      if (nxa.eq.1) then
	i = 1
	do j=1,ny
	  rhs(i,j) = phi(i,j)
	end do
      end if
      if (nxb.eq.1) then
	i = nx
	do j=1,ny
	  rhs(i,j) = phi(i,j)
	end do
      end if
      if (nyc.eq.1) then
	j = 1
	do i=1,nx
	  rhs(i,j) = phi(i,j)
	end do
      end if
      if (nyd.eq.1) then
	j = ny
	do i=1,nx
	  rhs(i,j) = phi(i,j)
	end do
      end if
      return
      end

      subroutine resmd2sp(nx,ny,phi,rhs,ncx,ncy,phic,rhsc,cofx,cofy,
     +                    resf)
c
c     restrict residual from fine to coarse mesh using fully weighted
c     residual restriction
c
      implicit none
      integer intl,nxa,nxb,nyc,nyd,ixp,jyq,iex,jey,nfx,nfy,iguess,
     +             maxcy,method,nwork,lwork,itero,ngrid,klevel,kcur,
     +             kcycle,iprer,ipost,intpol,kps
      integer nx,ny,ncx,ncy,i,j,ic,jc,ist,ifn,jst,jfn
      common/imud2/intl,nxa,nxb,nyc,nyd,ixp,jyq,iex,jey,nfx,nfy,iguess,
     +             maxcy,method,nwork,lwork,itero,ngrid,klevel,kcur,
     +             kcycle,iprer,ipost,intpol,kps
      real rhs(nx,ny),rhsc(ncx,ncy),resf(nx,ny)
      real phi(0:nx+1,0:ny+1),phic(0:ncx+1,0:ncy+1)
      real cofx(nx,3),cofy(ny,3)
c
c     set phic zero
c
      do jc=0,ncy+1
	do ic=0,ncx+1
	  phic(ic,jc) = 0.0
	end do
      end do
c
c     intialize residual to zero and set limits
c
      do j=1,ny
	do i=1,nx
	  resf(i,j) = 0.0
	end do
      end do
      ist = 1
      if (nxa.eq.1) ist = 2
      ifn = nx
      if (nxb.eq.1) ifn = nx-1
      jst = 1
      if (nyc.eq.1) jst = 2
      jfn = ny
      if (nyd.eq.1) jfn = ny-1
c
c     compute residual on fine mesh in resf
c
!$OMP PARALLEL DO SHARED(jst,jfn,ist,ifn,rhs,resf,cofx,cofy,phi)
!$OMP+PRIVATE(i,j)
      do j=jst,jfn
	do i=ist,ifn
	  resf(i,j) =  rhs(i,j)-(
     +                 cofx(i,1)*phi(i-1,j)+
     +                 cofx(i,2)*phi(i+1,j)+
     +                 cofy(j,1)*phi(i,j-1)+
     +                 cofy(j,2)*phi(i,j+1)+
     +                 (cofx(i,3)+cofy(j,3))*phi(i,j))
	end do
      end do
c
c     restrict resf to coarse mesh in rhsc
c
      call res2(nx,ny,resf,ncx,ncy,rhsc,nxa,nxb,nyc,nyd)
      return
      end

      subroutine relmd2sp(nx,ny,phi,rhs,cofx,cofy,tx,ty,sum)
c
c     relaxation for mud2sp
c
      implicit none
      integer nx,ny
      real phi(*),rhs(*),cofx(*),cofy(*),tx(*),ty(*),sum(*)
      integer intl,nxa,nxb,nyc,nyd,ixp,jyq,iex,jey,nfx,nfy,iguess,
     +             maxcy,method,nwork,lwork,itero,ngrid,klevel,kcur,
     +             kcycle,iprer,ipost,intpol,kps
      common/imud2/intl,nxa,nxb,nyc,nyd,ixp,jyq,iex,jey,nfx,nfy,iguess,
     +             maxcy,method,nwork,lwork,itero,ngrid,klevel,kcur,
     +             kcycle,iprer,ipost,intpol,kps
      if (method.eq.0) then                ! point relaxation
	call rmd2spp(nx,ny,phi,rhs,cofx,cofy)
      else if (method.eq.1) then           ! line x relaxation
	call slxmd2sp(nx,ny,phi,rhs,cofx,cofy,tx,sum)
      else if (method.eq.2) then           ! line y relaxation
	call slymd2sp(nx,ny,phi,rhs,cofx,cofy,ty,sum)
      else if (method.eq.3) then           ! line x&y relaxation
	call slxmd2sp(nx,ny,phi,rhs,cofx,cofy,tx,sum)
	call slymd2sp(nx,ny,phi,rhs,cofx,cofy,ty,sum)
      end if
      return
      end

      subroutine rmd2spp(nx,ny,phi,rhs,cofx,cofy)
c
c     gauss-seidel red/black point relaxation
c
      implicit none
      integer nx,ny,i,j,ist,ifn,jst,jfn
      integer intl,nxa,nxb,nyc,nyd,ixp,jyq,iex,jey,nfx,nfy,iguess,
     +             maxcy,method,nwork,lwork,itero,ngrid,klevel,kcur,
     +             kcycle,iprer,ipost,intpol,kps
      common/imud2/intl,nxa,nxb,nyc,nyd,ixp,jyq,iex,jey,nfx,nfy,iguess,
     +             maxcy,method,nwork,lwork,itero,ngrid,klevel,kcur,
     +             kcycle,iprer,ipost,intpol,kps
      real phi(0:nx+1,0:ny+1),rhs(nx,ny),cofx(nx,3),cofy(ny,3)
c
c     set loop limits to avoid specified boundaries
c     in red/black sweeps
c
      ist = 1
      if (nxa.eq.1) ist = 3
      ifn = nx
      if (nxb.eq.1) ifn = nx-1
      jst = 1
      if (nyc.eq.1) jst = 3
      jfn = ny
      if (nyd.eq.1) jfn = ny-1
c
c    periodic adjustment bypass block
c
      if (nxa*nyc.ne.0) then
c
c     relax on red grid points
c
!$OMP PARALLEL DO SHARED(rhs,cofx,cofy,phi,ist,ifn,jst,jfn) PRIVATE(i,j)
	do i=ist,ifn,2
	  do j=jst,jfn,2
	    phi(i,j) = (rhs(i,j) - (
     +                  cofx(i,1)*phi(i-1,j)+cofx(i,2)*phi(i+1,j) +
     +                  cofy(j,1)*phi(i,j-1)+cofy(j,2)*phi(i,j+1)))/
     +                 (cofx(i,3)+cofy(j,3))
	  end do
	end do
!$OMP PARALLEL DO SHARED(rhs,cofx,cofy,phi,ifn,jfn) PRIVATE(i,j)
	do i=2,ifn,2
	  do j=2,jfn,2
	    phi(i,j) = (rhs(i,j) -
     +                 (cofx(i,1)*phi(i-1,j)+cofx(i,2)*phi(i+1,j) +
     +                 cofy(j,1)*phi(i,j-1)+cofy(j,2)*phi(i,j+1)))/
     +                 (cofx(i,3)+cofy(j,3))
	  end do
	end do
c
c     relax on black grid points
c
!$OMP PARALLEL DO SHARED(rhs,cofx,cofy,phi,ist,ifn,jfn) PRIVATE(i,j)
	do i=ist,ifn,2
	  do j=2,jfn,2
	    phi(i,j) = (rhs(i,j) -
     +                 (cofx(i,1)*phi(i-1,j)+cofx(i,2)*phi(i+1,j) +
     +                 cofy(j,1)*phi(i,j-1)+cofy(j,2)*phi(i,j+1)))/
     +                 (cofx(i,3)+cofy(j,3))
	  end do
	end do
!$OMP PARALLEL DO SHARED(rhs,cofx,cofy,phi,ifn,jst,jfn) PRIVATE(i,j)
	do i=2,ifn,2
	  do j=jst,jfn,2
	    phi(i,j) = (rhs(i,j) -
     +                 (cofx(i,1)*phi(i-1,j)+cofx(i,2)*phi(i+1,j) +
     +                 cofy(j,1)*phi(i,j-1)+cofy(j,2)*phi(i,j+1)))/
     +                 (cofx(i,3)+cofy(j,3))
	  end do
	end do
	return
      end if
c
c    set periodic virtual boundaries
c
      if (nxa.eq.0) then
	do j=1,ny
	  phi(0,j) = phi(nx-1,j)
	  phi(nx+1,j) = phi(2,j)
	end do
      end if
      if (nyc.eq.0) then
	do i=1,nx
	  phi(i,0) = phi(i,ny-1)
	  phi(i,ny+1) = phi(i,2)
	end do
      end if
c
c     relax on red grid points
c
!$OMP PARALLEL DO SHARED(rhs,cofx,cofy,phi,ist,ifn,jst,jfn) PRIVATE(i,j)
      do i=ist,ifn,2
	do j=jst,jfn,2
	    phi(i,j) = (rhs(i,j) -
     +                 (cofx(i,1)*phi(i-1,j)+cofx(i,2)*phi(i+1,j) +
     +                 cofy(j,1)*phi(i,j-1)+cofy(j,2)*phi(i,j+1)))/
     +                 (cofx(i,3)+cofy(j,3))
	end do
      end do
!$OMP PARALLEL DO SHARED(rhs,cofx,cofy,phi,ifn,jfn) PRIVATE(i,j)
      do i=2,ifn,2
	do j=2,jfn,2
	    phi(i,j) = (rhs(i,j) -
     +                 (cofx(i,1)*phi(i-1,j)+cofx(i,2)*phi(i+1,j) +
     +                 cofy(j,1)*phi(i,j-1)+cofy(j,2)*phi(i,j+1)))/
     +                 (cofx(i,3)+cofy(j,3))
	end do
      end do
c
c    ensure periodic virtual boundary red points are set
c
      if (nxa.eq.0) then
	do j=1,ny
	  phi(0,j) = phi(nx-1,j)
	  phi(nx+1,j) = phi(2,j)
	end do
      end if
      if (nyc.eq.0) then
	do i=1,nx
	  phi(i,0) = phi(i,ny-1)
	  phi(i,ny+1) = phi(i,2)
	end do
      end if
c
c     relax on black grid points
c
c
!$OMP PARALLEL DO SHARED(rhs,cofx,cofy,phi,ist,ifn,jfn) PRIVATE(i,j)
      do i=ist,ifn,2
	do j=2,jfn,2
	    phi(i,j) = (rhs(i,j) -
     +                 (cofx(i,1)*phi(i-1,j)+cofx(i,2)*phi(i+1,j) +
     +                 cofy(j,1)*phi(i,j-1)+cofy(j,2)*phi(i,j+1)))/
     +                 (cofx(i,3)+cofy(j,3))
	end do
      end do
!$OMP PARALLEL DO SHARED(rhs,cofx,cofy,phi,ifn,jst,jfn) PRIVATE(i,j)
      do i=2,ifn,2
	do j=jst,jfn,2
	    phi(i,j) = (rhs(i,j) -
     +                 (cofx(i,1)*phi(i-1,j)+cofx(i,2)*phi(i+1,j) +
     +                 cofy(j,1)*phi(i,j-1)+cofy(j,2)*phi(i,j+1)))/
     +                 (cofx(i,3)+cofy(j,3))
	end do
      end do
c
c     final set of periodic virtual boundaries
c
      if (nxa.eq.0) then
	do j=1,ny
	  phi(0,j) = phi(nx-1,j)
	  phi(nx+1,j) = phi(2,j)
	end do
      end if
      if (nyc.eq.0) then
	do i=1,nx
	  phi(i,0) = phi(i,ny-1)
	  phi(i,ny+1) = phi(i,2)
	end do
      end if
      return
      end

      subroutine slxmd2sp(nx,ny,phi,rhs,cofx,cofy,tx,sum)
c
c     line relaxation in the x direction (periodic or nonperiodic)
c
      implicit none
      integer nx,ny,i,ib,j,ist,ifn,jst,jfn,jm1,jp1
      integer intl,nxa,nxb,nyc,nyd,ixp,jyq,iex,jey,nfx,nfy,iguess,
     +             maxcy,method,nwork,lwork,itero,ngrid,klevel,kcur,
     +             kcycle,iprer,ipost,intpol,kps
      common/imud2/intl,nxa,nxb,nyc,nyd,ixp,jyq,iex,jey,nfx,nfy,iguess,
     +             maxcy,method,nwork,lwork,itero,ngrid,klevel,kcur,
     +             kcycle,iprer,ipost,intpol,kps
      real phi(0:nx+1,0:ny+1),cofx(nx,3),cofy(ny,3),tx(nx,ny,*),sum(ny)
      real rhs(nx,ny)
c
c     replace line x with point gauss-seidel if
c     x direction is periodic and nx = 3 (coarsest)
c
      if (nxa .eq. 0 .and. nx .eq. 3) then
	call rmd2spp(nx,ny,phi,rhs,cofx,cofy)
	return
      end if
      jst = 3
      jfn = ny-1
      ist = 2
      if (nxa.ne.1) ist = 1
      ifn = nx-1
      if (nxb.ne.1) ifn = nx
      if (nxa.ne.0) then
c
c     non-periodic line relaxation in x direction
c     lines thru odd y points
c
	if (nyc.ne.1) then
	  jst = 1
	  j = 1
	  jm1 = ny-1
	  jp1 = 2
	  do i=ist,ifn
	    phi(i,j)=rhs(i,j)-(cofy(j,1)*phi(i,jm1)+cofy(j,2)*phi(i,jp1))
	  end do
	end if
!$OMP PARALLEL DO SHARED(ist,ifn,phi,rhs,cofy,ny) PRIVATE(i,j)
	do j=3,ny-1,2
	  do i=ist,ifn
	    phi(i,j)=rhs(i,j)-(cofy(j,1)*phi(i,j-1)+cofy(j,2)*phi(i,j+1))
	  end do
	end do
	if (mod(ny,2).ne.0.and.nyd.ne.1) then
	  jfn = ny
	  j = ny
	  jm1 = ny-1
	  jp1 = 2
	  do i=ist,ifn
	    phi(i,j)=rhs(i,j)-(cofy(j,1)*phi(i,jm1)+cofy(j,2)*phi(i,jp1))
	  end do
	end if
c
c     forward sweep
c
!$OMP PARALLEL DO SHARED(jst,jfn,phi,tx,nx) PRIVATE(i,ib,j)
      do j=jst,jfn,2
	do i=2,nx
	  phi(i,j) = phi(i,j)-tx(i-1,j,1)*phi(i-1,j)
	end do
c
c     backward sweep
c
	phi(nx,j) = phi(nx,j)/tx(nx,j,2)
	do ib=2,nx
	  i = nx-ib+1
	    phi(i,j) = (phi(i,j)-tx(i,j,3)*phi(i+1,j))/tx(i,j,2)
	end do
      end do
c
c     lines thru even y points
c
	jfn = ny-1
!$OMP PARALLEL DO SHARED(ist,ifn,cofy,rhs,phi,ny) PRIVATE(i,j)
	do j=2,ny-1,2
	  do i=ist,ifn
	    phi(i,j)=rhs(i,j)-(cofy(j,1)*phi(i,j-1)+cofy(j,2)*phi(i,j+1))
	  end do
	end do
	if (mod(ny,2).eq.0.and.nyd.ne.1) then
	  jfn = ny
	  j = ny
	  jm1 = ny-1
	  jp1 = 2
	  do i=ist,ifn
	    phi(i,j)=rhs(i,j)-cofy(j,1)*phi(i,jm1)-cofy(j,2)*phi(i,jp1)
	  end do
	end if
c
c     forward sweep
c
!$OMP PARALLEL DO SHARED(jfn,phi,tx,nx) PRIVATE(i,ib,j)
	do j=2,jfn,2
	  do i=2,nx
	    phi(i,j) = phi(i,j)-tx(i-1,j,1)*phi(i-1,j)
	  end do
c
c     backward sweep
c
	  phi(nx,j) = phi(nx,j)/tx(nx,j,2)
	  do ib=2,nx
	    i = nx-ib+1
	    phi(i,j) = (phi(i,j)-tx(i,j,3)*phi(i+1,j))/tx(i,j,2)
	  end do
	end do
c
c     adjust for periodic y if ny even (only possible at coarsest level)
c     red/black line will not allow phi(i,1) = phi(i,ny)
c
	if (nyc.eq.0.and.mod(ny,2).eq.0) then
	  do i=1,nx
	    phi(i,ny) = phi(i,1)
	  end do
	end if
      else
c
c     line periodic relaxation in x direction
c
	do j=1,ny
	  sum(j) = 0.0
	end do
c
c     set rhs and solve on lines thru odd y points
c
	jst = 3
	jfn = ny-1
	if (nyc.ne.1) then
	  jst = 1
	  j = 1
	  jm1 = ny-1
	  jp1 = 2
	  do i=1,nx-1
	    phi(i,j)=rhs(i,j)-cofy(j,1)*phi(i,jm1)-cofy(j,2)*phi(i,jp1)
	  end do
	end if
!$OMP PARALLEL DO SHARED(cofy,rhs,phi,nx,ny) PRIVATE(i,j)
	do j=3,ny-1,2
	  do i=1,nx-1
	    phi(i,j)=rhs(i,j)-cofy(j,1)*phi(i,j-1)-cofy(j,2)*phi(i,j+1)
	  end do
	end do
c
c     set last y point line if odd and non-specified
c
	if (mod(ny,2).ne.0.and.nyd.ne.1) then
	  jfn = ny
	  j = ny
	  jm1 = ny-1
	  jp1 = 2
	  do i=1,nx-1
	    phi(i,j)=rhs(i,j)-cofy(j,1)*phi(i,jm1)-cofy(j,2)*phi(i,jp1)
	  end do
	end if
c
c     forward sweep
c
!$OMP PARALLEL DO SHARED(sum,jst,jfn,phi,tx,nx) PRIVATE(i,ib,j)
	do j=jst,jfn,2
	  do i=2,nx-2
	    phi(i,j) = phi(i,j)-tx(i,j,1)*phi(i-1,j)
	  end do
	  do i=1,nx-2
	    sum(j) = sum(j)+tx(i,j,5)*phi(i,j)
	  end do
	  phi(nx-1,j) = phi(nx-1,j)-sum(j)
c
c     backward sweep
c
	  phi(nx-1,j) = phi(nx-1,j)/tx(nx-1,j,2)
	  phi(nx-2,j) = (phi(nx-2,j)-tx(nx-2,j,4)*phi(nx-1,j))
     +                   /tx(nx-2,j,2)
	  do ib=4,nx
	    i = nx-ib+1
	    phi(i,j) = (phi(i,j)-tx(i,j,3)*phi(i+1,j)-tx(i,j,4)*
     +                 phi(nx-1,j))/tx(i,j,2)
	  end do
	end do
c
c     set periodic point
c
	do j=jst,jfn,2
	  phi(nx,j) = phi(1,j)
	end do
c
c     set rhs and solve on lines thru even y points
c
	jfn = ny-1
!$OMP PARALLEL DO SHARED(cofy,phi,rhs,nx,ny) PRIVATE(i,j)
	do j=2,ny-1,2
	  do i=1,nx-1
	    phi(i,j)=rhs(i,j)-cofy(j,1)*phi(i,j-1)-cofy(j,2)*phi(i,j+1)
	  end do
	end do
c
c     set last y point if even and non-specified
c
	if (mod(ny,2).eq.0.and.nyd.ne.1) then
	  jfn = ny
	  j = ny
	  jm1 = ny-1
	  jp1 = 2
	  do i=1,nx-1
	    phi(i,j)=rhs(i,j)-cofy(j,1)*phi(i,jm1)-cofy(j,2)*phi(i,jp1)
	  end do
	end if
c
c     forward sweep
c
!$OMP PARALLEL DO SHARED(sum,jfn,phi,tx,nx) PRIVATE(i,ib,j)
	do j=2,jfn,2
	  do i=2,nx-2
	    phi(i,j) = phi(i,j)-tx(i,j,1)*phi(i-1,j)
	  end do
	  do i=1,nx-2
	    sum(j) = sum(j)+tx(i,j,5)*phi(i,j)
	  end do
	  phi(nx-1,j) = phi(nx-1,j)-sum(j)
c
c     backward sweep
c
	  phi(nx-1,j) = phi(nx-1,j)/tx(nx-1,j,2)
	  phi(nx-2,j) = (phi(nx-2,j)-tx(nx-2,j,4)*phi(nx-1,j))
     +                  /tx(nx-2,j,2)
	  do ib=4,nx
	    i = nx-ib+1
	    phi(i,j) = (phi(i,j)-tx(i,j,3)*phi(i+1,j)-tx(i,j,4)*
     +                 phi(nx-1,j))/tx(i,j,2)
	  end do
	end do
c
c     set periodic point
c
	do j=2,jfn,2
	  phi(nx,j) = phi(1,j)
	end do
c
c     adjust for periodic y if ny even (only possible at coarsest level)
c     red/black line will not allow phi(i,1) = phi(i,ny)
c
	if (nyc.eq.0.and.mod(ny,2).eq.0) then
	  do i=1,nx
	    phi(i,ny) = phi(i,1)
	  end do
	end if
      end if
c
c      final set of periodic and virtual x and y boundaries if necessary
c
      if (nyc.eq.0) then
	do i=1,nx
	  phi(i,0) = phi(i,ny-1)
	  phi(i,ny) = phi(i,1)
	  phi(i,ny+1) = phi(i,2)
	end do
      end if
      if (nxa.eq.0) then
	do j=1,ny
	  phi(0,j) = phi(nx-1,j)
	  phi(nx,j) = phi(1,j)
	  phi(nx+1,j) = phi(2,j)
	end do
      end if
      return
      end

      subroutine slymd2sp(nx,ny,phi,rhs,cofx,cofy,ty,sum)
      implicit none
      integer nx,ny,i,j,jb,ist,ifn,jst,jfn,im1,ip1
      integer intl,nxa,nxb,nyc,nyd,ixp,jyq,iex,jey,nfx,nfy,iguess,
     +             maxcy,method,nwork,lwork,itero,ngrid,klevel,kcur,
     +             kcycle,iprer,ipost,intpol,kps
      common/icud2/intl,nxa,nxb,nyc,nyd,ixp,jyq,iex,jey,nfx,nfy,iguess,
     +             maxcy,method,nwork,lwork,itero,ngrid,klevel,kcur,
     +             kcycle,iprer,ipost,intpol,kps
      real phi(0:nx+1,0:ny+1),rhs(nx,ny)
      real cofx(nx,3),cofy(ny,3),ty(ny,nx,*),sum(nx)
c
c     replace line y with point gauss-seidel if
c     y direction is periodic and ny = 3
c
      if (nyc .eq. 0 .and. ny .eq. 3) then
	call rcd2spp(nx,ny,phi,rhs,cofx,cofy)
	return
      end if
      ist = 3
      ifn = nx-1
      jst = 2
      if (nyc.ne.1) jst = 1
      jfn = ny-1
      if (nyd.ne.1) jfn = ny
      if (nyc.ne.0) then
c
c     non-periodic case
c     set adjusted rhs and solve on odd x points
c
	if (nxa.ne.1) then
	  i = 1
	  ist = 1
	  im1 = nx-1
	  ip1 = 2
	  do j=jst,jfn
	    phi(i,j)=rhs(i,j)-cofx(i,1)*phi(im1,j)-cofx(i,2)*phi(ip1,j)
	  end do
	end if
c
c     x interior
c
C$OMP PARALLEL DO SHARED(jst,jfn,phi,rhs,cofx,nx) PRIVATE(i,j)
	do i=3,nx-1,2
	  do j=jst,jfn
	    phi(i,j)=rhs(i,j)-cofx(i,1)*phi(i-1,j)-cofx(i,2)*phi(i+1,j)
	  end do
	end do
	if (mod(nx,2).ne.0.and.nxb.ne.1) then
	  ifn = nx
	  i = nx
	  im1 = nx-1
	  ip1 = 2
	  do j=jst,jfn
	    phi(i,j)=rhs(i,j)-cofx(i,1)*phi(im1,j)-cofx(i,2)*phi(ip1,j)
	  end do
	end if
c
c     forward sweep
c
C$OMP PARALLEL DO SHARED(ist,ifn,phi,ty,ny) PRIVATE(i,jb,j)
	do i=ist,ifn,2
	  do j=2,ny
	    phi(i,j) = phi(i,j)-ty(j-1,i,1)*phi(i,j-1)
	  end do
c
c     backward sweep
c
	  phi(i,ny) = phi(i,ny)/ty(ny,i,2)
	  do jb=2,ny
	    j = ny-jb+1
	    phi(i,j) = (phi(i,j)-ty(j,i,3)*phi(i,j+1))/ty(j,i,2)
	  end do
	end do
c
c     repeat for even x points
c
	ifn = nx-1
C$OMP PARALLEL DO SHARED(jst,jfn,phi,rhs,cofx,nx) PRIVATE(i,j)
	do i=2,nx-1,2
	  do j=jst,jfn
	    phi(i,j)=rhs(i,j)-cofx(i,1)*phi(i-1,j)-cofx(i,2)*phi(i+1,j)
	  end do
	end do
	if (mod(nx,2).eq.0.and.nxb.ne.1) then
	  ifn = nx
	  i = nx
	  im1 = nx-1
	  ip1 = 2
	  do j=jst,jfn
	    phi(i,j)=rhs(i,j)-cofx(i,1)*phi(im1,j)-cofx(i,2)*phi(ip1,j)
	  end do
	end if
c
c     forward sweep
c
C$OMP PARALLEL DO SHARED(ifn,phi,ty,ny) PRIVATE(i,jb,j)
	do i=2,ifn,2
	  do j=2,ny
	    phi(i,j) = phi(i,j)-ty(j-1,i,1)*phi(i,j-1)
	  end do
c
c     backward sweep
c
	  phi(i,ny) = phi(i,ny)/ty(ny,i,2)
	  do jb=2,ny
	    j = ny-jb+1
	    phi(i,j) = (phi(i,j)-ty(j,i,3)*phi(i,j+1))/ty(j,i,2)
	  end do
	end do
c
c     adjust x boundaries if nxa=0 and nx is even (only possible at coar
c     red/black line ordering above will not allow phi(1,j) = phi(nx,j)
c
	if (nxa.eq.0.and.mod(nx,2).eq.0) then
	  do j=1,ny
	    phi(nx,j) = phi(1,j)
	  end do
	end if
      else
c
c     line periodic relaxation in y direction
c
	do i=1,nx
	  sum(i) = 0.0
	end do
c
c     set rhs and solve on odd x points
c
	ist = 3
	ifn = nx-1
	if (nxa.ne.1) then
	  ist = 1
	  i = 1
	  im1 = nx-1
	  ip1 = 2
	  do j=1,ny-1
	    phi(i,j)=rhs(i,j)-cofx(i,1)*phi(im1,j)-cofx(i,2)*phi(ip1,j)
	  end do
	end if
C$OMP PARALLEL DO SHARED(phi,rhs,cofx,nx,ny) PRIVATE(i,j)
	do i=3,nx-1,2
	  do j=1,ny-1
	    phi(i,j)=rhs(i,j)-cofx(i,1)*phi(i-1,j)-cofx(i,2)*phi(i+1,j)
	  end do
	end do
c
c     set last x point if odd and non-specified
c
	if (mod(nx,2).ne.0.and.nxb.ne.1) then
	  ifn = nx
	  i = nx
	  im1 = nx-1
	  ip1 = 2
	  do j=1,ny-1
	    phi(i,j)=rhs(i,j)-cofx(i,1)*phi(im1,j)-cofx(i,2)*phi(ip1,j)
	  end do
	end if
c
c     forward sweep
c
C$OMP PARALLEL DO SHARED(ist,ifn,phi,ty,ny) PRIVATE(i,jb,j)
	do i=ist,ifn,2
	  do j=2,ny-2
	    phi(i,j) = phi(i,j)-ty(j,i,1)*phi(i,j-1)
	  end do
	  do j=1,ny-2
	    sum(i) = sum(i)+ty(j,i,5)*phi(i,j)
	  end do
	  phi(i,ny-1) = phi(i,ny-1)-sum(i)
c
c     backward sweep
c
	  phi(i,ny-1) = phi(i,ny-1)/ty(ny-1,i,2)
	  phi(i,ny-2) = (phi(i,ny-2)-ty(ny-2,i,4)*phi(i,ny-1))/
     +                   ty(ny-2,i,2)
	  do jb=4,ny
	    j = ny-jb+1
	    phi(i,j) = (phi(i,j)-ty(j,i,3)*phi(i,j+1)-ty(j,i,4)*
     +                 phi(i,ny-1))/ty(j,i,2)
	  end do
	end do
c
c     set periodic point
c
	do i=ist,ifn,2
	  phi(i,ny) = phi(i,1)
	end do
c
c     set rhs and solve on even x points
c
	ifn = nx-1
C$OMP PARALLEL DO SHARED(phi,rhs,cofx,nx,ny) PRIVATE(i,j)
	do i=2,nx-1,2
	  do j=1,ny-1
	    phi(i,j)=rhs(i,j)-cofx(i,1)*phi(i-1,j)-cofx(i,2)*phi(i+1,j)
	  end do
	end do
c
c     set last x line if even and non-specified
c
	if (mod(nx,2).eq.0.and.nxb.ne.1) then
	  ifn = nx
	  i = nx
	  im1 = nx-1
	  ip1 = 2
	  do j=1,ny-1
	    phi(i,j)=rhs(i,j)-cofx(i,1)*phi(im1,j)-cofx(i,2)*phi(ip1,j)
	  end do
	end if
c
c     forward sweep
c
C$OMP PARALLEL DO SHARED(ifn,phi,ty,ny) PRIVATE(i,jb,j)
	do i=2,ifn,2
	  do j=2,ny-2
	    phi(i,j) = phi(i,j)-ty(j,i,1)*phi(i,j-1)
	  end do
	  do j=1,ny-2
	    sum(i) = sum(i)+ty(j,i,5)*phi(i,j)
	  end do
	  phi(i,ny-1) = phi(i,ny-1)-sum(i)
c
c     backward sweep
c
	  phi(i,ny-1) = phi(i,ny-1)/ty(ny-1,i,2)
	  phi(i,ny-2) = (phi(i,ny-2)-ty(ny-2,i,4)*phi(i,ny-1))/
     +                   ty(ny-2,i,2)
	  do jb=4,ny
	    j = ny-jb+1
	    phi(i,j) = (phi(i,j)-ty(j,i,3)*phi(i,j+1)-ty(j,i,4)*
     +                  phi(i,ny-1))/ty(j,i,2)
	  end do
	end do
c
c     set periodic point
c
	do i=2,ifn,2
	  phi(i,ny) = phi(i,1)
	end do
c
c     adjust for periodic x if nx even (only possible at coarsest level)
c     red/black line will not allow phi(1,j) = phi(nx,j)
c
	if (nxa.eq.0.and.mod(nx,2).eq.0) then
	  do j=1,ny
	    phi(nx,j) = phi(1,j)
	  end do
	end if
      end if
c
c      set periodic and virtual x and y boundaries if necessary
c
      if (nyc.eq.0) then
	do i=1,nx
	  phi(i,0) = phi(i,ny-1)
	  phi(i,ny) = phi(i,1)
	  phi(i,ny+1) = phi(i,2)
	end do
      end if
      if (nxa.eq.0) then
	do j=1,ny
	  phi(0,j) = phi(nx-1,j)
	  phi(nx,j) = phi(1,j)
	  phi(nx+1,j) = phi(2,j)
	end do
      end if
      return
      end

c
c     file mudcom.f
c
c  . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
c  .                                                             .
c  .                  copyright (c) 1999 by UCAR                 .
c  .                                                             .
c  .       UNIVERSITY CORPORATION for ATMOSPHERIC RESEARCH       .
c  .                                                             .
c  .                      all rights reserved                    .
c  .                                                             .
c  .                                                             .
c  .                      MUDPACK version 5.0                    .
c  .                                                             .
c  . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
c
c ... author and specialist
c
c          John C. Adams (National Center for Atmospheric Research)
c          email: johnad@ucar.edu, phone: 303-497-1213

c ... For MUDPACK 5.0 information, visit the website:
c     (http://www.scd.ucar.edu/css/software/mudpack)
c
c ... purpose
c
c     mudcom.f is a common subroutines file containing subroutines
c     called by some or all of the real two- and three-dimensional
c     mudpack solvers.  mudcom.f must be loaded with any real mudpack
c     solver.
c
      subroutine swk2(nfx,nfy,phif,rhsf,phi,rhs)
c
c     set phif,rhsf input in arrays which include
c     virtual boundaries for phi (for all 2-d real codes)
c
      implicit none
      integer nfx,nfy,i,j
      real phif(nfx,nfy),rhsf(nfx,nfy)
      real phi(0:nfx+1,0:nfy+1),rhs(nfx,nfy)
      do j=1,nfy
	do i=1,nfx
	  phi(i,j) = phif(i,j)
	  rhs(i,j) = rhsf(i,j)
	end do
      end do
c
c     set virtual boundaries in phi to zero
c
      do j=0,nfy+1
	phi(0,j) = 0.0
	phi(nfx+1,j) = 0.0
      end do
      do i=0,nfx+1
	phi(i,0) = 0.0
	phi(i,nfy+1) = 0.0
      end do
      return
      end

      subroutine trsfc2(nx,ny,phi,rhs,ncx,ncy,phic,rhsc)
c
c     transfer fine grid to coarse grid
c
      implicit none
      integer nx,ny,ncx,ncy,i,j,ic,jc
      real phi(0:nx+1,0:ny+1),rhs(nx,ny)
      real phic(0:ncx+1,0:ncy+1),rhsc(ncx,ncy)
c
c     set virtual boundaries in phic to zero
c
      do jc=0,ncy+1
	phic(0,jc) = 0.0
	phic(ncx+1,jc) = 0.0
      end do
      do ic=0,ncx+1
	phic(ic,0) = 0.0
	phic(ic,ncy+1) = 0.0
      end do
      if (ncx.lt.nx .and. ncy.lt.ny) then
c
c     coarsening in both x and y
c
	do jc=1,ncy
	  j = jc+jc-1
	  do ic=1,ncx
	    i = ic+ic-1
	    phic(ic,jc) = phi(i,j)
	    rhsc(ic,jc) = rhs(i,j)
	  end do
	end do
      else if (ncx.lt.nx .and. ncy.eq.ny) then
c
c     coarsening in x only
c
	do jc=1,ncy
	  j = jc
	  do ic=1,ncx
	    i = ic+ic-1
	    phic(ic,jc) = phi(i,j)
	    rhsc(ic,jc) = rhs(i,j)
	  end do
	end do
      else
c
c     coarsening in y only
c
	do jc=1,ncy
	  j = jc+jc-1
	  do ic=1,ncx
	    i = ic
	    phic(ic,jc) = phi(i,j)
	    rhsc(ic,jc) = rhs(i,j)
	  end do
	end do
      end if
      return
      end

      subroutine res2(nx,ny,resf,ncx,ncy,rhsc,nxa,nxb,nyc,nyd)
      implicit none
      integer nx,ny,ncx,ncy,nxa,nxb,nyc,nyd
      integer i,j,ic,jc,im1,ip1,jm1,jp1,ix,jy
c
c     restrict fine grid residual in resf to coarse grid in rhsc
c     using full weighting for all real 2d codes
c
      real resf(nx,ny),rhsc(ncx,ncy)
c
c     set x,y coarsening integer subscript scales
c
      ix = 1
      if (ncx.eq.nx) ix = 0
      jy = 1
      if (ncy.eq.ny) jy = 0
c
c     restrict on interior
c
      if (ncy.lt.ny .and. ncx.lt.nx) then
c
c     coarsening in both directions
c
!$OMP PARALLEL DO PRIVATE(i,j,ic,jc), SHARED(resf,rhsc,ncx,ncy)
	do jc=2,ncy-1
	  j = jc+jc-1
	  do ic=2,ncx-1
	    i = ic+ic-1
	    rhsc(ic,jc) = (resf(i-1,j-1)+resf(i+1,j-1)+resf(i-1,j+1)+
     +                     resf(i+1,j+1)+2.*(resf(i-1,j)+resf(i+1,j)+
     +                     resf(i,j-1)+resf(i,j+1))+4.*resf(i,j))*.0625
	  end do
	end do
      else if (ncy.eq.ny) then
c
c     no coarsening in y but coarsening in x
c
!$OMP PARALLEL DO PRIVATE(i,j,ic,jc), SHARED(resf,rhsc,ncx,ncy)
	do jc=2,ncy-1
	  j = jc
	  do ic=2,ncx-1
	    i = ic+ic-1
	    rhsc(ic,jc) = (resf(i-1,j-1)+resf(i+1,j-1)+resf(i-1,j+1)+
     +                     resf(i+1,j+1)+2.*(resf(i-1,j)+resf(i+1,j)+
     +                     resf(i,j-1)+resf(i,j+1))+4.*resf(i,j))*.0625
	  end do
	end do
      else
c
c     no coarsening in x but coarsening in y
c
!$OMP PARALLEL DO PRIVATE(i,j,ic,jc), SHARED(resf,rhsc,ncx,ncy)
	do jc=2,ncy-1
	  j = jc+jc-1
	  do ic=2,ncx-1
	    i = ic
	    rhsc(ic,jc) = (resf(i-1,j-1)+resf(i+1,j-1)+resf(i-1,j+1)+
     +                     resf(i+1,j+1)+2.*(resf(i-1,j)+resf(i+1,j)+
     +                     resf(i,j-1)+resf(i,j+1))+4.*resf(i,j))*.0625
	  end do
	end do
      end if
c
c     set residual on boundaries
c
      do jc=1,ncy,ncy-1
c
c     y=yc,yd boundaries
c
	j = jc+jy*(jc-1)
	jm1 = max0(j-1,2)
	jp1 = min0(j+1,ny-1)
	if (j.eq.1 .and. nyc.eq.0) jm1 = ny-1
	if (j.eq.ny .and. nyc.eq.0) jp1 = 2
c
c     y=yc,yd and x=xa,xb cornors
c
	do ic=1,ncx,ncx-1
	  i = ic+ix*(ic-1)
	  im1 = max0(i-1,2)
	  ip1 = min0(i+1,nx-1)
	  if (i.eq.1 .and. nxa.eq.0) im1 = nx-1
	  if (i.eq.nx .and. nxa.eq.0) ip1 = 2
	  rhsc(ic,jc) = (resf(im1,jm1)+resf(ip1,jm1)+resf(im1,jp1)+
     +                   resf(ip1,jp1)+2.*(resf(im1,j)+resf(ip1,j)+
     +                   resf(i,jm1)+resf(i,jp1))+4.*resf(i,j))*.0625
	end do
c
c     set y=yc,yd interior edges
c
	do ic=2,ncx-1
	  i = ic+ix*(ic-1)
	  rhsc(ic,jc) = (resf(i-1,jm1)+resf(i+1,jm1)+resf(i-1,jp1)+
     +                   resf(i+1,jp1)+2.*(resf(i-1,j)+resf(i+1,j)+
     +                   resf(i,jm1)+resf(i,jp1))+4.*resf(i,j))*.0625
	end do
      end do
c
c     set x=xa,xb interior edges
c
      do ic=1,ncx,ncx-1
	i = ic+ix*(ic-1)
	im1 = max0(i-1,2)
	ip1 = min0(i+1,nx-1)
	if (i.eq.1 .and. nxa.eq.0) im1 = nx-1
	if (i.eq.nx .and. nxa.eq.0) ip1 = 2
	do jc=2,ncy-1
	  j = jc+jy*(jc-1)
	  rhsc(ic,jc) = (resf(im1,j-1)+resf(ip1,j-1)+resf(im1,j+1)+
     +                   resf(ip1,j+1)+2.*(resf(im1,j)+resf(ip1,j)+
     +                   resf(i,j-1)+resf(i,j+1))+4.*resf(i,j))*.0625
	end do
      end do
c
c     set coarse grid residual zero on specified boundaries
c
      if (nxa.eq.1) then
	do jc=1,ncy
	  rhsc(1,jc) = 0.0
	end do
      end if
      if (nxb.eq.1) then
	do jc=1,ncy
	  rhsc(ncx,jc) = 0.0
	end do
      end if
      if (nyc.eq.1) then
	do ic=1,ncx
	  rhsc(ic,1) = 0.0
	end do
      end if
      if (nyd.eq.1) then
	do ic=1,ncx
	  rhsc(ic,ncy) = 0.0
	end do
      end if
      return
      end
c
c     prolon2 modified from rgrd2u 11/20/97
c
      subroutine prolon2(ncx,ncy,p,nx,ny,q,nxa,nxb,nyc,nyd,intpol)
      implicit none
      integer ncx,ncy,nx,ny,intpol,nxa,nxb,nyc,nyd
      real p(0:ncx+1,0:ncy+1),q(0:nx+1,0:ny+1)
      integer i,j,jc,ist,ifn,jst,jfn,joddst,joddfn
      ist = 1
      ifn = nx
      jst = 1
      jfn = ny
      joddst = 1
      joddfn = ny
      if (nxa.eq.1) then
	ist = 2
      end if
      if (nxb.eq.1) then
	ifn = nx-1
      end if
      if (nyc.eq.1) then
	jst = 2
	joddst = 3
      end if
      if (nyd.eq.1) then
	jfn = ny-1
	joddfn = ny-2
      end if
      if (intpol.eq.1 .or. ncy.lt.4) then
c
c     linearly interpolate in y
c
	if (ncy .lt. ny) then
c
c     ncy grid is an every other point subset of ny grid
c     set odd j lines interpolating in x and then set even
c     j lines by averaging odd j lines
c
	  do j=joddst,joddfn,2
	    jc = j/2+1
	    call prolon1(ncx,p(0,jc),nx,q(0,j),nxa,nxb,intpol)
	  end do
	  do j=2,jfn,2
	    do i=ist,ifn
	      q(i,j) = 0.5*(q(i,j-1)+q(i,j+1))
	    end do
	  end do
c
c     set periodic virtual boundaries if necessary
c
	  if (nyc.eq.0) then
	    do i=ist,ifn
	      q(i,0) = q(i,ny-1)
	      q(i,ny+1) = q(i,2)
	    end do
	  end if
	  return
	else
c
c     ncy grid equals ny grid so interpolate in x only
c
	  do j=jst,jfn
	    jc = j
	    call prolon1(ncx,p(0,jc),nx,q(0,j),nxa,nxb,intpol)
	  end do
c
c     set periodic virtual boundaries if necessary
c
	  if (nyc.eq.0) then
	    do i=ist,ifn
	      q(i,0) = q(i,ny-1)
	      q(i,ny+1) = q(i,2)
	    end do
	  end if
	  return
	end if
      else
c
c     cubically interpolate in y
c
	if (ncy .lt. ny) then
c
c     set every other point of ny grid by interpolating in x
c
	  do j=joddst,joddfn,2
	    jc = j/2+1
	    call prolon1(ncx,p(0,jc),nx,q(0,j),nxa,nxb,intpol)
	  end do
c
c     set deep interior of ny grid using values just
c     generated and symmetric cubic interpolation in y
c
	  do j=4,ny-3,2
	    do i=ist,ifn
	    q(i,j)=(-q(i,j-3)+9.*(q(i,j-1)+q(i,j+1))-q(i,j+3))*.0625
	    end do
	  end do
c
c     interpolate from q at j=2 and j=ny-1
c
	  if (nyc.ne.0) then
c
c     asymmetric formula near nonperiodic y boundaries
c
	    do i=ist,ifn
	      q(i,2)=(5.*q(i,1)+15.*q(i,3)-5.*q(i,5)+q(i,7))*.0625
	      q(i,ny-1)=(5.*q(i,ny)+15.*q(i,ny-2)-5.*q(i,ny-4)+
     +                    q(i,ny-6))*.0625
	    end do
	  else
c
c     periodicity in y alows symmetric formula near bndys
c
	    do i=ist,ifn
	      q(i,2) = (-q(i,ny-2)+9.*(q(i,1)+q(i,3))-q(i,5))*.0625
	      q(i,ny-1)=(-q(i,ny-4)+9.*(q(i,ny-2)+q(i,ny))-q(i,3))*.0625
	      q(i,ny+1) = q(i,2)
	      q(i,0) = q(i,ny-1)
	    end do
	  end if
	  return
	else
c
c     ncy grid is equals ny grid so interpolate in x only
c
	  do j=jst,jfn
	    jc = j
	    call prolon1(ncx,p(0,jc),nx,q(0,j),nxa,nxb,intpol)
	  end do
c
c     set periodic virtual boundaries if necessary
c
	  if (nyc.eq.0) then
	    do i=ist,ifn
	      q(i,0) = q(i,ny-1)
	      q(i,ny+1) = q(i,2)
	    end do
	  end if
	  return
	end if
      end if
      end

c
c     11/20/97  modification of rgrd1u.f for mudpack
c
      subroutine prolon1(ncx,p,nx,q,nxa,nxb,intpol)
      implicit none
      integer intpol,nxa,nxb,ncx,nx,i,ic,ist,ifn,ioddst,ioddfn
      real p(0:ncx+1),q(0:nx+1)
      ist = 1
      ioddst = 1
      ifn = nx
      ioddfn = nx
      if (nxa.eq.1) then
	ist = 2
	ioddst = 3
      end if
      if (nxb.eq.1) then
	ifn = nx-1
	ioddfn = nx-2
      end if
      if (intpol.eq.1 .or. ncx.lt.4) then
c
c     linear interpolation in x
c
	if (ncx .lt. nx) then
c
c     every other point of nx grid is ncx grid
c
	  do i=ioddst,ioddfn,2
	    ic = (i+1)/2
	    q(i) = p(ic)
	  end do
	  do i=2,ifn,2
	    q(i) = 0.5*(q(i-1)+q(i+1))
	  end do
	else
c
c     nx grid equals ncx grid
c
	  do i=ist,ifn
	    q(i) = p(i)
	  end do
	end if
c
c     set virtual end points if periodic
c
	if (nxa.eq.0) then
	  q(0) = q(nx-1)
	  q(nx+1) = q(2)
	end if
	return
      else
c
c     cubic interpolation in x
c
	if (ncx .lt. nx) then
	  do i=ioddst,ioddfn,2
	    ic = (i+1)/2
	    q(i) = p(ic)
	  end do
c
c      set deep interior with symmetric formula
c
	  do i=4,nx-3,2
	    q(i)=(-q(i-3)+9.*(q(i-1)+q(i+1))-q(i+3))*.0625
	  end do
c
c     interpolate from q at i=2 and i=nx-1
c
	  if (nxa.ne.0) then
c
c     asymmetric formula near nonperiodic bndys
c
	    q(2)=(5.*q(1)+15.*q(3)-5.*q(5)+q(7))*.0625
	    q(nx-1)=(5.*q(nx)+15.*q(nx-2)-5.*q(nx-4)+q(nx-6))*.0625
	  else
c
c     periodicity in x alows symmetric formula near bndys
c
	    q(2) = (-q(nx-2)+9.*(q(1)+q(3))-q(5))*.0625
	    q(nx-1) = (-q(nx-4)+9.*(q(nx-2)+q(nx))-q(3))*.0625
	    q(nx+1) = q(2)
	    q(0) = q(nx-1)
	  end if
	  return
	else
c
c     ncx grid equals nx grid
c
	  do i=ist,ifn
	    q(i) = p(i)
	  end do
	  if (nxa.eq.0) then
	    q(0) = q(nx-1)
	    q(nx+1) = q(2)
	  end if
	  return
	end if
      end if
      end


      subroutine cor2(nx,ny,phif,ncx,ncy,phic,nxa,nxb,nyc,nyd,intpol,
     +                phcor)
c
c     add coarse grid correction in phic to fine grid approximation
c     in phif using linear or cubic interpolation
c
      implicit none
      integer i,j,nx,ny,ncx,ncy,nxa,nxb,nyc,nyd,intpol,ist,ifn,jst,jfn
      real phif(0:nx+1,0:ny+1),phic(0:ncx+1,0:ncy+1)
      real phcor(0:nx+1,0:ny+1)
      do j=0,ny+1
	do i=0,nx+1
	  phcor(i,j) = 0.0
	end do
      end do
c
c     lift correction in phic to fine grid in phcor
c
      call prolon2(ncx,ncy,phic,nx,ny,phcor,nxa,nxb,nyc,nyd,intpol)
c
c     add correction in phcor to phif on nonspecified boundaries
c
      ist = 1
      ifn = nx
      jst = 1
      jfn = ny
      if (nxa.eq.1) ist = 2
      if (nxb.eq.1) ifn = nx-1
      if (nyc.eq.1) jst = 2
      if (nyd.eq.1) jfn = ny-1
      do j=jst,jfn
	do i=ist,ifn
	  phif(i,j) = phif(i,j) + phcor(i,j)
	end do
      end do
c
c     add periodic points if necessary
c
      if (nyc.eq.0) then
	do i=ist,ifn
	  phif(i,0) = phif(i,ny-1)
	  phif(i,ny+1) = phif(i,2)
	end do
      end if
      if (nxa.eq.0) then
	do j=jst,jfn
	  phif(0,j) = phif(nx-1,j)
	  phif(nx+1,j) = phif(2,j)
	end do
      end if
      end

      subroutine pde2(nx,ny,u,i,j,ux3,ux4,uy3,uy4,nxa,nyc)
      implicit none
      integer nx,ny,i,j,nxa,nyc
      real u(nx,ny),dlx,dly,dlxx,dlyy,tdlx3,tdly3,dlx4,dly4
      common/pde2com/dlx,dly,dlxx,dlyy,tdlx3,tdly3,dlx4,dly4
      real ux3,ux4,uy3,uy4
c
c     use second order approximation in u to estimate (second order)
c     third and fourth partial derivatives in the x and y direction
c     non-symmetric difference formula (derived from the  routine
c     finpdf,findif) are used at and one point in from mixed boundaries.
c
      if (nxa.ne.0) then
c
c     nonperiodic in x
c
	if(i.gt.2 .and. i.lt.nx-1) then
	  ux3 = (-u(i-2,j)+2.0*u(i-1,j)-2.0*u(i+1,j)+u(i+2,j))/tdlx3
	  ux4 = (u(i-2,j)-4.0*u(i-1,j)+6.0*u(i,j)-4.0*u(i+1,j)+u(i+2,j))
     +           /dlx4
	else if (i.eq.1) then
	  ux3 = (-5.0*u(1,j)+18.0*u(2,j)-24.0*u(3,j)+14.0*u(4,j)-
     +           3.0*u(5,j))/tdlx3
	  ux4 = (3.0*u(1,j)-14.0*u(2,j)+26.0*u(3,j)-24.0*u(4,j)+
     +           11.0*u(5,j)-2.0*u(6,j))/dlx4
	else if (i.eq.2) then
	  ux3 = (-3.0*u(1,j)+10.0*u(2,j)-12.0*u(3,j)+6.0*u(4,j)-u(5,j))
     +           /tdlx3
	  ux4 = (2.0*u(1,j)-9.0*u(2,j)+16.0*u(3,j)-14.0*u(4,j)+
     +           6.0*u(5,j)-u(6,j))/dlx4
	else if (i.eq.nx-1) then
	  ux3 = (u(nx-4,j)-6.0*u(nx-3,j)+12.0*u(nx-2,j)-10.0*u(nx-1,j)+
     +           3.0*u(nx,j))/tdlx3
	 ux4 = (-u(nx-5,j)+6.0*u(nx-4,j)-14.0*u(nx-3,j)+16.0*u(nx-2,j)-
     +           9.0*u(nx-1,j)+2.0*u(nx,j))/dlx4
	else if (i.eq.nx) then
	  ux3 = (3.0*u(nx-4,j)-14.0*u(nx-3,j)+24.0*u(nx-2,j)-
     +           18.0*u(nx-1,j)+5.0*u(nx,j))/tdlx3
	  ux4 = (-2.0*u(nx-5,j)+11.0*u(nx-4,j)-24.0*u(nx-3,j)+
     +           26.0*u(nx-2,j)-14.0*u(nx-1,j)+3.0*u(nx,j))/dlx4
	end if
      else
c
c     periodic in x
c
	if(i.gt.2 .and. i.lt.nx-1) then
	  ux3 = (-u(i-2,j)+2.0*u(i-1,j)-2.0*u(i+1,j)+u(i+2,j))/tdlx3
	  ux4 = (u(i-2,j)-4.0*u(i-1,j)+6.0*u(i,j)-4.0*u(i+1,j)+u(i+2,j))
     +           /dlx4
	else if (i.eq.1) then
	  ux3 = (-u(nx-2,j)+2.0*u(nx-1,j)-2.0*u(2,j)+u(3,j))/tdlx3
	  ux4 = (u(nx-2,j)-4.0*u(nx-1,j)+6.0*u(1,j)-4.0*u(2,j)+u(3,j))
     +          /dlx4
	else if (i.eq.2) then
	  ux3 = (-u(nx-1,j)+2.0*u(1,j)-2.0*u(3,j)+u(4,j))/(tdlx3)
	  ux4 = (u(nx-1,j)-4.0*u(1,j)+6.0*u(2,j)-4.0*u(3,j)+u(4,j))/dlx4
	else if (i.eq.nx-1) then
	  ux3 = (-u(nx-3,j)+2.0*u(nx-2,j)-2.0*u(1,j)+u(2,j))/tdlx3
	  ux4 = (u(nx-3,j)-4.0*u(nx-2,j)+6.0*u(nx-1,j)-4.0*u(1,j)+
     +           u(2,j))/dlx4
	else if (i.eq.nx) then
	  ux3 = (-u(nx-2,j)+2.0*u(nx-1,j)-2.0*u(2,j)+u(3,j))/tdlx3
	  ux4 = (u(nx-2,j)-4.0*u(nx-1,j)+6.0*u(nx,j)-4.0*u(2,j)+u(3,j))
     +          /dlx4
	end if
      end if
c
c     y partial derivatives
c
      if (nyc.ne.0) then
c
c     not periodic in y
c
	if (j.gt.2 .and. j.lt.ny-1) then
	  uy3 = (-u(i,j-2)+2.0*u(i,j-1)-2.0*u(i,j+1)+u(i,j+2))/tdly3
	  uy4 = (u(i,j-2)-4.0*u(i,j-1)+6.0*u(i,j)-4.0*u(i,j+1)+u(i,j+2))
     +          /dly4
	else if (j.eq.1) then
	  uy3 = (-5.0*u(i,1)+18.0*u(i,2)-24.0*u(i,3)+14.0*u(i,4)-
     +            3.0*u(i,5))/tdly3
	  uy4 = (3.0*u(i,1)-14.0*u(i,2)+26.0*u(i,3)-24.0*u(i,4)+
     +           11.0*u(i,5)-2.0*u(i,6))/dly4
	else if (j.eq.2) then
	  uy3 = (-3.0*u(i,1)+10.0*u(i,2)-12.0*u(i,3)+6.0*u(i,4)-u(i,5))
     +          /tdly3
	  uy4 = (2.0*u(i,1)-9.0*u(i,2)+16.0*u(i,3)-14.0*u(i,4)+
     +           6.0*u(i,5)-u(i,6))/dly4
	else if (j.eq.ny-1) then
	  uy3 = (u(i,ny-4)-6.0*u(i,ny-3)+12.0*u(i,ny-2)-10.0*u(i,ny-1)+
     +           3.0*u(i,ny))/tdly3
	  uy4 = (-u(i,ny-5)+6.0*u(i,ny-4)-14.0*u(i,ny-3)+16.0*u(i,ny-2)-
     +           9.0*u(i,ny-1)+2.0*u(i,ny))/dly4
	else if (j.eq.ny) then
	  uy3 = (3.0*u(i,ny-4)-14.0*u(i,ny-3)+24.0*u(i,ny-2)-
     +           18.0*u(i,ny-1)+5.0*u(i,ny))/tdly3
	  uy4 = (-2.0*u(i,ny-5)+11.0*u(i,ny-4)-24.0*u(i,ny-3)+
     +           26.0*u(i,ny-2)-14.0*u(i,ny-1)+3.0*u(i,ny))/dly4
	end if
      else
c
c     periodic in y
c
	if (j.gt.2 .and. j.lt.ny-1) then
	  uy3 = (-u(i,j-2)+2.0*u(i,j-1)-2.0*u(i,j+1)+u(i,j+2))/tdly3
	  uy4 = (u(i,j-2)-4.0*u(i,j-1)+6.0*u(i,j)-4.0*u(i,j+1)+u(i,j+2))
     +           /dly4
	else if (j.eq.1) then
	  uy3 = (-u(i,ny-2)+2.0*u(i,ny-1)-2.0*u(i,2)+u(i,3))/tdly3
	  uy4 = (u(i,ny-2)-4.0*u(i,ny-1)+6.0*u(i,1)-4.0*u(i,2)+u(i,3))
     +          /dly4
	else if (j.eq.2) then
	  uy3 = (-u(i,ny-1)+2.0*u(i,1)-2.0*u(i,3)+u(i,4))/(tdly3)
	  uy4 = (u(i,ny-1)-4.0*u(i,1)+6.0*u(i,2)-4.0*u(i,3)+u(i,4))/dly4
	else if (j.eq.ny-1) then
	  uy3 = (-u(i,ny-3)+2.0*u(i,ny-2)-2.0*u(i,1)+u(i,2))/tdly3
	  uy4 = (u(i,ny-3)-4.0*u(i,ny-2)+6.0*u(i,ny-1)-4.0*u(i,1)+
     +           u(i,2))/dly4
	else if (j.eq.ny) then
	  uy3 = (-u(i,ny-2)+2.0*u(i,ny-1)-2.0*u(i,2)+u(i,3))/tdly3
	  uy4 = (u(i,ny-2)-4.0*u(i,ny-1)+6.0*u(i,ny)-4.0*u(i,2)+u(i,3))
     +          /dly4
	end if
      end if
      return
      end

      subroutine swk3(nfx,nfy,nfz,phif,rhsf,phi,rhs)
c
c     set phif,rhsf input in arrays which include
c     virtual boundaries for phi (for all 2-d real codes)
c
      implicit none
      integer nfx,nfy,nfz,i,j,k
      real phif(nfx,nfy,nfz),rhsf(nfx,nfy,nfz)
      real phi(0:nfx+1,0:nfy+1,0:nfz+1),rhs(nfx,nfy,nfz)
      do k=1,nfz
	do j=1,nfy
	  do i=1,nfx
	    phi(i,j,k) = phif(i,j,k)
	    rhs(i,j,k) = rhsf(i,j,k)
	  end do
	end do
      end do
c
c     set virtual boundaries in phi to zero
c
      do k=0,nfz+1
	do j=0,nfy+1
	  phi(0,j,k) = 0.0
	  phi(nfx+1,j,k) = 0.0
	end do
      end do
      do k=0,nfz+1
	do i=0,nfx+1
	  phi(i,0,k) = 0.0
	  phi(i,nfy+1,k) = 0.0
	end do
      end do
      do j=0,nfy+1
	do i=0,nfx+1
	  phi(i,j,0) = 0.0
	  phi(i,j,nfz+1) = 0.0
	end do
      end do
      return
      end

      subroutine trsfc3(nx,ny,nz,phi,rhs,ncx,ncy,ncz,phic,rhsc)
c
c     transfer fine grid to coarse grid
c
      implicit none
      integer nx,ny,nz,ncx,ncy,ncz,i,j,k,ic,jc,kc,ix,jy,kz
      real phi(0:nx+1,0:ny+1,0:nz+1),rhs(nx,ny,nz)
      real phic(0:ncx+1,0:ncy+1,0:ncz+1),rhsc(ncx,ncy,ncz)
c
c     set virtual boundaries in phic to zero
c
      do kc=0,ncz+1
	do jc=0,ncy+1
	  phic(0,jc,kc) = 0.0
	  phic(ncx+1,jc,kc) = 0.0
	end do
      end do
      do kc=0,ncz+1
	do ic=0,ncx+1
	  phic(ic,0,kc) = 0.0
	  phic(ic,ncy+1,kc) = 0.0
	end do
      end do
      do jc=0,ncy+1
	do ic=0,ncx+1
	  phic(ic,jc,0) = 0.0
	  phic(ic,jc,ncz+1) = 0.0
	end do
      end do
      if (ncx.lt.nx .and. ncy.lt.ny .and. ncz.lt.nz) then
c
c     coarsening in x,y,z (usually the case?)
c
	do kc=1,ncz
	k = kc+kc-1
	do jc=1,ncy
	  j = jc+jc-1
	  do ic=1,ncx
	    i = ic+ic-1
	    phic(ic,jc,kc) = phi(i,j,k)
	    rhsc(ic,jc,kc) = rhs(i,j,k)
	  end do
	end do
	end do
      else
c
c     no coarsening in at least one dimension
c
	ix = 1
	if (ncx.eq.nx) ix = 0
	jy = 1
	if (ncy.eq.ny) jy = 0
	kz = 1
	if (ncz.eq.nz) kz = 0

	do kc=1,ncz
	  k = kc+kz*(kc-1)
	  do jc=1,ncy
	    j = jc+jy*(jc-1)
	    do ic=1,ncx
	      i = ic+ix*(ic-1)
	      phic(ic,jc,kc) = phi(i,j,k)
	      rhsc(ic,jc,kc) = rhs(i,j,k)
	    end do
	  end do
	end do
      end if
      return
      end

      subroutine res3(nx,ny,nz,resf,ncx,ncy,ncz,rhsc,
     +                nxa,nxb,nyc,nyd,nze,nzf)
      implicit none
      integer nx,ny,nz,ncx,ncy,ncz,nxa,nxb,nyc,nyd,nze,nzf
      integer ix,jy,kz,i,j,k,ic,jc,kc,im1,ip1,jm1,jp1,km1,kp1
      real rm,rk,rp
c
c     restrict fine grid residual in resf to coarse grid in rhsc
c     using full weighting
c
      real resf(nx,ny,nz),rhsc(ncx,ncy,ncz)
c
c     set x,y,z coarsening integer subscript scales
c
      ix = 1
      if (ncx.eq.nx) ix = 0
      jy = 1
      if (ncy.eq.ny) jy = 0
      kz = 1
      if (ncz.eq.nz) kz = 0
c
c     restrict on interior
c
      if (ncz.lt.nz .and. ncy.lt.ny .and. ncx.lt.nx) then
c
c     coarsening in x,y,z
c
!$OMP PARALLEL DO PRIVATE(i,j,k,ic,jc,kc,rm,rk,rp)
!$OMP+SHARED(resf,rhsc,ncx,ncy,ncz)
      do kc=2,ncz-1
	k = kc+kc-1
	do jc=2,ncy-1
	  j = jc+jc-1
	  do ic=2,ncx-1
	    i = ic+ic-1
c
c     weight on k-1,k,k+1 z planes in rm,rk,rp
c
	    rm=(resf(i-1,j-1,k-1)+resf(i+1,j-1,k-1)+resf(i-1,j+1,k-1)+
     +      resf(i+1,j+1,k-1)+2.*(resf(i-1,j,k-1)+resf(i+1,j,k-1)+
     +      resf(i,j-1,k-1)+resf(i,j+1,k-1))+4.*resf(i,j,k-1))*.0625
	    rk=(resf(i-1,j-1,k)+resf(i+1,j-1,k)+resf(i-1,j+1,k)+
     +      resf(i+1,j+1,k)+2.*(resf(i-1,j,k)+resf(i+1,j,k)+
     +      resf(i,j-1,k)+resf(i,j+1,k))+4.*resf(i,j,k))*.0625
	    rp=(resf(i-1,j-1,k+1)+resf(i+1,j-1,k+1)+resf(i-1,j+1,k+1)+
     +      resf(i+1,j+1,k+1)+2.*(resf(i-1,j,k+1)+resf(i+1,j,k+1)+
     +      resf(i,j-1,k+1)+resf(i,j+1,k+1))+4.*resf(i,j,k+1))*.0625
c
c     weight in z direction for final result
c
	    rhsc(ic,jc,kc) = 0.25*(rm+2.*rk+rp)
	  end do
	end do
      end do
      else
c
c     allow for noncoarsening in any of x,y,z
c
!$OMP PARALLEL DO PRIVATE(i,j,k,ic,jc,kc,rm,rk,rp)
!$OMP+SHARED(ix,jy,kz,resf,rhsc,ncx,ncy,ncz)
      do kc=2,ncz-1
	k = kc+kz*(kc-1)
	do jc=2,ncy-1
	  j = jc+jy*(jc-1)
	  do ic=2,ncx-1
	    i = ic+ix*(ic-1)
c
c     weight on k-1,k,k+1 z planes in rm,rk,rp
c
	    rm=(resf(i-1,j-1,k-1)+resf(i+1,j-1,k-1)+resf(i-1,j+1,k-1)+
     +      resf(i+1,j+1,k-1)+2.*(resf(i-1,j,k-1)+resf(i+1,j,k-1)+
     +      resf(i,j-1,k-1)+resf(i,j+1,k-1))+4.*resf(i,j,k-1))*.0625
	    rk=(resf(i-1,j-1,k)+resf(i+1,j-1,k)+resf(i-1,j+1,k)+
     +      resf(i+1,j+1,k)+2.*(resf(i-1,j,k)+resf(i+1,j,k)+
     +      resf(i,j-1,k)+resf(i,j+1,k))+4.*resf(i,j,k))*.0625
	    rp=(resf(i-1,j-1,k+1)+resf(i+1,j-1,k+1)+resf(i-1,j+1,k+1)+
     +      resf(i+1,j+1,k+1)+2.*(resf(i-1,j,k+1)+resf(i+1,j,k+1)+
     +      resf(i,j-1,k+1)+resf(i,j+1,k+1))+4.*resf(i,j,k+1))*.0625
c
c     weight in z direction for final result
c
	    rhsc(ic,jc,kc) = 0.25*(rm+2.*rk+rp)
	  end do
	end do
      end do
      end if
c
c     set residual on boundaries
c
      do ic=1,ncx,ncx-1
c
c     x=xa and x=xb
c
	i = ic+ix*(ic-1)
	im1 = max0(i-1,2)
	ip1 = min0(i+1,nx-1)
	if (i.eq.1 .and. nxa.eq.0) im1 = nx-1
	if (i.eq.nx .and. nxb.eq.0) ip1 = 2
c
c    (y,z) interior
c
!$OMP PARALLEL DO PRIVATE(j,k,jc,kc,rm,rk,rp)
!$OMP+SHARED(kz,jy,ic,im1,i,ip1,resf,rhsc,ncy,ncz)
	do kc=2,ncz-1
	  k = kc+kz*(kc-1)
	  do jc=2,ncy-1
	    j = jc+jy*(jc-1)
	    rm=(resf(im1,j-1,k-1)+resf(ip1,j-1,k-1)+resf(im1,j+1,k-1)+
     +      resf(ip1,j+1,k-1)+2.*(resf(im1,j,k-1)+resf(ip1,j,k-1)+
     +      resf(i,j-1,k-1)+resf(i,j+1,k-1))+4.*resf(i,j,k-1))*.0625
	    rk=(resf(im1,j-1,k)+resf(ip1,j-1,k)+resf(im1,j+1,k)+
     +      resf(ip1,j+1,k)+2.*(resf(im1,j,k)+resf(ip1,j,k)+
     +      resf(i,j-1,k)+resf(i,j+1,k))+4.*resf(i,j,k))*.0625
	    rp=(resf(im1,j-1,k+1)+resf(ip1,j-1,k+1)+resf(im1,j+1,k+1)+
     +      resf(ip1,j+1,k+1)+2.*(resf(im1,j,k+1)+resf(ip1,j,k+1)+
     +      resf(i,j-1,k+1)+resf(i,j+1,k+1))+4.*resf(i,j,k+1))*.0625
	    rhsc(ic,jc,kc) = 0.25*(rm+2.*rk+rp)
	  end do
	end do
c
c     x=xa,xb and y=yc,yd interior edges
c
	do jc=1,ncy,ncy-1
	  j = jc+jy*(jc-1)
	  jm1 = max0(j-1,2)
	  jp1 = min0(j+1,ny-1)
	  if (j.eq.1 .and. nyc.eq.0) jm1 = ny-1
	  if (j.eq.ny .and. nyc.eq.0) jp1 = 2
	  do kc=2,ncz-1
	    k = kc+kz*(kc-1)
	    rm=(resf(im1,jm1,k-1)+resf(ip1,jm1,k-1)+resf(im1,jp1,k-1)+
     +      resf(ip1,jp1,k-1)+2.*(resf(im1,j,k-1)+resf(ip1,j,k-1)+
     +      resf(i,jm1,k-1)+resf(i,jp1,k-1))+4.*resf(i,j,k-1))*.0625
	    rk=(resf(im1,jm1,k)+resf(ip1,jm1,k)+resf(im1,jp1,k)+
     +      resf(ip1,jp1,k)+2.*(resf(im1,j,k)+resf(ip1,j,k)+
     +      resf(i,jm1,k)+resf(i,jp1,k))+4.*resf(i,j,k))*.0625
	    rp=(resf(im1,jm1,k+1)+resf(ip1,jm1,k+1)+resf(im1,jp1,k+1)+
     +      resf(ip1,jp1,k+1)+2.*(resf(im1,j,k+1)+resf(ip1,j,k+1)+
     +      resf(i,jm1,k+1)+resf(i,jp1,k+1))+4.*resf(i,j,k+1))*.0625
	    rhsc(ic,jc,kc) = 0.25*(rm+2.*rk+rp)
	  end do
c     x=xa,xb; y=yc,yd; z=ze,zf cornors
	  do kc=1,ncz,ncz-1
	  k = kc+kz*(kc-1)
	  km1 = max0(k-1,2)
	  kp1 = min0(k+1,nz-1)
	  if (k.eq.1 .and. nze.eq.0) km1 = nz-1
	  if (k.eq.nz .and. nzf.eq.0) kp1 = 2
	  rm=(resf(im1,jm1,km1)+resf(ip1,jm1,km1)+resf(im1,jp1,km1)+
     +    resf(ip1,jp1,km1)+2.*(resf(im1,j,km1)+resf(ip1,j,km1)+
     +    resf(i,jm1,km1)+resf(i,jp1,km1))+4.*resf(i,j,km1))*.0625
	  rk=(resf(im1,jm1,k)+resf(ip1,jm1,k)+resf(im1,jp1,k)+
     +    resf(ip1,jp1,k)+2.*(resf(im1,j,k)+resf(ip1,j,k)+
     +    resf(i,jm1,k)+resf(i,jp1,k))+4.*resf(i,j,k))*.0625
	  rp=(resf(im1,jm1,kp1)+resf(ip1,jm1,kp1)+resf(im1,jp1,kp1)+
     +    resf(ip1,jp1,kp1)+2.*(resf(im1,j,kp1)+resf(ip1,j,kp1)+
     +    resf(i,jm1,kp1)+resf(i,jp1,kp1))+4.*resf(i,j,kp1))*.0625
	  rhsc(ic,jc,kc) = 0.25*(rm+2.*rk+rp)
	  end do
	end do
c
c      x=xa,xb and z=ze,zf edges
c
	do kc=1,ncz,ncz-1
	  k = kc+kz*(kc-1)
	  km1 = max0(k-1,2)
	  kp1 = min0(k+1,nz-1)
	  if (k.eq.1 .and. nze.eq.0) km1 = nz-1
	  if (k.eq.nz .and. nzf.eq.0) kp1 = 2
	   do jc=2,ncy-1
	    j = jc+jy*(jc-1)
	    rm=(resf(im1,j-1,km1)+resf(ip1,j-1,km1)+resf(im1,j+1,km1)+
     +      resf(ip1,j+1,km1)+2.*(resf(im1,j,km1)+resf(ip1,j,km1)+
     +      resf(i,j-1,km1)+resf(i,j+1,km1))+4.*resf(i,j,km1))*.0625
	    rk=(resf(im1,j-1,k)+resf(ip1,j-1,k)+resf(im1,j+1,k)+
     +      resf(ip1,j+1,k)+2.*(resf(im1,j,k)+resf(ip1,j,k)+
     +      resf(i,j-1,k)+resf(i,j+1,k))+4.*resf(i,j,k))*.0625
	    rp=(resf(im1,j-1,kp1)+resf(ip1,j-1,kp1)+resf(im1,j+1,kp1)+
     +      resf(ip1,j+1,kp1)+2.*(resf(im1,j,kp1)+resf(ip1,j,kp1)+
     +      resf(i,j-1,kp1)+resf(i,j+1,kp1))+4.*resf(i,j,kp1))*.0625
	    rhsc(ic,jc,kc) = 0.25*(rm+2.*rk+rp)
	  end do
	end do
      end do
c
c     y boundaries y=yc and y=yd
c
      do jc=1,ncy,ncy-1
	j = jc+jy*(jc-1)
	jm1 = max0(j-1,2)
	jp1 = min0(j+1,ny-1)
	if (j.eq.1 .and. nyc.eq.0) jm1 = ny-1
	if (j.eq.ny .and. nyd.eq.0) jp1 = 2
c
c     (x,z) interior
c
!$OMP PARALLEL DO PRIVATE(i,k,ic,kc,rm,rk,rp)
!$OMP+SHARED(ix,kz,jc,jm1,j,jp1,resf,rhsc,ncx,ncz)
	do kc=2,ncz-1
	  k = kc+kz*(kc-1)
	  do ic=2,ncx-1
	    i = ic+ix*(ic-1)
	    rm=(resf(i-1,jm1,k-1)+resf(i+1,jm1,k-1)+resf(i-1,jp1,k-1)+
     +      resf(i+1,jp1,k-1)+2.*(resf(i-1,j,k-1)+resf(i+1,j,k-1)+
     +      resf(i,jm1,k-1)+resf(i,jp1,k-1))+4.*resf(i,j,k-1))*.0625
	    rk=(resf(i-1,jm1,k)+resf(i+1,jm1,k)+resf(i-1,jp1,k)+
     +      resf(i+1,jp1,k)+2.*(resf(i-1,j,k)+resf(i+1,j,k)+
     +      resf(i,jm1,k)+resf(i,jp1,k))+4.*resf(i,j,k))*.0625
	    rp=(resf(i-1,jm1,k+1)+resf(i+1,jm1,k+1)+resf(i-1,jp1,k+1)+
     +      resf(i+1,jp1,k+1)+2.*(resf(i-1,j,k+1)+resf(i+1,j,k+1)+
     +      resf(i,jm1,k+1)+resf(i,jp1,k+1))+4.*resf(i,j,k+1))*.0625
	    rhsc(ic,jc,kc) = 0.25*(rm+2.*rk+rp)
	  end do
	end do
c
c     y=yc,yd and z=ze,zf edges
c
	do kc=1,ncz,ncz-1
	  k = kc+kz*(kc-1)
	  km1 = max0(k-1,2)
	  kp1 = min0(k+1,nz-1)
	  if (k.eq.1 .and. nze.eq.0) km1 = nz-1
	  if (k.eq.nz .and. nzf.eq.0) kp1 = 2
c
c     interior in x
c
	  do ic=2,ncx-1
	    i = ic+ix*(ic-1)
	    rm=(resf(i-1,jm1,km1)+resf(i+1,jm1,km1)+resf(i-1,jp1,km1)+
     +      resf(i+1,jp1,km1)+2.*(resf(i-1,j,km1)+resf(i+1,j,km1)+
     +      resf(i,jm1,km1)+resf(i,jp1,km1))+4.*resf(i,j,km1))*.0625
	    rk=(resf(i-1,jm1,k)+resf(i+1,jm1,k)+resf(i-1,jp1,k)+
     +      resf(i+1,jp1,k)+2.*(resf(i-1,j,k)+resf(i+1,j,k)+
     +      resf(i,jm1,k)+resf(i,jp1,k))+4.*resf(i,j,k))*.0625
	    rp=(resf(i-1,jm1,kp1)+resf(i+1,jm1,kp1)+resf(i-1,jp1,kp1)+
     +      resf(i+1,jp1,kp1)+2.*(resf(i-1,j,kp1)+resf(i+1,j,kp1)+
     +      resf(i,jm1,kp1)+resf(i,jp1,kp1))+4.*resf(i,j,kp1))*.0625
	    rhsc(ic,jc,kc) = 0.25*(rm+2.*rk+rp)
	  end do
	end do
      end do
c
c     z=ze,zf boundaries
c
      do kc=1,ncz,ncz-1
	k = kc+kz*(kc-1)
	km1 = max0(k-1,2)
	kp1 = min0(k+1,nz-1)
	if (k.eq.1 .and. nze.eq.0) km1 = nz-1
	if (k.eq.nz .and. nzf.eq.0) kp1 = 2
c
c     (x,y) interior
c
!$OMP PARALLEL DO PRIVATE(i,j,ic,jc,rm,rk,rp)
!$OMP+SHARED(ix,jy,kc,km1,k,kp1,resf,rhsc,ncx,ncz)
	do jc=2,ncy-1
	  j = jc+jy*(jc-1)
	  do ic=2,ncx-1
	    i = ic+ix*(ic-1)
	    rm=(resf(i-1,j-1,km1)+resf(i+1,j-1,km1)+resf(i-1,j+1,km1)+
     +      resf(i+1,j+1,km1)+2.*(resf(i-1,j,km1)+resf(i+1,j,km1)+
     +      resf(i,j-1,km1)+resf(i,j+1,km1))+4.*resf(i,j,km1))*.0625
	    rk=(resf(i-1,j-1,k)+resf(i+1,j-1,k)+resf(i-1,j+1,k)+
     +      resf(i+1,j+1,k)+2.*(resf(i-1,j,k)+resf(i+1,j,k)+
     +      resf(i,j-1,k)+resf(i,j+1,k))+4.*resf(i,j,k))*.0625
	    rp=(resf(i-1,j-1,kp1)+resf(i+1,j-1,kp1)+resf(i-1,j+1,kp1)+
     +      resf(i+1,j+1,kp1)+2.*(resf(i-1,j,kp1)+resf(i+1,j,kp1)+
     +      resf(i,j-1,kp1)+resf(i,j+1,kp1))+4.*resf(i,j,kp1))*.0625
	    rhsc(ic,jc,kc) = 0.25*(rm+2.*rk+rp)
	  end do
	end do
      end do
c
c     set coarse grid residual to zero at specified boundaries
c
      if (nxa.eq.1) then
	ic = 1
	do kc=1,ncz
	  do jc=1,ncy
	    rhsc(ic,jc,kc) = 0.0
	  end do
	end do
      end if
      if (nxb.eq.1) then
	ic = ncx
	do kc=1,ncz
	  do jc=1,ncy
	    rhsc(ic,jc,kc) = 0.0
	  end do
	end do
      end if
      if (nyc.eq.1) then
	jc = 1
	do kc=1,ncz
	  do ic=1,ncx
	    rhsc(ic,jc,kc) = 0.0
	  end do
	end do
      end if
      if (nyd.eq.1) then
	jc = ncy
	do kc=1,ncz
	  do ic=1,ncx
	    rhsc(ic,jc,kc) = 0.0
	  end do
	end do
      end if
      if (nze.eq.1) then
	kc = 1
	do jc=1,ncy
	  do ic=1,ncx
	    rhsc(ic,jc,kc) = 0.0
	  end do
	end do
      end if
      if (nzf.eq.1) then
	kc = ncz
	do jc=1,ncy
	  do ic=1,ncx
	    rhsc(ic,jc,kc) = 0.0
	  end do
	end do
      end if
      return
      end

c
c     prolon3 modified from prolon2 11/25/97
c
      subroutine prolon3(ncx,ncy,ncz,p,nx,ny,nz,q,nxa,nxb,nyc,nyd,
     +                   nze,nzf,intpol)
      implicit none
      integer ncx,ncy,ncz,nx,ny,nz,intpol,nxa,nxb,nyc,nyd,nze,nzf
      real p(0:ncx+1,0:ncy+1,0:ncz+1),q(0:nx+1,0:ny+1,0:nz+1)
      integer i,j,k,kc,ist,ifn,jst,jfn,kst,kfn,koddst,koddfn
      ist = 1
      ifn = nx
      jst = 1
      jfn = ny
      kst = 1
      kfn = nz
      koddst = 1
      koddfn = nz
      if (nxa.eq.1) then
	ist = 2
      end if
      if (nxb.eq.1) then
	ifn = nx-1
      end if
      if (nyc.eq.1) then
	jst = 2
      end if
      if (nyd.eq.1) then
	jfn = ny-1
      end if
      if (nze.eq.1) then
	kst = 2
	koddst = 3
      end if
      if (nzf.eq.1) then
	kfn = nz-1
	koddfn = nz-2
      end if
      if (intpol.eq.1 .or. ncz.lt.4) then
c
c     linearly interpolate in z
c
	if (ncz .lt. nz) then
c
c     ncz grid is an every other point subset of nz grid
c     set odd k planes interpolating in x&y and then set even
c     k planes by averaging odd k planes
c
	  do k=koddst,koddfn,2
	    kc = k/2+1
	    call prolon2(ncx,ncy,p(0,0,kc),nx,ny,q(0,0,k),nxa,nxb,nyc,
     +                   nyd,intpol)
	  end do
	  do k=2,kfn,2
	    do j=jst,jfn
	      do i=ist,ifn
		q(i,j,k) = 0.5*(q(i,j,k-1)+q(i,j,k+1))
	      end do
	    end do
	  end do
c
c     set periodic virtual boundaries if necessary
c
	  if (nze.eq.0) then
	    do j=jst,jfn
	      do i=ist,ifn
		q(i,j,0) = q(i,j,nz-1)
		q(i,j,nz+1) = q(i,j,2)
	      end do
	    end do
	  end if
	  return
	else
c
c     ncz grid is equals nz grid so interpolate in x&y only
c
	  do k=kst,kfn
	    kc = k
	    call prolon2(ncx,ncy,p(0,0,kc),nx,ny,q(0,0,k),nxa,nxb,nyc,
     +                   nyd,intpol)
	  end do
c
c     set periodic virtual boundaries if necessary
c
	  if (nze.eq.0) then
	    do j=jst,jfn
	    do i=ist,ifn
	      q(i,j,0) = q(i,j,nz-1)
	      q(i,j,nz+1) = q(i,j,2)
	    end do
	    end do
	  end if
	  return
	end if
      else
c
c     cubically interpolate in z
c
	if (ncz .lt. nz) then
c
c     set every other point of nz grid by interpolating in x&y
c
	  do k=koddst,koddfn,2
	    kc = k/2+1
	    call prolon2(ncx,ncy,p(0,0,kc),nx,ny,q(0,0,k),nxa,nxb,nyc,
     +                   nyd,intpol)
	  end do
c
c     set deep interior of nz grid using values just
c     generated and symmetric cubic interpolation in z
c
	  do k=4,nz-3,2
	    do j=jst,jfn
	    do i=ist,ifn
	    q(i,j,k)=(-q(i,j,k-3)+9.*(q(i,j,k-1)+q(i,j,k+1))-q(i,j,k+3))
     +                *.0625
	    end do
	    end do
	  end do
c
c     interpolate from q at k=2 and k=nz-1
c
	  if (nze.ne.0) then
c
c     asymmetric formula near nonperiodic z boundaries
c
	    do j=jst,jfn
	    do i=ist,ifn
	      q(i,j,2)=(5.*q(i,j,1)+15.*q(i,j,3)-5.*q(i,j,5)+q(i,j,7))
     +                  *.0625
	      q(i,j,nz-1)=(5.*q(i,j,nz)+15.*q(i,j,nz-2)-5.*q(i,j,nz-4)+
     +                    q(i,j,nz-6))*.0625
	    end do
	    end do
	  else
c
c     periodicity in y alows symmetric formula near bndys
c
	    do j=jst,jfn
	    do i=ist,ifn
	      q(i,j,2) = (-q(i,j,nz-2)+9.*(q(i,j,1)+q(i,j,3))-q(i,j,5))
     +                   *.0625
	      q(i,j,nz-1)=(-q(i,j,nz-4)+9.*(q(i,j,nz-2)+q(i,j,nz))-
     +                      q(i,j,3))*.0625
	      q(i,j,nz+1) = q(i,j,2)
	      q(i,j,0) = q(i,j,nz-1)
	    end do
	    end do
	  end if
	  return
	else
c
c     ncz grid is equals nx grid so interpolate in x&y only
c
	  do k=kst,kfn
	    kc = k
	    call prolon2(ncx,ncy,p(0,0,kc),nx,ny,q(0,0,k),nxa,nxb,nyc,
     +                   nyd,intpol)
	  end do
c
c     set periodic virtual boundaries if necessary
c
	  if (nze.eq.0) then
	    do j=jst,jfn
	    do i=ist,ifn
	      q(i,j,0) = q(i,j,nz-1)
	      q(i,j,nz+1) = q(i,j,2)
	    end do
	    end do
	  end if
	  return
	end if
      end if
      end

      subroutine cor3(nx,ny,nz,phif,ncx,ncy,ncz,phic,nxa,nxb,nyc,nyd,
     +                nze,nzf,intpol,phcor)
      implicit none
      integer nx,ny,nz,ncx,ncy,ncz,nxa,nxb,nyc,nyd,nze,nzf,intpol
      integer i,j,k,ist,ifn,jst,jfn,kst,kfn
c
c     add coarse grid correction in phic to fine grid approximation
c     in phif using linear or cubic interpolation
c
      real phif(0:nx+1,0:ny+1,0:nz+1),phic(0:ncx+1,0:ncy+1,0:ncz+1)
      real phcor(0:nx+1,0:ny+1,0:nz+1)
      do k=0,nz+1
	do j=0,ny+1
	  do i=0,nx+1
	    phcor(i,j,k) = 0.0
	  end do
	end do
      end do
c
c     lift correction in phic to fine grid in phcor
c
      call prolon3(ncx,ncy,ncz,phic,nx,ny,nz,phcor,nxa,nxb,nyc,nyd,
     +             nze,nzf,intpol)
c
c     add correction in phcor to phif on nonspecified boundaries
c
      ist = 1
      ifn = nx
      jst = 1
      jfn = ny
      kst = 1
      kfn = nz
      if (nxa.eq.1) ist = 2
      if (nxb.eq.1) ifn = nx-1
      if (nyc.eq.1) jst = 2
      if (nyd.eq.1) jfn = ny-1
      if (nze.eq.1) kst = 2
      if (nzf.eq.1) kfn = nz-1
      do k=kst,kfn
	do j=jst,jfn
	  do i=ist,ifn
	    phif(i,j,k) = phif(i,j,k) + phcor(i,j,k)
	  end do
	end do
      end do
c
c     add periodic points if necessary
c
      if (nze.eq.0) then
	do j=jst,jfn
	  do i=ist,ifn
	    phif(i,j,0) = phif(i,j,nz-1)
	    phif(i,j,nz+1) = phif(i,j,2)
	  end do
	end do
      end if
      if (nyc.eq.0) then
	do k=kst,kfn
	  do i=ist,ifn
	    phif(i,0,k) = phif(i,ny-1,k)
	    phif(i,ny+1,k) = phif(i,2,k)
	  end do
	end do
      end if
      if (nxa.eq.0) then
	do k=kst,kfn
	  do j=jst,jfn
	    phif(0,j,k) = phif(nx-1,j,k)
	    phif(nx+1,j,k) = phif(2,j,k)
	  end do
	end do
      end if
      end

      subroutine per3vb(nx,ny,nz,phi,nxa,nyc,nze)
c
c     set virtual periodic boundaries from interior values
c     in three dimensions (for all 3-d solvers)
c
      implicit none
      integer nx,ny,nz,nxa,nyc,nze,j,k,i
      real phi(0:nx+1,0:ny+1,0:nz+1)
      if (nxa.eq.0) then
	do k=1,nz
	  do j=1,ny
	    phi(0,j,k) = phi(nx-1,j,k)
	    phi(nx,j,k) = phi(1,j,k)
	    phi(nx+1,j,k) = phi(2,j,k)
	  end do
	end do
      end if
      if (nyc.eq.0) then
	do k=1,nz
	  do i=1,nx
	    phi(i,0,k) = phi(i,ny-1,k)
	    phi(i,ny,k) = phi(i,1,k)
	    phi(i,ny+1,k) = phi(i,2,k)
	  end do
	end do
      end if
      if (nze.eq.0) then
	do j=1,ny
	  do i=1,nx
	    phi(i,j,0) = phi(i,j,nz-1)
	    phi(i,j,nz) = phi(i,j,1)
	    phi(i,j,nz+1) = phi(i,j,2)
	  end do
	end do
      end if
      return
      end

      subroutine pde2cr(nx,ny,u,i,j,ux3y,uxy3,ux2y2)
c
c     compute mixed partial derivative approximations
c
      implicit none
      integer nx,ny,i,j,n1,n2,n3,n4,m1,m2,m3,m4
      real u(nx,ny),ux3y,uxy3,ux2y2
      integer intl,nxa,nxb,nyc,nyd,ixp,jyq,iex,jey,nfx,nfy,iguess,
     +             maxcy,method,nwork,lwork,itero,ngrid,klevel,kcur,
     +             kcycle,iprer,ipost,intpol,kps
      real xa,xb,yc,yd,tolmax,relmax
      common/imud2cr/intl,nxa,nxb,nyc,nyd,ixp,jyq,iex,jey,nfx,nfy,
     +               iguess, maxcy,method,nwork,lwork,itero,ngrid,
     +               klevel,kcur,kcycle,iprer,ipost,intpol,kps
      common/fmud2cr/xa,xb,yc,yd,tolmax,relmax
      real dlx,dly,dyox,dxoy,dlx2,dly2,dlxx,dlxy,dlyy,dlxy2,
     +             dlxy4,dxxxy4,dxyyy4,dxxyy,tdlx3,tdly3,dlx4,dly4,
     +             dlxxx,dlyyy
      common/com2dcr/dyox,dxoy,dlx2,dly2,dlxy,dlxy2,dlxy4,
     +               dxxxy4,dxyyy4,dxxyy,dlxxx,dlyyy
      common/pde2com/dlx,dly,dlxx,dlyy,tdlx3,tdly3,dlx4,dly4
      n1=ny-1
      n2=ny-2
      n3=ny-3
      n4=ny-4
      m1=nx-1
      m2=nx-2
      m3=nx-3
      m4=nx-4

      if (i.eq.1) then

      if ((j.gt.2.and.j.lt.ny-1)) then
c     x=xa, yinterior
      ux3y=(5*u(1,j-1)-18*u(2,j-1)+24*u(3,j-1)-14*u(4,j-1)+3*u(5,j-1)
     +     -5*u(1,j+1)+18*u(2,j+1)-24*u(3,j+1)+14*u(4,j+1)-3*u(5,j+1))
     + /dxxxy4
      uxy3=(3*u(1,j-2)-4*u(2,j-2)+u(3,j-2)
     +     -6*u(1,j-1)+8*u(2,j-1)-2*u(3,j-1)
     +     +6*u(1,j+1)-8*u(2,j+1)+2*u(3,j+1)
     +    -3*u(1,j+2)+4*u(2,j+2)-u(3,j+2))/dxyyy4
      else if (j.eq.1) then
c     (xa,yc)
      ux3y=(15*u(1,1)-54*u(2,1)+72*u(3,1)-42*u(4,1)+9*u(5,1)
     + -20*u(1,2)+72*u(2,2)-96*u(3,2)+56*u(4,2)-12*u(5,2)
     + +5*u(1,3)-18*u(2,3)+24*u(3,3)-14*u(4,3)+3*u(5,3))
     + /dxxxy4
      uxy3=(15*u(1,1)-20*u(2,1)+5*u(3,1)
     +           -54*u(1,2)+72*u(2,2)-18*u(3,2)
     +           +72*u(1,3)-96*u(2,3)+24*u(3,3)
     +           -42*u(1,4)+56*u(2,4)-14*u(3,4)
     +           +9*u(1,5)-12*u(2,5)+3*u(3,5))
     + /dxyyy4
      ux2y2=(4*u(1,1)-10*u(2,1)+8*u(3,1)-2*u(4,1)
     +           -10*u(1,2)+25*u(2,2)-20*u(3,2)+5*u(4,2)
     +           +8*u(1,3)-20*u(2,3)+16*u(3,3)-4*u(4,3)
     +           -2*u(1,4)+5*u(2,4)-4*u(3,4)+u(4,4))
     + /dxxyy
      else if (j.eq.2) then
c     (xa,yc+dly)
      ux3y=(5*u(1,1)-18*u(2,1)+24*u(3,1)-14*u(4,1)+3*u(5,1)
     +           -5*u(1,3)+18*u(2,3)-24*u(3,3)+14*u(4,3)-3*u(5,3))
     + /dxxxy4
      uxy3=(9*u(1,1)-12*u(2,1)+3*u(3,1)
     +      -30*u(1,2)+40*u(2,2)-10*u(3,2)
     +      +36*u(1,3)-48*u(2,3)+12*u(3,3)
     +      -18*u(1,4)+24*u(2,4)-6*u(3,4)
     +      +3*u(1,5)-4*u(2,5)+u(3,5))
     + /dxyyy4
      else if (j.eq.ny-1) then
c     x=xa,y=yd-dly
      ux3y=(5*u(1,j-1)-18*u(2,j-1)+24*u(3,j-1)-14*u(4,j-1)+3*u(5,j-1)
     +     -5*u(1,j+1)+18*u(2,j+1)-24*u(3,j+1)+14*u(4,j+1)-3*u(5,j+1))
      uxy3=(5*u(1,n2)-18*u(2,n2)+24*u(3,n2)-14*u(4,n2)+3*u(5,n2)
     +   -5*u(1,ny)+18*u(2,ny)-24*u(3,ny)+14*u(4,ny)-3*u(5,ny))
     + /dxyyy4
      else if (j.eq.ny) then
c     x=xa, y=yd
      ux3y=(-5*u(1,n2)+18*u(2,n2)-24*u(3,n2)+14*u(4,n2)-3*u(5,n2)
     +     +20*u(1,n1)-72*u(2,n1)+96*u(3,n1)-56*u(4,n1)+12*u(5,n1)
     +   -15*u(1,ny)+54*u(2,ny)-72*u(3,ny)+42*u(4,ny)-9*u(5,ny))
     + /dxxxy4
      uxy3=(-9*u(1,n4)+12*u(2,n4)-3*u(3,n4)
     +      +42*u(1,n3)-56*u(2,n3)+14*u(3,n3)
     +      -72*u(1,n2)+96*u(2,n2)-24*u(3,n2)
     +      +54*u(1,n1)-72*u(2,n1)+18*u(3,n1)
     +      -15*u(1,ny)+20*u(2,ny)-5*u(3,ny))
     + /dxyyy4
      ux2y2=(-2*u(1,n3)+5*u(2,n3)-4*u(3,n3)+u(4,n3)
     +      +8*u(1,n2)-20*u(2,n2)+16*u(3,n2)-4*u(4,n2)
     +           -10*u(1,n1)+25*u(2,n1)-20*u(3,n1)+5*u(4,n1)
     +      +4*u(1,ny)-10*u(2,ny)+8*u(3,ny)-2*u(4,ny))
     + /dxxyy
      end if

      else if (i.eq.2) then

      if ((j.gt.2.and.j.lt.ny-1)) then
c     x=xa+dlx, y interior
      ux3y=(3*u(1,j-1)-10*u(2,j-1)+12*u(3,j-1)-6*u(4,j-1)+u(5,j-1)
     +-3*u(1,j+1)+10*u(2,j+1)-12*u(3,j+1)+6*u(4,j+1)-u(5,j+1))/dxxxy4
      uxy3=(u(1,j-2)-u(3,j-2)-2*u(1,j-1)+2*u(3,j-1)
     +     +2*u(1,j+1)-2*u(3,j+1)-u(1,j+2)+u(3,j+2))/dxyyy4
      else if (j.eq.1) then
c     x=xa+dlx, y=yc
      ux3y=(9*u(1,1)-30*u(2,1)+36*u(3,1)-18*u(4,1)+3*u(5,1)
     +      -12*u(1,2)+40*u(2,2)-48*u(3,2)+24*u(4,2)-4*u(5,2)
     +      +3*u(1,3)-10*u(2,3)+12*u(3,3)-6*u(4,3)+u(5,3))
     + /dxxxy4
      uxy3=(5*u(1,1)-5*u(3,1)-18*u(1,2)+18*u(3,2)
     +      +24*u(1,3)-24*u(3,3)-14*u(1,4)
     +      +14*u(3,4)+3*u(1,5)-3*u(3,5))
     + /dxyyy4
      else if (j.eq.2) then
c     at x=xa+dlx,y=yc+dly
      ux3y=(3*u(1,1)-10*u(2,1)+12*u(3,1)-6*u(4,1)+u(5,1)
     +      -3*u(1,3)+10*u(2,3)-12*u(3,3)+6*u(4,3)-u(5,3))
     + /dxxxy4
      uxy3=(3*u(1,1)-3*u(3,1)-10*u(1,2)+10*u(3,2)
     +      +12*u(1,3)-12*u(3,3)-6*u(1,4)+6*u(3,4)
     +      +u(1,5)-u(3,5))
     + /dxyyy4
      else if (j.eq.ny-1) then
c     x=xa+dlx,y=yd-dly
      ux3y=(3*u(1,n2)-10*u(2,n2)+12*u(3,n2)-6*u(4,n2)+u(5,n2)
     +      -3*u(1,ny)+10*u(2,ny)-12*u(3,ny)+6*u(4,ny)-u(5,ny))
     + /dxxxy4
      uxy3=(-u(1,n4)+u(3,n4)+6*u(1,n3)-6*u(3,n3)
     +      -12*u(1,n2)+12*u(3,n2)+10*u(1,n1)-10*u(3,n1)
     +      -3*u(1,ny)+3*u(3,ny))
     + /dxyyy4
      else if (j.eq.ny) then
c     at x=xa+dlx,y=yd
      ux3y=(-3*u(1,n2)+10*u(2,n2)-12*u(3,n2)+6*u(4,n2)-u(5,n2)
     +  +12*u(1,n1)-40*u(2,n1)+48*u(3,n1)-24*u(4,n1)+4*u(5,n1)
     +  -9*u(1,ny)+30*u(2,ny)-36*u(3,ny)+18*u(4,ny)-3*u(5,ny))
     + /dxxxy4
      uxy3=(-3*u(1,n4)+3*u(3,n4)+14*u(1,n3)-14*u(3,n3)
     +      -24*u(1,n2)+24*u(3,n2)+18*u(1,n1)-18*u(3,n1)
     +      -5*u(1,ny)+5*u(3,ny))
     + /dxyyy4
      end if

      else if (i.gt.2 .and. i.lt.nx-1) then

      if (j.eq.1) then
c     y=yc,x interior
      ux3y=(3.0*u(i-2,1)-6.0*u(i-1,1)+6.0*u(i+1,1)-3.0*u(i+2,1)
     +       -4.0*u(i-2,2)+8.0*u(i-1,2)-8.0*u(i+1,2)+4.0*u(i+2,2)
     +       +u(i-2,3)-2.0*u(i-1,3)+2.0*u(i+1,3)-u(i+2,3))
     + /dxxxy4
      uxy3=(5.0*u(i-1,1)-5.0*u(i+1,1)-18.0*u(i-1,2)+18.0*u(i+1,2)
     +   +24.0*u(i-1,3)-24.0*u(i+1,3)-14.0*u(i-1,4)+14.0*u(i+1,4)
     +   +3.0*u(i-1,5)-3.0*u(i+1,5))
     + /dxyyy4
      else if (j.eq.2) then
c     y=yc+dly,x interior
      ux3y=(u(i-2,1)-2.0*u(i-1,1)+2.0*u(i+1,1)-u(i+2,1)
     +      -u(i-2,3)+2.0*u(i-1,3)-2.0*u(i+1,3)+u(i+2,3))
     + /dxxxy4
      uxy3=(u(i-1,1)-u(i+1,1)-2.0*u(i-1,2)+2.0*u(i+1,2)
     +      +2.0*u(i-1,4)-2.0*u(i+1,4)-u(i-1,5)+u(i+1,5))
     + /dxyyy4
      else if (j.eq.ny-1) then
c     y=yd-dly, x interior
      ux3y=(u(i-2,n2)-2.0*u(i-1,n2)+2.0*u(i+1,n2)-u(i+2,n2)
     +      -u(i-2,ny)+2.0*u(i-1,ny)-2.0*u(i+1,ny)+u(i+2,ny))
     + /dxxxy4
      uxy3=(-u(i-1,n4)+u(i+1,n4)+6.0*u(i-1,n3)-6.0*u(i+1,n3)
     +  -12.0*u(i-1,n2)+12.0*u(i+1,n2)+10.0*u(i-1,n1)-10.0*u(i+1,n1)
     +  -3.0*u(i-1,ny)+3.0*u(i+1,ny))
     + /dxyyy4
      else if (j.eq.ny) then
c     at y=yd, x interior
      ux3y=(-u(i-2,n2)+2.0*u(i-1,n2)-2.0*u(i+1,n2)+u(i+2,n2)
     + +4.0*u(i-2,n1)-8.0*u(i-1,n1)+8.0*u(i+1,n1)-4.0*u(i+2,n1)
     + -3.0*u(i-2,ny)+6.0*u(i-1,ny)-6.0*u(i+1,ny)+3.0*u(i+2,ny))
     + /dxxxy4
      uxy3=(-3.0*u(i-1,n4)+3.0*u(i+1,n4)+14.0*u(i-1,n3)-14.0*u(i+1,n3)
     +  -24.0*u(i-1,n2) +24.0*u(i+1,n2)+18.0*u(i-1,n1)-18.0*u(i+1,n1)
     +  -5.0*u(i-1,ny)+5.0*u(i+1,ny))
     +/dxyyy4
      end if

      else if (i.eq.nx-1) then

      if ((j.gt.2.and.j.lt.ny-1)) then
c     x=xb-dlx,y interior
      ux3y=(-u(m4,j-1)+6.*u(m3,j-1)-12.*u(m2,j-1)+10.*u(m1,j-1)-3.*u(nx
     +,j-1)+u(m4,j+1)-6.*u(m3,j+1)+12.*u(m2,j+1)-10.*u(m1,j+1)+3.*u(nx,j
     ++1)) /dxxxy4
      uxy3=(u(m2,j-2)-u(nx,j-2)-2.*u(m2,j-1)+2.*u(nx,j-1)
     + +2.*u(m2,j+1)-2.*u(nx,j+1)-u(m2,j+2)+u(nx,j+2)) /dxyyy4
      else if (j.eq.1) then
c     at x=xb-dlx, y=yc
      ux3y=(-3.0*u(m4,1)+18.0*u(m3,1)-36.0*u(m2,1)+30.0*u(m1,1)-9.0*u(
     +nx,1)+4.0*u(m4,2)-24.0*u(m3,2)+48.0*u(m2,2)-40.0*u(m1,2)+12.0*u(nx
     +,2)-u(m4,3)+6.0*u(m3,3)-12.0*u(m2,3)+10.0*u(m1,3)-3.0*u(nx,3))
     + /dxxxy4
      uxy3=(5.0*u(m2,1)-5.0*u(nx,1)-18.0*u(m2,2)+18.0*u(nx,2)
     + +24.0*u(m2,3)-24.0*u(nx,3)-14.0*u(m2,4)+14.0*u(nx,4)
     + +3.0*u(m2,5)-3.0*u(nx,5))
     + /dxyyy4
      else if (j.eq.2) then
c     x=xb-dlx,y=yc+dly
      ux3y=(-u(m4,1)+6.0*u(m3,1)-12.0*u(m2,1)+10.*u(m1,1)-3.*u(nx,1)
     +         +u(m4,3)-6.0*u(m3,3)+12.0*u(m2,3)-10.*u(m1,3)+3.*u(nx,3))
     + /dxxxy4
      uxy3=(3.0*u(m2,1)-3.*u(nx,1)-10.*u(m2,2)+10.*u(nx,2)
     + +12.*u(m2,3)-12.*u(nx,3)-6.*u(m2,4)+6.*u(nx,4)
     + +u(m2,5)-u(nx,5)) / dxyyy4
      else if (j.eq.ny-1) then
c     at x=xb-dlx,y=yd-dly
      ux3y=(-u(m4,n2)+6.*u(m3,n2)-12.*u(m2,n2)+10.*u(m1,n2)-3.*u(nx,n2)
     + +u(m4,ny)-6.*u(m3,ny)+12.*u(m2,ny)-10.*u(m1,ny)+3.*u(nx,ny))
     + /dxxxy4
      uxy3=(-u(m2,n4)+u(nx,n4)+6*u(m2,n3)-6.*u(nx,n3)
     + -12.*u(m2,n2)+12.*u(nx,n2)+10.*u(m2,n1)-10.*u(nx,n1)
     + -3.*u(m2,ny)+3.*u(nx,ny)) / dxyyy4
      else if (j.eq.ny) then
c     at x=xb.dlx,y=yd
      ux3y=(u(m4,n2)-6.*u(m3,n2)+12.*u(m2,n2)-10.*u(m1,n2)+3.*u(nx,n2)
     + -4.*u(m4,n1)+24.*u(m3,n1)-48.*u(m2,n1)+40.*u(m1,n1)-12.*u(nx,n1)
     + +3.*u(m4,ny)-18.*u(m3,ny)+36.*u(m2,ny)-30.*u(m1,ny)+9.*u(nx,ny))
     + / dxxxy4
      uxy3=(-3.*u(m2,n4)+3.*u(nx,n4)+14.*u(m2,n3)-14.*u(nx,n3)
     + -24.*u(m2,n2)+24.*u(nx,n2)+18.*u(m2,n1)-18.*u(nx,n1)
     + -5.*u(m2,ny)+5.*u(nx,ny)) / dxyyy4
      end if

      else if (i.eq.nx) then

      if ((j.gt.2.and.j.lt.ny-1)) then
c     x=xb,y interior
      ux3y=(-3.*u(m4,j-1)+14.*u(m3,j-1)-24.*u(m2,j-1)+18.*u(m1,j-1)-5.*
     +u(nx,j-1)+3.*u(m4,j+1)-14.*u(m3,j+1)+24.*u(m2,j+1)-18.*u(m1,j+1)+5
     +.*u(nx,j+1)) / dxxxy4
      uxy3=(-u(m2,j-2)+4.*u(m1,j-2)-3.*u(nx,j-2)
     +  +2.*u(m2,j-1)-8.*u(m1,j-1)+6.*u(nx,j-1)
     +  -2.*u(m2,j+1)+8.*u(m1,j+1)-6.*u(nx,j+1)
     +  +u(m2,j+2)-4.*u(m1,j+2)+3.*u(nx,j+2)) / dxyyy4
      else if (j.eq.1) then
c     x=xb,y=yc
      ux3y=(-9.*u(m4,1)+42.*u(m3,1)-72.*u(m2,1)+54.*u(m1,1)-15.*u(nx,1)
     + +12.*u(m4,2)-56.*u(m3,2)+96.*u(m2,2)-72.*u(m1,2)+20.*u(nx,2)
     + -3.*u(m4,3)+14.*u(m3,3)-24.*u(m2,3)+18.*u(m1,3)-5.*u(nx,3))
     + /dxxxy4
      uxy3=(-5.*u(m2,1)+20.*u(m1,1)-15.*u(nx,1)
     +  +18.*u(m2,2)-72.*u(m1,2)+54.*u(nx,2)
     +  -24.*u(m2,3)+96.*u(m1,3)-72.*u(nx,3)
     +  +14.*u(m2,4)-56.*u(m1,4)+42.*u(nx,4)
     +  -3.*u(m2,5)+12.*u(m1,5)-9.*u(nx,5)) / dxyyy4
      ux2y2=(-2.*u(m3,1)+8.*u(m2,1)-10.*u(m1,1)+4.*u(nx,1)
     + +5.*u(m3,2)-20.*u(m2,2)+25.*u(m1,2)-10.*u(nx,2)
     + -4.*u(m3,3)+16.*u(m2,3)-20.*u(m1,3)+8.*u(nx,3)
     + +u(m3,4)-4.*u(m2,4)+5.*u(m1,4)-2.*u(nx,4)) / dxxyy
      else if (j.eq.2) then
c     x=xb,y=yc+dly
      ux3y=(-3.*u(m4,1)+14.*u(m3,1)-24.*u(m2,1)+18.*u(m1,1)-5.*u(nx,1)
     + +3.*u(m4,3)-14.*u(m3,3)+24.*u(m2,3)-18.*u(m1,3)+5.*u(nx,3))
     + / dxxxy4
      uxy3=(-3.*u(m2,1)+12.*u(m1,1)-9.*u(nx,1)
     + +10.*u(m2,2)-40.*u(m1,2)+30.*u(nx,2)
     + -12.*u(m2,3)+48.*u(m1,3)-36.*u(nx,3)
     + +6.*u(m2,4)-24.*u(m1,4)+18.*u(nx,4)
     + -u(m2,5)+4.*u(m1,5)-3.*u(nx,5)) / dxyyy4
      else if (j.eq.ny-1) then
c     x=xb,y=yd-dly
      ux3y=(-3.*u(m4,n2)+14.*u(m3,n2)-24.*u(m2,n2)+18.*u(m1,n2)-5.*u(nx
     +,n2)+3.*u(m4,ny)-14.*u(m3,ny)+24.*u(m2,ny)-18.*u(m1,ny)+5.*u(nx,ny
     +)) / dxxxy4
      uxy3=(u(m2,n4)-4.*u(m1,n4)+3.*u(nx,n4)
     + -6.*u(m2,n3)+24.*u(m1,n3)-18.*u(nx,n3)
     + +12.*u(m2,n2)-48.*u(m1,n2)+36.*u(nx,n2)
     + -10.*u(m2,n1)+40.*u(m1,n1)-30.*u(nx,n1)
     + +3.*u(m2,ny)-12.*u(m1,ny)+9.*u(nx,ny)) / dxyyy4
      else if (j.eq.ny) then
c     x=xb,y=yd
      ux3y=(3.*u(m4,n2)-14.*u(m3,n2)+24.*u(m2,n2)-18.*u(m1,n2)+5.*u(nx,
     +n2)-12.*u(m4,n1)+56.*u(m3,n1)-96.*u(m2,n1)+72.*u(m1,n1)-20.*u(nx,
     +n1)+9.*u(m4,ny)-42.*u(m3,ny)+72.*u(m2,ny)-54.*u(m1,ny)+15.*u(nx,ny
     +)) / dxxxy4
      uxy3=(3.*u(m2,n4)-12.*u(m1,n4)+9.*u(nx,n4)
     + -14.*u(m2,n3)+56.*u(m1,n3)-42.*u(nx,n3)
     + +24.*u(m2,n2)-96.*u(m1,n2)+72.*u(nx,n2)
     + -18.*u(m2,n1)+72.*u(m1,n1)-54.*u(nx,n1)
     + +5.*u(m2,ny)-20.*u(m1,ny)+15.*u(nx,ny)) / dxyyy4
      ux2y2=(u(m3,n3)-4.*u(m2,n3)+5.*u(m1,n3)-2.*u(nx,n3)
     + -4.*u(m3,n2)+16.*u(m2,n2)-20.*u(m1,n2)+8.*u(nx,n2)
     + +5.0*u(m3,n1)-20.*u(m2,n1)+25.*u(m1,n1)-10.*u(nx,n1)
     + -2.*u(m3,ny)+8.*u(m2,ny)-10.*u(m1,ny)+4.*u(nx,ny))
     + / dxxyy
      end if

      end if

      return
      end
      subroutine pde3(nx,ny,nz,u,i,j,k,ux3,ux4,uy3,uy4,uz3,uz4,
     +                nxa,nyc,nze)
c
c     estimate third and fourth partial derivatives in x,y,z
c
      implicit none
      integer nx,ny,nz,i,j,k,nxa,nyc,nze
      real u(nx,ny,nz)
      real dlx,dly,dlz,dlxx,dlyy,dlzz,tdlx3,tdly3,tdlz3,dlx4,dly4,dlz4
      common/pde3com/dlx,dly,dlz,dlxx,dlyy,dlzz,tdlx3,tdly3,tdlz3,
     +               dlx4,dly4,dlz4
      real ux3,ux4,uy3,uy4,uz3,uz4
c
c     x,y partial derivatives
c
      call p3de2(nx,ny,u(1,1,k),i,j,ux3,ux4,uy3,uy4,nxa,nyc)
c
c     z partial derivatives
c
      if (nze.ne.0) then
c
c     nonperiodic in z
c
      if(k.gt.2 .and. k.lt.nz-1) then
      uz3=(-u(i,j,k-2)+2.0*u(i,j,k-1)-2.0*u(i,j,k+1)+u(i,j,k+2))/tdlz3
      uz4=(u(i,j,k-2)-4.0*u(i,j,k-1)+6.0*u(i,j,k)-4.0*u(i,j,k+1)+
     +         u(i,j,k+2))/dlz4
      else if (k.eq.1) then
      uz3=(-5.0*u(i,j,1)+18.0*u(i,j,2)-24.0*u(i,j,3)+14.0*u(i,j,4)-
     +          3.0*u(i,j,5))/tdlz3
      uz4 = (3.0*u(i,j,1)-14.0*u(i,j,2)+26.0*u(i,j,3)-24.0*u(i,j,4)+
     +           11.0*u(i,j,5)-2.0*u(i,j,6))/dlz4
      else if (k.eq.2) then
      uz3 = (-3.0*u(i,j,1)+10.0*u(i,j,2)-12.0*u(i,j,3)+6.0*u(i,j,4)-
     +            u(i,j,5))/tdlz3
      uz4 = (2.0*u(i,j,1)-9.0*u(i,j,2)+16.0*u(i,j,3)-14.0*u(i,j,4)+6.0*
     +           u(i,j,5)-u(i,j,6))/dlz4
      else if (k.eq.nz-1) then
      uz3 = (u(i,j,nz-4)-6.0*u(i,j,nz-3)+12.0*u(i,j,nz-2)-10.0*
     +           u(i,j,nz-1)+3.0*u(i,j,nz))/tdlz3
      uz4 = (-u(i,j,nz-5)+6.0*u(i,j,nz-4)-14.0*u(i,j,nz-3)+16.0*
     +            u(i,j,nz-2)-9.0*u(i,j,nz-1)+2.0*u(i,j,nz))/dlz4
      else if (k.eq.nz) then
      uz3 = (3.0*u(i,j,nz-4)-14.0*u(i,j,nz-3)+24.0*u(i,j,nz-2)-18.0*
     +           u(i,j,nz-1)+5.0*u(i,j,nz))/tdlz3
      uz4 = (-2.0*u(i,j,nz-5)+11.0*u(i,j,nz-4)-24.0*u(i,j,nz-3)+26.0*
     +           u(i,j,nz-2)-14.0*u(i,j,nz-1)+3.0*u(i,j,nz))/dlz4
      end if
      else
c
c     periodic in z so use symmetric formula even "near" z boundaies
c
      if(k.gt.2 .and. k.lt.nz-1) then
      uz3=(-u(i,j,k-2)+2.0*u(i,j,k-1)-2.0*u(i,j,k+1)+u(i,j,k+2))/tdlz3
      uz4=(u(i,j,k-2)-4.0*u(i,j,k-1)+6.0*u(i,j,k)-4.0*u(i,j,k+1)+
     +     u(i,j,k+2))/dlz4
      else if (k.eq.1) then
      uz3 = (-u(i,j,nz-2)+2.0*u(i,j,nz-1)-2.0*u(i,j,2)+u(i,j,3))/tdlz3
      uz4 = (u(i,j,nz-2)-4.0*u(i,j,nz-1)+6.0*u(i,j,1)-4.0*u(i,j,2)+
     +       u(i,j,3))/dlz4
      else if (k.eq.2) then
      uz3 = (-u(i,j,nz-1)+2.0*u(i,j,1)-2.0*u(i,j,3)+u(i,j,4))/(tdlz3)
      uz4 = (u(i,j,nz-1)-4.0*u(i,j,1)+6.0*u(i,j,2)-4.0*u(i,j,3)+
     +       u(i,j,4))/dlz4
      else if (k.eq.nz-1) then
      uz3 = (-u(i,j,nz-3)+2.0*u(i,j,nz-2)-2.0*u(i,j,1)+u(i,j,2))/tdlz3
      uz4 = (u(i,j,nz-3)-4.0*u(i,j,nz-2)+6.0*u(i,j,nz-1)-4.0*u(i,j,1)+
     +       u(i,j,2))/ dlz4
      else if (k.eq.nz) then
      uz3 = (-u(i,j,nz-2)+2.0*u(i,j,nz-1)-2.0*u(i,j,2)+u(i,j,3))/tdlz3
      uz4 = (u(i,j,nz-2)-4.0*u(i,j,nz-1)+6.0*u(i,j,nz)-4.0*u(i,j,2)+
     +       u(i,j,3))/dlz4
      end if
      end if
      return
      end

      subroutine p3de2(nx,ny,u,i,j,ux3,ux4,uy3,uy4,nxa,nyc)
c
c     third and fourth partial derivatives in x and y
c
      implicit none
      integer nx,ny,i,j,nxa,nyc,l
      real u(nx,ny)
      real dlx,dly,dlz,dlxx,dlyy,dlzz,tdlx3,tdly3,tdlz3,dlx4,dly4,dlz4
      common/pde3com/dlx,dly,dlz,dlxx,dlyy,dlzz,tdlx3,tdly3,tdlz3,
     +               dlx4,dly4,dlz4
      real ux3,ux4,uy3,uy4
      l=ny
c
c     x partial derivatives
c
      call p3de1(nx,u(1,j),i,ux3,ux4,nxa)
c
c     y partial derivatives
c
      if (nyc.ne.0) then
c
c     not periodic in y
c
      if (j.gt.2 .and. j.lt.ny-1) then
      uy3 = (-u(i,j-2)+2.0*u(i,j-1)-2.0*u(i,j+1)+u(i,j+2))/tdly3
      uy4 = (u(i,j-2)-4.0*u(i,j-1)+6.0*u(i,j)-4.0*u(i,j+1)+u(i,j+2))/
     +       dly4
      else if (j.eq.1) then
      uy3 = (-5.0*u(i,1)+18.0*u(i,2)-24.0*u(i,3)+14.0*u(i,4)-
     +        3.0*u(i,5))/tdly3
      uy4 = (3.0*u(i,1)-14.0*u(i,2)+26.0*u(i,3)-24.0*u(i,4)+
     +       11.0*u(i,5)-2.0*u(i,6))/dly4
      else if (j.eq.2) then
      uy3 = (-3.0*u(i,1)+10.0*u(i,2)-12.0*u(i,3)+6.0*u(i,4)-u(i,5))/
     +       tdly3
      uy4 = (2.0*u(i,1)-9.0*u(i,2)+16.0*u(i,3)-14.0*u(i,4)+6.0*u(i,5)-
     +       u(i,6))/dly4
      else if (j.eq.ny-1) then
      uy3 = (u(i,l-4)-6.0*u(i,l-3)+12.0*u(i,l-2)-10.0*u(i,l-1)+
     +       3.0*u(i,l))/tdly3
      uy4 = (-u(i,l-5)+6.0*u(i,l-4)-14.0*u(i,l-3)+16.0*u(i,l-2)-
     +       9.0*u(i,l-1)+2.0*u(i,l))/dly4
      else if (j.eq.ny) then
      uy3 = (3.0*u(i,l-4)-14.0*u(i,l-3)+24.0*u(i,l-2)-18.0*u(i,l-1)+
     +       5.0*u(i,l))/tdly3
      uy4 = (-2.0*u(i,l-5)+11.0*u(i,l-4)-24.0*u(i,l-3)+26.0*u(i,l-2)-
     +       14.0*u(i,l-1)+3.0*u(i,l))/dly4
      end if
      else
c
c     periodic in y
c
      if (j.gt.2 .and. j.lt.ny-1) then
      uy3 = (-u(i,j-2)+2.0*u(i,j-1)-2.0*u(i,j+1)+u(i,j+2))/tdly3
      uy4 = (u(i,j-2)-4.0*u(i,j-1)+6.0*u(i,j)-4.0*u(i,j+1)+u(i,j+2))/
     +        dly4
      else if (j.eq.1) then
      uy3 = (-u(i,l-2)+2.0*u(i,l-1)-2.0*u(i,2)+u(i,3))/tdly3
      uy4 = (u(i,l-2)-4.0*u(i,l-1)+6.0*u(i,1)-4.0*u(i,2)+u(i,3))/dly4
      else if (j.eq.2) then
      uy3 = (-u(i,l-1)+2.0*u(i,1)-2.0*u(i,3)+u(i,4))/(tdly3)
      uy4 = (u(i,l-1)-4.0*u(i,1)+6.0*u(i,2)-4.0*u(i,3)+u(i,4))/dly4
      else if (j.eq.ny-1) then
      uy3 = (-u(i,l-3)+2.0*u(i,l-2)-2.0*u(i,1)+u(i,2))/tdly3
      uy4 = (u(i,l-3)-4.0*u(i,l-2)+6.0*u(i,l-1)-4.0*u(i,1)+u(i,2))/
     +        dly4
      else if (j.eq.ny) then
      uy3 = (-u(i,l-2)+2.0*u(i,l-1)-2.0*u(i,2)+u(i,3))/tdly3
      uy4 = (u(i,l-2)-4.0*u(i,l-1)+6.0*u(i,l)-4.0*u(i,2)+u(i,3))/dly4
      end if
      end if
      return
      end

      subroutine p3de1(nx,u,i,ux3,ux4,nxa)
c
c     third and fourth derivatives in x
c
      implicit none
      integer nx,i,nxa,k
      real u(nx)
      real dlx,dly,dlz,dlxx,dlyy,dlzz,tdlx3,tdly3,tdlz3,dlx4,dly4,dlz4
      common/pde3com/dlx,dly,dlz,dlxx,dlyy,dlzz,tdlx3,tdly3,tdlz3,
     +               dlx4,dly4,dlz4
      real ux3,ux4
      k = nx
      if (nxa.ne.0) then
c
c     nonperiodic in x
c
      if(i.gt.2 .and. i.lt.nx-1) then
      ux3 = (-u(i-2)+2.0*u(i-1)-2.0*u(i+1)+u(i+2))/tdlx3
      ux4 = (u(i-2)-4.0*u(i-1)+6.0*u(i)-4.0*u(i+1)+u(i+2))/dlx4
      else if (i.eq.1) then
      ux3 = (-5.0*u(1)+18.0*u(2)-24.0*u(3)+14.0*u(4)-3.0*u(5))/tdlx3
      ux4 = (3.0*u(1)-14.0*u(2)+26.0*u(3)-24.0*u(4)+11.0*u(5)-2.0*u(6))
     +       /dlx4
      else if (i.eq.2) then
      ux3 = (-3.0*u(1)+10.0*u(2)-12.0*u(3)+6.0*u(4)-u(5))/tdlx3
      ux4 = (2.0*u(1)-9.0*u(2)+16.0*u(3)-14.0*u(4)+6.0*u(5)-u(6))/dlx4
      else if (i.eq.nx-1) then
      ux3 = (u(k-4)-6.0*u(k-3)+12.0*u(k-2)-10.0*u(k-1)+3.0*u(k))/tdlx3
      ux4 = (-u(k-5)+6.0*u(k-4)-14.0*u(k-3)+16.0*u(k-2)-9.0*u(k-1)+
     +      2.0*u(k))/dlx4
      else if (i.eq.nx) then
      ux3 = (3.0*u(k-4)-14.0*u(k-3)+24.0*u(k-2)-18.0*u(k-1)+5.0*u(k))/
     +       tdlx3
      ux4 = (-2.0*u(k-5)+11.0*u(k-4)-24.0*u(k-3)+26.0*u(k-2)-
     +       14.0*u(k-1)+3.0*u(k))/dlx4
      end if
      else
c
c     periodic in x
c
      if(i.gt.2 .and. i.lt.nx-1) then
      ux3 = (-u(i-2)+2.0*u(i-1)-2.0*u(i+1)+u(i+2))/tdlx3
      ux4 = (u(i-2)-4.0*u(i-1)+6.0*u(i)-4.0*u(i+1)+u(i+2))/dlx4
      else if (i.eq.1) then
      ux3 = (-u(k-2)+2.0*u(k-1)-2.0*u(2)+u(3))/tdlx3
      ux4 = (u(k-2)-4.0*u(k-1)+6.0*u(1)-4.0*u(2)+u(3))/dlx4
      else if (i.eq.2) then
      ux3 = (-u(k-1)+2.0*u(1)-2.0*u(3)+u(4))/(tdlx3)
      ux4 = (u(k-1)-4.0*u(1)+6.0*u(2)-4.0*u(3)+u(4))/dlx4
      else if (i.eq.nx-1) then
      ux3 = (-u(k-3)+2.0*u(k-2)-2.0*u(1)+u(2))/tdlx3
      ux4 = (u(k-3)-4.0*u(k-2)+6.0*u(k-1)-4.0*u(1)+u(2))/dlx4
      else if (i.eq.nx) then
      ux3 = (-u(k-2)+2.0*u(k-1)-2.0*u(2)+u(3))/tdlx3
      ux4 = (u(k-2)-4.0*u(k-1)+6.0*u(k)-4.0*u(2)+u(3))/dlx4
      end if
      end if
      return
      end
c
c
c     factri and factrip are:
c     subroutines called by any real mudpack solver which uses line
c     relaxation(s) within multigrid iteration.  these subroutines do
c     a vectorized factorization of m simultaneous tridiagonal systems
c     of order n arising from nonperiodic or periodic discretizations
c
      subroutine factri(m,n,a,b,c)
c
c     factor the m simultaneous tridiagonal systems of order n
c
      implicit none
      integer m,n,i,j
      real a(n,m),b(n,m),c(n,m)
      do i=2,n
	do j=1,m
	  a(i-1,j) = a(i-1,j)/b(i-1,j)
	  b(i,j) = b(i,j)-a(i-1,j)*c(i-1,j)
       end do
      end do
      return
      end

      subroutine factrp(m,n,a,b,c,d,e,sum)
c
c     factor the m simultaneous "tridiagonal" systems of order n
c     from discretized periodic system (leave out periodic n point)
c     (so sweeps below only go from i=1,2,...,n-1) n > 3 is necessary
c
      implicit none
      integer m,n,i,j
      real a(n,m),b(n,m),c(n,m),d(n,m),e(n,m),sum(m)
      do j=1,m
	d(1,j) = a(1,j)
      end do
      do i=2,n-2
	do j=1,m
	  a(i,j) = a(i,j)/b(i-1,j)
	  b(i,j) = b(i,j)-a(i,j)*c(i-1,j)
	  d(i,j) = -a(i,j)*d(i-1,j)
       end do
      end do
c
c     correct computation of last d element
c
      do j=1,m
	d(n-2,j) = c(n-2,j)+d(n-2,j)
      end do
      do j=1,m
	e(1,j) = c(n-1,j)/b(1,j)
      end do
      do i=2,n-3
	do j=1,m
	  e(i,j) = -e(i-1,j)*c(i-1,j)/b(i,j)
	end do
      end do
      do j=1,m
	e(n-2,j) = (a(n-1,j)-e(n-3,j)*c(n-3,j))/b(n-2,j)
      end do
c
c     compute  inner product (e,d) for each j in sum(j)
c
      do j=1,m
	sum(j) = 0.
      end do
      do i=1,n-2
	do j=1,m
	  sum(j) = sum(j)+e(i,j)*d(i,j)
	end do
      end do
c
c     set last diagonal element
c
      do j=1,m
	b(n-1,j) = b(n-1,j)-sum(j)
      end do
      return
      end

      subroutine transp(n,amat)
c
c     transpose n by n real matrix
c
      implicit none
      integer n,i,j
      real amat(n,n),temp
      do i=1,n-1
	do j=i+1,n
	  temp = amat(i,j)
	  amat(i,j) = amat(j,i)
	  amat(j,i) = temp
	end do
      end do
      return
      end

      subroutine sgfa (a,lda,n,ipvt,info)
      integer lda,n,ipvt(1),info                                                
      real a(lda,1)                                                             
      real t                                                                    
      integer isfmax,j,k,kp1,l,nm1
      info = 0                                                                  
      nm1 = n - 1                                                               
      if (nm1 .lt. 1) go to 70                                                  
      do 60 k = 1, nm1                                                          
         kp1 = k + 1                                                            
	 l = isfmax(n-k+1,a(k,k),1) + k - 1
         ipvt(k) = l                                                            
         if (a(l,k) .eq. 0.0e0) go to 40                                        
            if (l .eq. k) go to 10                                              
               t = a(l,k)                                                       
               a(l,k) = a(k,k)                                                  
               a(k,k) = t                                                       
   10       continue                                                            
            t = -1.0e0/a(k,k)                                                   
	    call sscl(n-k,t,a(k+1,k),1)
            do 30 j = kp1, n                                                    
               t = a(l,j)                                                       
               if (l .eq. k) go to 20                                           
                  a(l,j) = a(k,j)                                               
                  a(k,j) = t                                                    
   20          continue                                                         
	       call sxpy(n-k,t,a(k+1,k),1,a(k+1,j),1)
   30       continue                                                            
         go to 50                                                               
   40    continue                                                               
            info = k                                                            
   50    continue                                                               
   60 continue                                                                  
   70 continue                                                                  
      ipvt(n) = n                                                               
      if (a(n,n) .eq. 0.0e0) info = n                                           
      return                                                                    
      end                                                                       
                                                                                
      subroutine sgsl (a,lda,n,ipvt,b,job)
      integer lda,n,ipvt(1),job                                                 
      real a(lda,1),b(1)                                                        
      real sdt,t
      integer k,kb,l,nm1                                                        
      nm1 = n - 1                                                               
      if (job .ne. 0) go to 50                                                  
         if (nm1 .lt. 1) go to 30                                               
         do 20 k = 1, nm1                                                       
            l = ipvt(k)                                                         
            t = b(l)                                                            
            if (l .eq. k) go to 10                                              
               b(l) = b(k)                                                      
               b(k) = t                                                         
   10       continue                                                            
	    call sxpy(n-k,t,a(k+1,k),1,b(k+1),1)
   20    continue                                                               
   30    continue                                                               
         do 40 kb = 1, n                                                        
            k = n + 1 - kb                                                      
            b(k) = b(k)/a(k,k)                                                  
            t = -b(k)                                                           
	    call sxpy(k-1,t,a(1,k),1,b(1),1)
   40    continue                                                               
      go to 100                                                                 
   50 continue                                                                  
         do 60 k = 1, n                                                         
	    t = sdt(k-1,a(1,k),1,b(1),1)
            b(k) = (b(k) - t)/a(k,k)                                            
   60    continue                                                               
         if (nm1 .lt. 1) go to 90                                               
         do 80 kb = 1, nm1                                                      
            k = n - kb                                                          
	    b(k) = b(k) + sdt(n-k,a(k+1,k),1,b(k+1),1)
            l = ipvt(k)                                                         
            if (l .eq. k) go to 70                                              
               t = b(l)                                                         
               b(l) = b(k)                                                      
               b(k) = t                                                         
   70       continue                                                            
   80    continue                                                               
   90    continue                                                               
  100 continue                                                                  
      return                                                                    
      end                                                                       
                                                                                
      real function sdt(n,sx,incx,sy,incy)
      real sx(1),sy(1),stemp                                                    
      integer i,incx,incy,ix,iy,m,mp1,n                                         
      stemp = 0.0e0                                                             
      sdt = 0.0e0
      if(n.le.0)return                                                          
      if(incx.eq.1.and.incy.eq.1)go to 20                                       
      ix = 1                                                                    
      iy = 1                                                                    
      if(incx.lt.0)ix = (-n+1)*incx + 1                                         
      if(incy.lt.0)iy = (-n+1)*incy + 1                                         
      do 10 i = 1,n                                                             
        stemp = stemp + sx(ix)*sy(iy)                                           
        ix = ix + incx                                                          
        iy = iy + incy                                                          
   10 continue                                                                  
      sdt = stemp
      return                                                                    
   20 m = mod(n,5)                                                              
      if( m .eq. 0 ) go to 40                                                   
      do 30 i = 1,m                                                             
        stemp = stemp + sx(i)*sy(i)                                             
   30 continue                                                                  
      if( n .lt. 5 ) go to 60                                                   
   40 mp1 = m + 1                                                               
      do 50 i = mp1,n,5                                                         
        stemp = stemp + sx(i)*sy(i) + sx(i + 1)*sy(i + 1) +                     
     +   sx(i + 2)*sy(i + 2) + sx(i + 3)*sy(i + 3) + sx(i + 4)*sy(i + 4)        
   50 continue                                                                  
   60 sdt = stemp
      return                                                                    
      end                                                                       
                                                                                
      integer function isfmax(n,sx,incx)
      real sx(1),smax                                                           
      integer i,incx,ix,n                                                       
      isfmax = 0
      if( n .lt. 1 ) return                                                     
      isfmax = 1
      if(n.eq.1)return                                                          
      if(incx.eq.1)go to 20                                                     
      ix = 1                                                                    
      smax = abs(sx(1))                                                         
      ix = ix + incx                                                            
      do 10 i = 2,n                                                             
         if(abs(sx(ix)).le.smax) go to 5                                        
	 isfmax = i
         smax = abs(sx(ix))                                                     
    5    ix = ix + incx                                                         
   10 continue                                                                  
      return                                                                    
   20 smax = abs(sx(1))                                                         
      do 30 i = 2,n                                                             
         if(abs(sx(i)).le.smax) go to 30                                        
	 isfmax = i
         smax = abs(sx(i))                                                      
   30 continue                                                                  
      return                                                                    
      end                                                                       

      subroutine sxpy(n,sa,sx,incx,sy,incy)
      real sx(1),sy(1),sa                                                       
      integer i,incx,incy,ix,iy,m,mp1,n                                         
      if(n.le.0)return                                                          
      if (sa .eq. 0.0) return                                                   
      if(incx.eq.1.and.incy.eq.1)go to 20                                       
      ix = 1                                                                    
      iy = 1                                                                    
      if(incx.lt.0)ix = (-n+1)*incx + 1                                         
      if(incy.lt.0)iy = (-n+1)*incy + 1                                         
      do 10 i = 1,n                                                             
        sy(iy) = sy(iy) + sa*sx(ix)                                             
        ix = ix + incx                                                          
        iy = iy + incy                                                          
   10 continue                                                                  
      return                                                                    
   20 m = mod(n,4)                                                              
      if( m .eq. 0 ) go to 40                                                   
      do 30 i = 1,m                                                             
        sy(i) = sy(i) + sa*sx(i)                                                
   30 continue                                                                  
      if( n .lt. 4 ) return                                                     
   40 mp1 = m + 1                                                               
      do 50 i = mp1,n,4                                                         
        sy(i) = sy(i) + sa*sx(i)                                                
        sy(i + 1) = sy(i + 1) + sa*sx(i + 1)                                    
        sy(i + 2) = sy(i + 2) + sa*sx(i + 2)                                    
        sy(i + 3) = sy(i + 3) + sa*sx(i + 3)                                    
   50 continue                                                                  
      return                                                                    
      end                                                                       

      subroutine sscl(n,sa,sx,incx)
      real sa,sx(1)                                                             
      integer i,incx,m,mp1,n,nincx                                              
      if(n.le.0)return                                                          
      if(incx.eq.1)go to 20                                                     
      nincx = n*incx                                                            
      do 10 i = 1,nincx,incx                                                    
        sx(i) = sa*sx(i)                                                        
   10 continue                                                                  
      return                                                                    
   20 m = mod(n,5)                                                              
      if( m .eq. 0 ) go to 40                                                   
      do 30 i = 1,m                                                             
        sx(i) = sa*sx(i)                                                        
   30 continue                                                                  
      if( n .lt. 5 ) return                                                     
   40 mp1 = m + 1                                                               
      do 50 i = mp1,n,5                                                         
        sx(i) = sa*sx(i)                                                        
        sx(i + 1) = sa*sx(i + 1)                                                
        sx(i + 2) = sa*sx(i + 2)                                                
        sx(i + 3) = sa*sx(i + 3)                                                
        sx(i + 4) = sa*sx(i + 4)                                                
   50 continue                                                                  
      return                                                                    
      end                                                                       

      subroutine rcd2spp(nx,ny,phi,rhs,cofx,cofy)
c
c     gauss-seidel red/black point relaxation
c
      implicit none
      integer nx,ny,i,j,ist,ifn,jst,jfn
      integer intl,nxa,nxb,nyc,nyd,ixp,jyq,iex,jey,nfx,nfy,iguess,
     +             maxcy,method,nwork,lwork,itero,ngrid,klevel,kcur,
     +             kcycle,iprer,ipost,intpol,kps
      common/icud2/intl,nxa,nxb,nyc,nyd,ixp,jyq,iex,jey,nfx,nfy,iguess,
     +             maxcy,method,nwork,lwork,itero,ngrid,klevel,kcur,
     +             kcycle,iprer,ipost,intpol,kps
      complex phi(0:nx+1,0:ny+1),rhs(nx,ny),cofx(nx,3),cofy(ny,3)
c
c     set loop limits to avoid specified boundaries
c     in red/black sweeps
c
      ist = 1
      if (nxa.eq.1) ist = 3
      ifn = nx
      if (nxb.eq.1) ifn = nx-1
      jst = 1
      if (nyc.eq.1) jst = 3
      jfn = ny
      if (nyd.eq.1) jfn = ny-1
c
c    periodic adjustment bypass block
c
      if (nxa*nyc.ne.0) then
c
c     relax on red grid points
c
C$OMP PARALLEL DO SHARED(cofx,cofy,rhs,phi,ist,ifn,jst,jfn) PRIVATE(i,j)
	do i=ist,ifn,2
	  do j=jst,jfn,2
	    phi(i,j) = (rhs(i,j) - (
     +                  cofx(i,1)*phi(i-1,j)+cofx(i,2)*phi(i+1,j) +
     +                  cofy(j,1)*phi(i,j-1)+cofy(j,2)*phi(i,j+1)))/
     +                 (cofx(i,3)+cofy(j,3))
	  end do
	end do
C$OMP PARALLEL DO SHARED(cofx,cofy,rhs,phi,ist,ifn,jst,jfn) PRIVATE(i,j)
	do i=2,ifn,2
	  do j=2,jfn,2
	    phi(i,j) = (rhs(i,j) -
     +                 (cofx(i,1)*phi(i-1,j)+cofx(i,2)*phi(i+1,j) +
     +                 cofy(j,1)*phi(i,j-1)+cofy(j,2)*phi(i,j+1)))/
     +                 (cofx(i,3)+cofy(j,3))
	  end do
	end do
c
c     relax on black grid points
c
C$OMP PARALLEL DO SHARED(cofx,cofy,rhs,phi,ist,ifn,jst,jfn) PRIVATE(i,j)
	do i=ist,ifn,2
	  do j=2,jfn,2
	    phi(i,j) = (rhs(i,j) -
     +                 (cofx(i,1)*phi(i-1,j)+cofx(i,2)*phi(i+1,j) +
     +                 cofy(j,1)*phi(i,j-1)+cofy(j,2)*phi(i,j+1)))/
     +                 (cofx(i,3)+cofy(j,3))
	  end do
	end do
C$OMP PARALLEL DO SHARED(cofx,cofy,rhs,phi,ist,ifn,jst,jfn) PRIVATE(i,j)
	do i=2,ifn,2
	  do j=jst,jfn,2
	    phi(i,j) = (rhs(i,j) -
     +                 (cofx(i,1)*phi(i-1,j)+cofx(i,2)*phi(i+1,j) +
     +                 cofy(j,1)*phi(i,j-1)+cofy(j,2)*phi(i,j+1)))/
     +                 (cofx(i,3)+cofy(j,3))
	  end do
	end do
	return
      end if
c
c    set periodic virtual boundaries
c
      if (nxa.eq.0) then
	do j=1,ny
	  phi(0,j) = phi(nx-1,j)
	  phi(nx+1,j) = phi(2,j)
	end do
      end if
      if (nyc.eq.0) then
	do i=1,nx
	  phi(i,0) = phi(i,ny-1)
	  phi(i,ny+1) = phi(i,2)
	end do
      end if
c
c     relax on red grid points
c
C$OMP PARALLEL DO SHARED(cofx,cofy,rhs,phi,ist,ifn,jst,jfn) PRIVATE(i,j)
      do i=ist,ifn,2
	do j=jst,jfn,2
	    phi(i,j) = (rhs(i,j) -
     +                 (cofx(i,1)*phi(i-1,j)+cofx(i,2)*phi(i+1,j) +
     +                 cofy(j,1)*phi(i,j-1)+cofy(j,2)*phi(i,j+1)))/
     +                 (cofx(i,3)+cofy(j,3))
	end do
      end do
C$OMP PARALLEL DO SHARED(cofx,cofy,rhs,phi,ist,ifn,jst,jfn) PRIVATE(i,j)
      do i=2,ifn,2
	do j=2,jfn,2
	    phi(i,j) = (rhs(i,j) -
     +                 (cofx(i,1)*phi(i-1,j)+cofx(i,2)*phi(i+1,j) +
     +                 cofy(j,1)*phi(i,j-1)+cofy(j,2)*phi(i,j+1)))/
     +                 (cofx(i,3)+cofy(j,3))
	end do
      end do
c
c    ensure periodic virtual boundary red points are set
c
      if (nxa.eq.0) then
	do j=1,ny
	  phi(0,j) = phi(nx-1,j)
	  phi(nx+1,j) = phi(2,j)
	end do
      end if
      if (nyc.eq.0) then
	do i=1,nx
	  phi(i,0) = phi(i,ny-1)
	  phi(i,ny+1) = phi(i,2)
	end do
      end if
c
c     relax on black grid points
c
C$OMP PARALLEL DO SHARED(cofx,cofy,rhs,phi,ist,ifn,jst,jfn) PRIVATE(i,j)
      do i=ist,ifn,2
	do j=2,jfn,2
	    phi(i,j) = (rhs(i,j) -
     +                 (cofx(i,1)*phi(i-1,j)+cofx(i,2)*phi(i+1,j) +
     +                 cofy(j,1)*phi(i,j-1)+cofy(j,2)*phi(i,j+1)))/
     +                 (cofx(i,3)+cofy(j,3))
	end do
      end do
C$OMP PARALLEL DO SHARED(cofx,cofy,rhs,phi,ist,ifn,jst,jfn) PRIVATE(i,j)
      do i=2,ifn,2
	do j=jst,jfn,2
	    phi(i,j) = (rhs(i,j) -
     +                 (cofx(i,1)*phi(i-1,j)+cofx(i,2)*phi(i+1,j) +
     +                 cofy(j,1)*phi(i,j-1)+cofy(j,2)*phi(i,j+1)))/
     +                 (cofx(i,3)+cofy(j,3))
	end do
      end do
c
c     final set of periodic virtual boundaries
c
      if (nxa.eq.0) then
	do j=1,ny
	  phi(0,j) = phi(nx-1,j)
	  phi(nx+1,j) = phi(2,j)
	end do
      end if
      if (nyc.eq.0) then
	do i=1,nx
	  phi(i,0) = phi(i,ny-1)
	  phi(i,ny+1) = phi(i,2)
	end do
      end if
      return
      end

