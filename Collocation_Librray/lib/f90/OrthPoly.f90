Module Orth_Poly  ! fundamental calculations with Jacobi polynomials
   Use Array_Print, only : ArrayPrint, VectorPrint
   include 'defs.fi'
   Private
   ! Public functions (many should be made private after testing ---------------
   Public Jacobi_Pm           ! calculate polynomials & derivatives
   Public Jacobi_Rec          ! calculate range of polynomials with recurrence formula
   Public Jacobi_All          ! calculate all polynomials from 0 to n by recurrence
   Public Jacobi_Quadrature   ! nodes & weights
   Public Jacobi_Quad         ! nodes & weights
   Public Jacobi_Weights      ! new quadrature and barycentric weight calculation
   Public Jacobi_Xroot        ! root calculation by Newton-Raphson
   Public Jacobi_Iter         ! calculate root correction
   Public Jacobi_Pn_update    ! applies Taylor series to P & P' (not currently used)
   Public RootEstimate0       ! simple root estimation
   Public RootEstimate        ! best root estimates
   Public Jacobi_Series       ! value of a Jacobi Polynomial Series
   Public Legendre            ! value of Legendre Polynomials & derivatives
   Public Legendre_Deriv      ! Legendre Polynomials derivative coefficients, modal
   Public Jacobi_Deriv1       ! (1-x^2)Pn' = [c(1)*x + c(2)]Pn + c(3)*Pn-1, calculate c(1:3)
   Public Jacobi_Deriv2       ! (1-x^2)Pn" = [c(1)*x + c(2)]Pn' + c(3)*Pn, calculate c(1:3)
   Public Jacobi_Deriv        ! dPn/dx = sum(d(n-1,0:n-1)*Pn(0:n-1)), calculate d(0:n-1,0:n-1)
   Public Jacobi_Ends         ! endpoint values
   Public Jacobi_Lead         ! lead coefficient
   Public Legendre_Transform   ! discrete Legendre Transform (not valid for Gauss pts)
   Public Jacobi_Transform     ! creates Legendre Transform via Jacobi Transform
   Public Polyset               ! sets a & b given type, geometry, symmetry (shortcut flag)
   Public x_Jacobi            ! recurrence coefficients, leading coef., norms
   Public d_Jacobi            ! integrals of Pn^2 (this is redundant due to x_Jacobi)
   Public FirstLT
   ! Public module parameters -----------------------------
   ! most of these should be made private, removed, or made parameters after testing
   ! >> Warning >> if Nboundary .ne. Nasymptotic, Jacobi_Pcomposite must be modified << Warning <<
   Integer, Parameter, Public :: Nboundary = 40          ! recurrent & boundary asymptotic (no interior)
   Integer, Parameter, Public :: Nasymptotic = Nboundary ! see Warning! use interior, recurrent & boundary
   logical, parameter         :: Default_monic = .false. ! default for optional monic arguments
   logical, parameter         :: Default_shift = .false. ! default for optional shift arguments
   logical, parameter         :: Always_short = .true.   ! use shortcut for roots & weights for ultraspherical
   Logical, Parameter         :: ShortOK = .true.        ! true if ok to use shortcut calcs
   Logical, Parameter         :: ShiftOK = .true.        ! true to use shifted shortcut calcs
   integer, private           :: Ultra_method = 1        ! controls version of Jacobi_AsymUltra
   integer, parameter         :: Deriv_meth = 0          ! controls how P' calculated, interior asymptotic
   real(float), private       :: ThetaFac                ! factor used in Jacobi_Pcalc
   integer, private           :: JacobiType              ! polynomial type 1-6,11,12 (see below)
   Integer, Parameter, Public :: Gauss=1, Lobatto=2, Chebyshev=3, RadauR=4, RadauL=5, Chebyshev1=6
   Integer, Parameter, Public :: ShortGauss=11, ShortLobatto=12
   Integer, parameter         :: Short=-1, Nonsymmetric=0, Planar=1, Cylindrical=2, Spherical=3
   Real(floatq),Parameter     :: qpi = 3.14159265358979323846264338327950288419716939937510_floatq
   Real(float), Parameter     :: pi = real(qpi,float)
   Real(float), Parameter     :: pi2 = (0.5_float)*pi
   Real(float), parameter     :: fdegree = pi/(180.0_float)
   Real(float), private       :: Debug = 0               ! controls extra output
   ! Formats ----------------------------------------------
   Character (len=*), Parameter :: fmt16 = '(25(f21.16,a1))'
   Character (len=*), Parameter :: fmt17 = '(25(f20.17,a1))'
   Character (len=*), Parameter :: fmt22 = '(25(es24.17,a1))'
   Character (len=*), Parameter :: fmtxx = '(25(es9.2,a1))'
   Character (len=*), Parameter :: fmt17i = '(i15,a1,25(i20,a1))'
   !-------------------------------------------------------------------------
   Contains
   !----------------------------------------------------------------------------
   Function Jacobi_Quadrature(n,ab,geom,wbfac,Nneg,MaxNR) result(w)
      integer, intent(in) :: n,geom
      Real(float), intent(in) :: ab(2)
      real(float), intent(out) :: wbfac
      integer, optional :: Nneg,MaxNR
      real(float) :: w(max(0,min(1,geom)):n+1,0:3)
      integer :: nx,ineg,g
      ! -----------------------------------------
      ! Calculate roots & weights numbered 0 to n+1 or 1 to n+1 (symmetric)
      ! W(:,0) - x roots
      ! W(:,1) - theta roots
      ! W(:,2) - barycentric weights
      ! W(:,3) - quadrature weights
      ! n - interior points (+2 for nonsymmetric, +1 for symmetric)
      ! a,b - polynomial alpha & beta values
      ! geom =-1 nonsymmetric ultraspherical using shortcut approach
      !      = 0 nonsymmetric
      !      = 1,2,3 planar,cylindrical & spherical symmetric
      ! wbfac - scaling factor for barycentric weights
      ! Nneg - number of negative x values (optional)
      !        correct values returned if not argument supplied
      !        absolute values returned
      !        theta = cos(x) is always returned, exceptions??
      ! MaxNR - maximum iterations (optional)
      ! -----------------------------------------
      nx = n;  ineg = 0;  g = geom
      g = max(0,g)
      if(ab(1) == ab(2))then
         nx = n/2          ! ultraspherical (cylindrical Gauss too)
         if(g == 0 .and. Always_short)g = -1  ! uses shortcut

      endif
      call Jacobi_Quad(n,nx,g,ab,w,wbfac,ineg,MaxNR)
      if(present(Nneg))then
         Nneg = ineg
      elseif(g <= Nonsymmetric)then
         w(:ineg-1,0) = -w(:ineg-1,0)
         w(:ineg-1,1) = pi - w(:ineg-1,1)
      else  ! symmetric, convert angles
         w(:ineg,1) = pi - w(:ineg,1)
      endif
   End Function Jacobi_Quadrature
   !----------------------------------------------------------------------------
   Subroutine Jacobi_Quad(nn,nx,g,ab,w,wbfac,nneg,MaxNR)
      integer, intent(in) :: nn,nx,g
      Real(float), intent(in) :: ab(2)
      real(float), intent(out) :: w(max(0,min(1,g)):,:),wbfac
      integer, intent(out) :: nneg
      integer, optional :: MaxNR
      real(float) :: a,b,xr(nx,0:1),pn(nx,0:1)
      Integer :: n
      logical theta
      ! -----------------------------------------
      ! Calculate roots and quadrature & barycentric weights
      !  nn - interior points, polynomial order (n = nn/2 for shortcut)
      !  nx - number unique root values
      !  g - geometry -1 to 3
      !  a,bb - Jacobi alpha & beta weight parameters (b = +/-0.5 for shortcut)
      !  Jacobi_Xroot - calulates roots & pn
      !     xr(:,0) - final accurate roots
      !     xr(:,1) - values at beginning of last iteration
      !     pn(:,0:1) - P and P' calculated at xr(:,1)
      !     theta = .true. when theta coordinate, i.e xr is ange, P' is w.r.t theta
      !  Jacobi_Weights - calculates weights
      !     w(:0:3) - x, theta, w - barycentric, w - quadrature
      !     wbfac - scaling factor for barycentric weights (returned)
      !  nneg - abs(x) returned, nneg is number that are negative
      !  MaxNR - optional, number of iterations performed
      ! -----------------------------------------
      n = nn;   a = ab(1);   b = ab(2)
      if(g < 0)then  ! shortcut ultraspherical
         b = (2*mod(n,2)-1)*0.5_float  ! -0.5 even, +0.5 odd
         n = n/2
      endif
      if(Debug > 0)write(iout,'("# Jacobi_Quad : n,nx=",2i7,a,3f5.1,a,i2)') &
         nn,nx,', ab,b =',ab,b,', geom =',g
      call Jacobi_Xroot(n,a,b,xr,pn,theta,nneg,MaxNR) ! calculate roots
      w = Jacobi_Weights(nn,g,a,b,xr,pn,wbfac,theta,nneg)
   End Subroutine Jacobi_Quad
   !----------------------------------------------------------------------------
   Subroutine Jacobi_Xroot(n,a,b,x,p,theta,nneg,MaxNR)
      integer, intent(in) :: n
      Real(float), intent(in) :: a,b
      Real(float), intent(out) :: x(:,0:),p(:,0:)
      logical, intent(out) :: theta
      integer, intent(out) :: nneg
      integer, optional :: MaxNR
      Integer :: niter,nb,nn
      ! ------------------------------
      ! Calculates Jacobi roots using higher order iteration.
      ! Uses higher order methods & vector code, i.e. solves for all roots at once
      !  n - order of polynomial
      !  a,b - Jacobi alpha & beta parameters
      !  x(:,0:1) - 1 accurate root, 0 next to last iteration
      !  p - polynomial & derivative values evaluated at x(:,0)
      !  theta = .true. if output values are in theta coordinates
      !  MaxNR - max iterations (optonal)
      !  nneg - number of values that are negative
      ! ------------------------------
      niter = 1;   if(epsilon(1.0_float) < 1.0e-20_float)niter = 2
      if(Present(MaxNR))niter = MaxNR
      nn = n;  nb = nint(b*2.0_float)
      if(abs(nb)==1)nn = 2*n + (nb+1)/2
      theta = nn > min(Nboundary,Nasymptotic)      ! calcs in theta if asymptotic used
      x(:,0) = RootEstimate(n,size(x,1),a,b,nneg)  ! for shortcut b=+/-.5 & n = nn/2
      if(Debug > 0)then
         write(iout,'("# Jacobi_Xroot: n =",i5,a,2f4.1,2(a,i2),a,L2)') &
            n,', a,b =',a,b,', niter =',niter,', nneg =',nneg,', Theta =',theta
         if(Debug > 2)call arrayprint('# Jacobi_Xroot: estimates',x(:,0),cos(x(:,0)))
      endif
      if(theta)then
         call Jacobi_Xroot_thet(n,a,b,x,p,niter,nneg)
      else
         call Jacobi_Xroot_X(n,a,b,x,p,niter,nneg)
      endif
      if(Debug > 2)call arrayprint('Jacobi_Xroot final x0,x1,p,p',x,p,fmtx=fmt22)
   End Subroutine Jacobi_Xroot
   !----------------------------------------------------------------------------
   Subroutine Jacobi_Xroot_thet(n,a,b,x,p,niter,nneg)
      integer, intent(in) :: n,niter
      Real(float), intent(in) :: a,b
      Real(float), intent(inout) :: x(:,0:)
      Real(float), intent(out) :: p(:,0:)
      integer, intent(in) :: nneg
      Real(float) :: dp(size(x,1))
      Integer :: i
      ! -------------------------------------------------
      ! Solve for root using theta angles, Newton-Raphson
      !  x(:,1) is root (theta)
      !  x(:,0) is root, previous to last iteration
      !  p(:,0:1) - Pn & Pn' evaluated at x(:,0)
      ! -------------------------------------------------
      if(Debug > 1)write(iout,'(a,10i5)')'# Jacobi_Xroot_thet: n,nneg =',n,nneg
      do i=1,niter               ! iteration, with debug & tolerance check
         p(:,0:1) = Jacobi_Pm(x(:,0),n,1,a,b,Theta=.true.,nneg=nneg)
         dp = p(:,0)/p(:,1)
         x(1:nneg,1) = x(1:nneg,0) + dp(1:nneg)
         x(nneg+1:,1) = x(nneg+1:,0) - dp(nneg+1:)
         if(Debug > 1)call Jacobi_Xroot_debug(i,x,dp)
         if(i < niter)x(:,0) = x(:,1)
      enddo
   End Subroutine Jacobi_Xroot_thet
   !----------------------------------------------------------------------------
   Subroutine Jacobi_Xroot_X(n,a,b,x,p,niter,kneg)
      integer, intent(in) :: n,niter,kneg
      Real(float), intent(in) :: a,b
      real(float), intent(inout) :: x(:,0:)
      real(float), intent(out) :: p(:,0:)
      real(float) :: dp(size(x,1))
      integer :: ord,i,nn,ib
      integer, parameter :: iord(3) = (/35,9,0/)
      ! -------------------------------------------------
      ! Solve for root using x coordinates, higher order iteration
      !  n is polynomial degree, nn is corresponding n for shortcut cases
      !  x(:,1) is current most accurate root
      !  x(:,0) is root at beginning of last iteration
      !  p(:,0:1) - Pn & Pn' evaluated at x(:,0)
      !  niter - max no. of iterations
      !  kneg - no. of negative values
      !  ord - 0 to 5 order of higher order, 0 is Newton-Raphson
      ! there is no real need to account for negatives for small n, but it was already coded
      ! -------------------------------------------------
      ib = nint(b*2.0_float)  ! -1 & 1 or shortcut, 0 & 2 are normal cases
      nn = n + mod(ib+2,2)*(n + (ib+1)/2) ! n, 2n, 2n+1 for b = 0,1,-/+0.5
      ord = IFirstLT(iord,nn)
      if(Debug > 1)write(iout,'(a,10i5)')'# Jacobi_Xroot_X: n,nn,ord,nneg =',n,nn,ord,kneg
      x(:,0) = cos(x(:,0))
      do i=1,niter               ! iteration, with debug
         dp = Jacobi_Iter(n,x(:,0),a,b,p,Ord,kneg)
         x(1:kneg,1)  = x(1:kneg,0)  + dp(1:kneg)
         x(kneg+1:,1) = x(kneg+1:,0) - dp(kneg+1:)
         if(Debug > 1)call Jacobi_Xroot_debug(i,x,dp)
         if(i < niter)x(:,0) = x(:,1)
      enddo
   End Subroutine Jacobi_Xroot_X
   !----------------------------------------------------------------------------
   Function Jacobi_Iter(n,x,a,b,p,ord,kneg)  result(dp)
      Real(float) :: x(:),a,b,p(:,0:),dp(size(x))
      integer :: n,ord,k,kneg
      Real(float) :: dpx(size(x)),dx(size(x),2:ord),r2(size(x)),xn(size(x))
      Real(float) :: c2(3),d3(0:2),d4(0:3)
      Real(float), parameter :: fac(2:5) = (1.0_float)/Real((/2,6,24,120/),float)
      ! -------------------------------
      ! This implements the higher order iteration scheme
      !  n - polynomial no.
      !  x - current root values
      !  a,b - Jacobi parameters, alpha and beta
      !  p(:,0:1) - P and P'
      !  ord = 1 to 5 (1 is N-R)
      !  kneg - number of negative values
      !  dp - calculated change in x
      ! -------------------------------
      p(:,0:1) = Jacobi_Pm(x,n,1,a,b,Theta=.false.,nneg=kneg)
      dp = p(:,0)/p(:,1)
      c2 = Jacobi_Deriv2(n,a,b)
      if(Ord > 1)then   ! first correction (3rd order) --------------
         r2 = 1.0_float - x*x
         dpx = dp/r2
         xn(:kneg) = -x(:kneg)   ! include negatives for higher order
         xn(kneg+1:) = x(kneg+1:)
         dx(:,2) = c2(1)*xn + c2(2)
      endif
      if(Ord > 2)then   ! 3rd term (4th order) ----------------------
         d3(0) = (c2(3)+c2(2)*c2(2))*2.0_float - c2(1)
         d3(1) = c2(2)*(c2(1)-0.5_float)*4.0_float
         d3(2) = (c2(1)*c2(1)- c2(1)*0.5_float - c2(3))*(2.0_float)
         dx(:,3) = d3(0) + xn*(d3(1) + xn*d3(2))
        !dx(:,3) = ((dx(:,2)*(dx(:,2)-x) + (c2(3)-c2(1)*0.5_float)*r2)*2.0_float)*fac(3)
      endif
      if(Ord > 3)then   ! 4th term (5th order) ----------------------
         d4(0) = (3.0_float)*c2(2)*(c2(3)*2.0_float + d3(0)) - d3(1)
         d4(1) = (3.0_float)*(c2(1)*c2(3)*2.0_float + c2(2)*d3(1)) &
                  + (c2(1)*3.0_float-4.0_float)*d3(0) - (2.0_float)*d3(2)
         d4(2) = (3.0_float)*c2(2)*(d3(2) - c2(3)*2.0_float) &
                  + (c2(1)*3.0_float-4.0_float)*d3(1) + d3(1)
         d4(3) = (-6.0_float)*c2(1)*c2(3) + (c2(1)*3.0_float-4.0_float)*d3(2) &
                  + (2.0_float)*d3(2)
         dx(:,4) = d4(0) + xn*(d4(1) + xn*(d4(2) + xn*d4(3)))
        !dx(:,4) = (r2*(c2(3)*6.0_float*dx(:,2)-d3(1)-d3(2)*x*2.0_float) &
        !          + ((c2(1)*x + c2(2))*3.0_float-x*4.0_float)*dx(:,3))*fac(4)
      endif
      if(Ord > 4)then   ! 5th term (6th order) ----------------------
         dx(:,5) = (12.0_float*c2(3)*dx(:,3) - (d4(1)+xn*(2.0_float*d4(2)+xn*3.0_float*d4(3))))*r2 &
               + (4.0_float*(c2(1)*xn+c2(2))-6.0_float*xn)*dx(:,4)
      endif
      do k = 2,Ord   ! calculate correction --------
         dx(:,k) = dx(:,k)*fac(k)
      enddo
      if(Debug > 2 .and. Ord > 1)then
         call arrayprint('Jacobi_Iter:',xn,dp,dpx,dx,fmtx=fmt22)
      endif
      Select Case(Ord)
      case(1)     ! Newton-Raphson
         continue
      case(2)     ! 3rd order
         dp = dp*(1.0_float + dpx*dx(:,2))
      case(3)     ! 4th order
         dp = dp*(1.0_float + dpx*(dx(:,2) + dpx*dx(:,3)))
      case(4)     ! 5th order
         dp = dp*(1.0_float + dpx*(dx(:,2) + dpx*(dx(:,3) + dpx*dx(:,4))))
      case default   ! 6th order
         dp = dp*(1.0_float + dpx*(dx(:,2) + dpx*(dx(:,3) + dpx*(dx(:,4) + dpx*dx(:,5)))))
      End Select
   End Function Jacobi_Iter
   !----------------------------------------------------------------------------
   subroutine Jacobi_Xroot_debug(i,x,dp)
      integer :: i,idp,nx
      real(float) :: dpmx,x(:,:),dp(:)
      ! debug output, root calculations
      nx = size(x,1)
      idp = maxloc(abs(dp),1);
      dpmx = abs(dp(idp))
      write(iout,'("# Jacobi_Xroot:",i2,i6,21es10.2)')i,idp,dpmx,dp(max(1,nx-19):nx)
      if(Debug > 3)call arrayprint('# x1,x0,dp:',x,dp,fmtx=fmt22)
   end subroutine Jacobi_Xroot_debug
   !----------------------------------------------------------------------------
   Function Jacobi_Weights(n,g,a,b,x,pn,wbfac,theta,nneg) result(w)
      integer, intent(in) :: n,g
      real(float), intent(in) :: a,b,x(:,0:),pn(:,0:)
      real(float), intent(out) :: wbfac
      logical, intent(in) :: Theta
      integer, intent(inout) :: nneg
      real(float) :: w(max(0,min(1,g)):n+1,0:3)
      real(float) :: pd(size(x,1)),xx(size(x,1)),q(size(x,1))
      real(float) :: c(3),cw(0:3,0:1),rx(size(x,1),0:1),dw(size(x,1),0:1)
      integer :: nn,nx,nodd,ns,ib,ord,i0,j,sym,iw(3),ex(0:1,0:1)
      real(float), parameter :: x0(0:2,2) = &
         reshape((/1.0_float,0.0_float,1.0_float,0.0_float,pi2,0.0_float/),(/3,2/))
      integer, parameter :: iord(4) = (/9999,50,7,0/)
      ! ----------------------------------------------------
      ! Calculates barycentric and quadrature weights using higher order method
      !  n - degree of polynomial (number interior points)
      !  g - geometry <= 0 or 1,2,3
      !  a,b - alpha,beta
      !  x(:,1) - accurate roots, x(:,0) less accurate
      !  pn(:,0:1) - p & p' at x(:,0)
      !  wbfac - scaling factor for barycentric weights
      !  theta = .true. x & pn' are in theta coordinates, .false. if not
      !  nneg - number of negative x values
      !  nn - corresponding full polynomial (for symmetric cases)
      !  ns - corresponding short cut polynomial (for short cut nonsymmetric cases)
      !  nx - number unique x points
      !
      ! returns W(sym:n+1,0:3), sym = 0/1
      !  0 - roots as x numbered 0/1 to n+1
      !  1 - roots as arccos(x) (theta) (not always when nneg)
      !  2 - barycentric weights
      !  3 - quadrature weights
      ! ----------------------------------------------------
      sym = max(0,min(1,g))
      nx = size(x,1);   nodd = mod(n,2)
      ib = nint(b*2.0_float)+1   ! 0,1,2,3 for -0.5,0,+0.5,1
      if(mod(ib,2) == 1)then     ! b = 0 or 1
         nn = n
         ns = n
      else                       ! b = +/-0.5
         nn = n + sym*(n + ib/2) ! n, or 2n or 2n+1 (symmetric problems only)
         ns = n/(2-sym)          ! n/2 for shortcut cases (not symmetric)
      endif
      ord = IFirstLT(iord,nn) ! order of weight calculation
      if(Debug > 0)then
         write(iout,'(a,5i5,2f5.1,5i3)')'Jacobi_Weights: nn,n,nx,ns,nneg,a,b,ord,s', &
            nn,n,nx,ns,nneg,a,b,ord,sym
      endif
      ! get Cn, w_0, w_ctr, w_n+1, r(x) exponents, Wb factor
      call jac_weight_values(n,sym,a,b,ex,cw,wbfac)
      w = 0.0_float
      i0 = n+1-nx             ! 1st unique: 1 or 1st to right of midpoint (ultraspherical)
      iw = (/sym,nx+1,n+1/)   ! 1st, center, last (some overwrites below
      w(iw,0:1) = x0          ! left,center,right x & theta
      w(iw,2:3) = cw(1:3,0:1) ! left,center,right wb & wq
      if(theta)then
         pd = -pn(:,1)/sin(x(:,0))
         q = -pn(:,0)/(pn(:,1)*xx)  ! P/((1-x^2)P')
         rx(:,0) = jac_wts_rxt(ex(:,0),x(:,0),nneg)
         rx(:,1) = jac_wts_rxt(ex(:,1),x(:,0),nneg)
         xx = cos(x(:,0))
         xx(:nneg) = -xx(:nneg)
         w(i0:n,0) = cos(x(:,1)) ! x right of midpoint
         w(i0:n,1) = x(:,1)
      else
         pd = pn(:,1)
         q = pn(:,0)/(pd*(1.0_float-x(:,0)**2))   ! P/((1-x^2)P')
         xx = x(:,0)
         xx(:nneg) = -xx(:nneg)
         rx(:,0) = jac_wts_rx(ex(:,0),xx)
         rx(:,1) = jac_wts_rx(ex(:,1),xx)
         w(i0:n,0) = x(:,1)      ! x right of midpoint
         w(i0:n,1) = acos(x(:,1))
      endif
      c = Jacobi_Deriv2(ns,a,b)   ! coefficient in derivative relationship
      dw(:,0) = jac_wt_ord(ex(:,0),c,xx,q,ord)
      w(i0:n,2) = cw(0,0)/(rx(:,0)*pd*(dw(:,0)+1.0_float))        ! Wb
      dw(:,1) = jac_wt_ord(ex(:,1),c,xx,q,ord)
      w(i0:n,3) = cw(0,1)/(rx(:,1)*pd*(dw(:,1)+1.0_float))**2     ! W
      if(debug > 2)then
         j = 1 ! max(1,nx-10)
         call arrayprint("Jacobi_Weights: x,p,p',rx,dw:",x(j:,:),pn(j:,:),rx(j:,:),dw(j:,:),fmtx=fmt22)
         j = sym ! max(1,nx-10)
         call arrayprint("Jacobi_Weights: w:",w(j:,:),fmtx=fmt22)
      endif
      call jac_wts_store(n,nx,nneg,g,theta,w)
   End Function Jacobi_Weights
   !----------------------------------------------------------------------------
   Subroutine jac_weight_values(nn,sym,a,b,ex,cw,wbfac)
      integer, intent(in) :: nn,sym
      real(float), intent(in) :: a,b
      integer, intent(out) :: ex(0:,0:)
      real(float), intent(out) :: cw(0:,0:),wbfac
      real(float) :: s0(0:1),p0(0:1),bs,ba,an,zn1,zn8,fshft
      integer :: icase,ia,ib,n
      logical :: Radau
      real(float), parameter :: sq2 = sqrt(2.0_float)
      ! -----------------------------------------------------
      ! returns exponents, constants, end & center weights and
      ! scaling factor for barycentric Weights
      ! parameters index is 0,1 for barycentric, quadrature weights
      ! nn - polynomial degree (must be 2n+ for full shortcut calcs)
      ! a - alpha
      ! b - beta (+/-0.5 for shortcut)
      ! ex - exponents in r(x) function
      ! cw[0,:] - constant in weight
      ! cw[1-3,:] = weight for left end, center & right end
      !             some of these are not used & overwritten
      ! wbfac - scaling factor in barycentric weights
      ! Note: the weights are:
      !  w = cw[0]*[r(x)*p']^(-1 or -2)
      !  r(x) = (1+x)^u*(1-x)^v, where u = ½ex(0) & v = ½ex(1) (see jac_wts_rx)
      ! -----------------------------------------------------
      cw = 0.0_float
      ib = nint(b*2.0_float)+2   ! 1,2,3,4 for b = -0.5,0,+0.5,1
      ia = nint(a)               ! 0,1
      icase = ia*10 + ib         ! 1 to 14
      if(sym > 0)icase = -icase  ! -14 to -1 symmetric
      Radau = icase==04 .or. icase==12
      n = nn;  bs = b;  ba = b
      if(mod(ib,2)==1)then       ! b = +/-0.5, shortcut,planar or spherical
         n = n/2
         if(sym==0)ba = a
      elseif(2*ia == ib-2)then   ! b == a, full ultraspherical
         n = n/2  ! corresponding shortcut
         bs = (2*mod(nn,2)-1)*0.5_float
      endif
      fshft = real(2**nn,float)
      wbfac = Jacobi_Lead(nn,a,ba)  ! wb normalizing factor
      p0 = Jacobi_Ends(nn,a,b)
      if(Radau)then
         s0 = 1.0_float; an = 1.0_float
      else
         s0 = Jacobi_Ends(n,a,bs)
         an = p0(1)/s0(1)
         cw(2,0) = -s0(1)/(s0(0)*p0(1))   ! Wb center
      endif
      if(Debug > 0)then ! -1)then
         write(iout,'(a,2(a1,i5),4(a1,f5.1))')'_weight_value: nn,n,a,b,bs,ba =', &
            tb,nn,tb,n,tb,a,tb,b,tb,bs,tb,ba
         write(iout,'(a,3i3)')'_weight_value: ia,ib,icase =',ia,ib,icase
         write(iout,'(a,8g14.7)')'_weight_value: wbfac,an,p0,s0 =',wbfac,an,p0,s0
      endif
      cw(1,0) = -0.5_float/p0(0)          ! Wb ends
      cw(3,0) = +0.5_float/p0(1)          ! Wb ends
      zn1 = real(nn+1,float)
      zn8 = (8.0_float)*zn1/(zn1 + 1.0_float)
      select case(icase)
      case(02)       ! Gauss full
         ex(:,0) = 2
         ex(:,1) = 1
         cw(0,0) = -1.0_float
         cw(0,1) = 2.0_float
         cw(2,1) = 2.0_float*cw(2,0)**2
      case(01)       ! Gauss shortcut even
         ex(0,0) = 1;  ex(1,0) = 2
         ex(:,1) = 1
         cw(0,0) = -1.0_float/(an*sq2)
         cw(0,1) = 0.5_float
      case(03)       ! Gauss shortcut odd
         ex(:,0) = 2
         ex(0,1) = 2;  ex(1,1) = 1
         cw(0,0) = -1.0_float/an
         cw(0,1) = 1.0_float
         cw(2,1) = 2.0_float*cw(2,0)**2
      case(14)       ! Lobatto full
         ex(:,:) = 2
         cw(0,0) = -1.0_float
         !cw(0,1) = zn8
         cw(:,1) = zn8*cw(:,0)**2
      case(11)       ! Lobatto shortcut even
         ex(0,:) = 1;  ex(1,:) = 2
         cw(0,0) = -1.0_float/(an*sq2)
         !cw(0,1) = zn8*cw(0,0)**2
         !cw(1:3,1) = zn8*cw(1:3,0)**2
         cw(:,1) = zn8*cw(:,0)**2
      case(13)       ! Lobatto shortcut odd
         ex(:,:) = 2
         cw(0,0) = -1.0_float/an
         !cw(0,1) = zn8*cw(0,0)**2
         !cw(1:3,1) = zn8*cw(1:3,0)**2
         cw(:,1) = zn8*cw(:,0)**2
      case(12)       ! Radau right
         ex(:,0) = 2
         ex(0,1) = 1;  ex(1,1) = 2
         cw(0,0) = -1.0_float
         cw(0,1) = 4.0_float
         cw(3,1) = 2.0_float/(p0(1)*p0(1))
      case(04)       ! Radau left
         ex(:,0) = 2
         ex(0,1) = 2;  ex(1,1) = 1
         cw(0,0) = -1.0_float
         cw(0,1) = 4.0_float
         cw(1,1) = 2.0_float/(p0(0)*p0(0))
      case(-10:-1)   ! Gauss symmetric
         ex(0,0) = 0;  ex(1,0) = 2
         ex(:,1) = 1
         cw(0,0) = -1.0_float
         cw(3,0) = (1.0_float)/p0(1)
         cw(0,1) = 0.5_float
         cw(2,1) = (0.5_float)*cw(2,0)**2
         wbfac = wbfac*fshft
      case(-20:-11)  ! Lobatto symmetric
         ex(0,0) = 0;  ex(1,0) = 2
         ex(0,1) = 1;  ex(1,1) = 2
         cw(0,0) = -1.0_float
         cw(3,0) = (1.0_float)/p0(1)
         cw(0,1) = zn1/(zn1 + b)
         cw(3,1) = cw(0,1)/(p0(1)*p0(1)*2.0_float)
         wbfac = wbfac*fshft
      end select
      if(debug > 0)then ! -1)then
         write(iout,'(a,4i3,a,2g14.7)')'_weight_value: exponents =',ex,', wbfac =',wbfac
         call arrayprint('_weight_value:',transpose(cw))
      endif
   End Subroutine jac_weight_values
   !----------------------------------------------------------------------------
   Subroutine jac_wts_store(n,nx,nneg,g,theta,w)
      integer, intent(in) :: n,nx,g
      integer, intent(inout) :: nneg
      logical, intent(in)  :: theta
      real(float), intent(inout) :: w(max(0,min(1,g)):n+1,0:3)
      integer :: i0,i1,nodd,k,sym
      ! ------------------------------------------------------------
      i0 = n+1-nx             ! 1st unique: 1 or 1st to right of midpoint
      nodd = mod(n,2)
      sym = max(0,min(1,g))
      if(Debug > 0)write(iout,'(a,10i5)')'jac_wts_store n,nx,nneg,i0,g =',n,nx,nneg,i0,g
      if(g == Cylindrical .and. n > nx)then  ! Gauss-cylindrical, special processing
         if(Debug > 0)call arrayprint('Cyl G: top',w(i0:,:),fmtx=fmt17)
         w(i0-1,1) = pi2
         w(1:nx,3) = w(n:i0:-1,3)
         if(theta)then
            w(i0-1,0) = sqrt(0.5_float)
            w(1:nx,1) = w(n:i0:-1,1)
            w(1:nx,0) = sin(w(1:nx,1)*0.5_float)
            w(i0:n,0) = cos(w(i0:n,1)*0.5_float)
            w(1:nx,2) = w(n:i0:-1,2)*tan(w(1:nx,1))**2   ! calculate lower from upper half
         else
            w(i0-1,0) = 0.5_float
            w(1:nx,1) = w(n:i0:-1,1)
            w(1:nx,0) = (1.0_float-w(n:i0:-1,0))*0.5_float
            w(i0:n,0) = (1.0_float+w(i0:n,0))*0.5_float
            w(1:nx,2) = w(n:i0:-1,2)*w(1:nx,0)/w(n:i0:-1,0)
            w(:,0) = sqrt(w(:,0))
         endif
         nneg = nx
         if(nodd == 0)w(1:nx,2) = -w(1:nx,2)  ! fix sign of barycentric weights
         if(debug > 0)call arrayprint('Jacobi_Weights: Gcyl w',w,fmtx=fmt17)
      elseif(g .ne. 0)then    ! other symmetric and shortcut
         i0 = n+1-nx;  i1 = i0+nneg-1
         if(theta)then                 ! will need some work to use
            w(i0:i1,0) = sin(w(i0:i1,1)*0.5_float)
            w(i1+1:,0) = cos(w(i1+1:,1)*0.5_float)
         else
            w(i0:i1,0) = sqrt((1.0_float-w(i0:i1,0))*0.5_float)
            w(i1+1:,0) = sqrt((1.0_float+w(i1+1:,0))*0.5_float)
         endif
         if(g < 0)then
            w(i0:,1) = w(i0:,1)*0.5_float ! angles still contain negatives
            w(i0:i1,1) = pi2 - w(i0:i1,1) ! shortcut, avoid folding twice
         endif
      endif
      if(n > nx .and. g < 1)then ! reflect ultraspherical
         nneg = nx + 1
         i0 = n+1-nx;  i1 = n+1-sym
         do k = 0,size(w,2)-1
            w(sym:nx,k) = w(i1:i0:-1,k)
         enddo
         if(nodd == 0)w(sym:nx,2) = -w(sym:nx,2)
      elseif(nneg == 0 .and. g < 1)then   ! not tracking negatives
         w(0,0) = -1.0_float
         w(0,1) = pi
      elseif(g <= 0)then   ! +1 for left end, nonsymmetric
         nneg = nneg + 1
      endif
   end subroutine jac_wts_store
   !----------------------------------------------------------------------------
   function jac_wts_rx(ex,x) result(rx)
      real(float) :: x(:),rx(size(x))
      integer :: ex(:)
      ! calculate rx = [(1+x)**ex1/2]*[(1-x)**ex2/2]
      Select case(ex(1)*10+ex(2))
      case(02)
         rx = 1.0_float - x
      case(11)
         rx = sqrt(1.0_float - x*x)
      case(12)
         rx = sqrt((1.0_float - x)*(1.0_float - x*x))
      case(21)
         rx = sqrt((1.0_float + x)*(1.0_float - x*x))
      case(22)
         rx = 1.0_float - x*x
      End Select
   end function jac_wts_rx
   !----------------------------------------------------------------------------
   Function jac_wts_rxt(ex,xt,kneg)  result(rx)
      integer, intent(in) :: ex(:),kneg
      real(float), intent(in) :: xt(:)
      real(float) :: rx(size(xt))
      real(float), parameter :: sq2 = sqrt(2.0_float)
      ! return r(x) parameters for Weights
      ! calculate rx = [(1+x)**ex1/2]*[(1-x)**ex2/2]
      ! xt is folded theta angle
      Select case(ex(1)*10+ex(2))
      case(02)
         rx(1:kneg) = cos(xt(:kneg)*0.5_float)    ! sqrt(1 - x) (folded)
         rx(kneg+1:) = sin(xt(kneg+1:)*0.5_float) ! sqrt(1 - x)
         rx = (2.0_float)*rx**2                   ! 1 - x
      case(11)
         rx = sin(xt)                             ! sqrt(1-x^2)
      case(12)
         rx(1:kneg) = cos(xt(:kneg)*0.5_float)    ! sqrt(1 - x) (folded)
         rx(kneg+1:) = sin(xt(kneg+1:)*0.5_float) ! sqrt(1 - x)
         rx = sq2*rx*sin(xt)                 ! sqrt[(1-x^2)(1-x)]
      case(21)
         rx(1:kneg) = sin(xt(:kneg)*0.5_float)    ! sqrt(1 + x) (folded)
         rx(kneg+1:) = cos(xt(kneg+1:)*0.5_float) ! sqrt(1 + x)
         rx = sq2*rx*sin(xt)                      ! sqrt[(1-x^2)(1+x)]
      case(22)
         rx = sin(xt)**2                          ! 1 - x^2
      End Select
   End Function jac_wts_rxt
      ! ------------------------------------------------
   Function jac_wt_ord(ex,c,x,q,ord)   result(w)
      integer, intent(in) :: ex(:),ord
      real(float), intent(in) :: c(3),x(:),q(:)
      real(float) :: w(size(x))
      real(float) :: u,v,d(3,0:3),dx(size(x),3)
      ! Expansion of Weight function for [(1+x)^u*(1-x)^v]*P'
      d = 0.0_float
      u = real(ex(1),float)*0.5_float
      v = real(ex(2),float)*0.5_float
      if(ord > 1)then
         d(1,0) = v - u - c(2) ! v - u - (alpha-beta)
         d(1,1) = v + u - c(1) ! v + u - (alpha+beta+2)
      endif
      if(ord > 2)then
         d(2,0) = (v - u)*d(1,0) - d(1,1) - c(3)
         d(2,1) = (v - u)*d(1,1) + (u + v - 2.0_float)*d(1,0)
         d(2,2) = (u + v - 1.0_float)*d(1,1) + c(3)
      endif
      if(ord > 3)then
         d(3,0) =(c(2)*d(1,0)+d(2,0))*(v-u) - (c(3)+d(1,1))*c(2)-d(2,1)
         d(3,1) =(c(1)*d(1,0)+c(2)*d(1,1)+d(2,1))*(v-u)+(c(2)*d(1,0)+d(2,0))*(v+u)  &
                -(c(3)+d(1,1))*c(1)-2.0_float*(c(2)*d(1,0)+d(2,2))-4.0_float*d(2,0)
         d(3,2) =(c(1)*d(1,1)+d(2,2))*(v-u)+(c(1)*d(1,0)+c(2)*d(1,1)+d(2,1))*(u+v)  &
               - 2.0_float*c(1)*d(1,0)+(c(3)-d(1,1))*c(2)-3.0_float*d(2,1)
         d(3,3) = (c(1)*d(1,1)+d(2,2))*(u+v) + (c(3)-d(1,1))*c(1) - 2.0_float*d(2,2)
      endif
      d(2,:) = d(2,:)*0.5_float
      d(3,:) = d(3,:)/6.0_float
      select case(ord)
      case(1)
         w = 0.0_float
      case(2)
         w = q*(d(1,0) + d(1,1)*x)
      case(3)
         w = q*(d(1,0)+d(1,1)*x + q*(d(2,0)+x*(d(2,1)+x*d(2,2))))
      case default
         w = q*(d(1,0)+d(1,1)*x + q*(d(2,0)+x*(d(2,1)+x*d(2,2))  &
            + q*(d(3,0)+x*(d(3,1)+x*(d(3,2)+x*d(3,3))))))
      end select
      if(Debug > +2)then
         dx(:,1) = d(1,0) + d(1,1)*x
         dx(:,2) = d(2,0)+x*(d(2,1)+x*d(2,2))
         dx(:,3) = d(3,0)+x*(d(3,1)+x*(d(3,2)+x*d(3,3)))
         write(iout,'(a,i2,a,2i2,a,3g12.5)')'jac_wt_ord: ord =',ord,' ex =',ex,' Deriv2 =',c
         call ArrayPrint('jac_wt_ord: d',d(1:ord-1,0:ord-1))
         Call ArrayPrint('jac_wt_ord: x,v,q,dx',x,w,q,dx(:,1:ord-1))
      endif
      !w = w + 1.0_float
   end function jac_wt_ord
   !----------------------------------------------------------------------------
   Subroutine Jacobi_Pn_update(n,a,b,xt,dxt,pn,nneg,theta,OrdP)
      integer, intent(in) :: n,nneg
      Real(float), intent(in) :: a,b,xt(:),dxt(:)
      Real(float), intent(inout) :: pn(:,0:)
      logical, intent(in) :: theta
      integer, optional :: OrdP
      integer :: ord,ip,ib,nn
      real(float) :: x(size(xt)),dx(size(xt)),r2(size(xt))
      integer, parameter :: iord(3) = (/35,9,0/)
      ! --------------------------------------------
      ! This function corrects approximate values of Pn & Pn' using Taylor series.
      !  n,a,b - polynomial pn^(a,b)
      !  xt - which were used to calculate pn, i.e. pn(xt) are supplied
      !  dxt - pn values are calculated at xt + dxt, i.e. pn(xt+dxt) are returned
      !  Ordp - 1,2 or 3 terms used
      !  p = pn + pn'*dx + pn''*dx^2/2 ...
      !  p'= pn' + pn''*dx + pn'''*dx^2/2 ...
      ! --------------------------------------------
      ib = nint(b*2.0_float)  ! -1 & 1 or shortcut, 0 & 2 are normal cases
      nn = n + mod(ib+2,2)*(n + (ib+1)/2) ! n, 2n, 2n+1 for b = 0,1,-/+0.5
      ord = IFirstLT(iord,nn)
      if(present(OrdP))ord = min(max(1,OrdP),3)
      if(Debug > 0)write(iout,'(a,3i5,i2,2f7.1,L2)') &
         '_Pn_update: n,nn,nneg,ord,a,b,theta =',n,nn,nneg,ord,a,b,theta
      if(theta)then
         dx = -2.0_float*sin(xt+dxt*0.5_float)*sin(dxt*0.5_float)
         r2 = sin(xt);  x = cos(xt)
         pn(:,1) = -pn(:,1)/r2
         r2 = r2*r2
         !call arrayprint('_update',x,dxt,dx,r2,fmtx=fmt22)
      else
         r2 = 1.0_float - xt*xt
         x = xt;  dx = dxt
      endif
      if(nneg > 0)then
         ip = mod(n+1,2)   ! point to the odd one p or p'
         pn(:nneg,ip) = -pn(:nneg,ip)  ! flip it
         call Jacobi_Pn_update_X(n,b,a,x(1:nneg),dx(1:nneg),r2(1:nneg),pn(1:nneg,:),ord)
         pn(:nneg,ip) = -pn(:nneg,ip)  ! flip it back
         call Jacobi_Pn_update_X(n,a,b,x(nneg+1:),dx(nneg+1:),r2(nneg+1:),pn(nneg+1:,:),ord)
      else
         call Jacobi_Pn_update_X(n,a,b,x,dx,r2,pn,ord)
      endif
      if(theta)then
         pn(:,1) = -pn(:,1)*sin(xt+dxt)
      endif
      End Subroutine Jacobi_Pn_update
   !----------------------------------------------------------------------------
   Subroutine Jacobi_Pn_update_X(n,a,b,x,dx,dr2,pn,ord)
      integer, intent(in) :: n,ord
      Real(float), intent(in) :: a,b,x(:),dx(:)
      Real(float), intent(inout) :: dr2(:),pn(:,0:)
      real(float) :: c2(0:2),dpn(size(x),0:2)
      real(float), parameter :: div3 = 1.0_float/3.0_float
      ! --------------------------------------------
      ! This function corrects approximate values of Pn & Pn' using Taylor
      ! series.
      !  n,a,b - polynomial pn^(a,b)
      !  x - values used to calculate pn
      !  dx - change in x
      !  dr2 = 1 - x^2, then dx/(1-x^2)
      !  Ordp - 1,2 or 3 terms used
      !  p and p' evaluated at x - dx
      !  p = pn - pn'*dx + pn''*dx^2/2 ...
      !  p'= pn' - pn''*dx + pn'''*dx^2/2 ...
      ! --------------------------------------------
      c2 = Jacobi_Deriv2(n,a,b)
      if(Debug > 1)write(iout,'(a,i5,i2,2f5.1,2f7.2,g14.6)') &
         '#_Pn_update_X: n,ord,a,b,deriv2',n,ord,a,b,c2
      !call arrayprint('_Pn_update_X x,dx,dpr,pn',x,dx,dpr,pn,fmtx=fmt22)
      dpn(:,0) = ((c2(0)*x + c2(1))*pn(:,1) + c2(2)*pn(:,0))
      if(ord > 1)then
         c2(2) = c2(2) + c2(0);   c2(0) = c2(0) + 2.0_float
         dpn(:,1) = ((c2(0)*x + c2(1))*dpn(:,0) + c2(2)*dr2*pn(:,1))
      endif
      if(Ord > 2)then
         c2(2) = c2(2) + c2(0);   c2(0) = c2(0) + 2.0_float
         dpn(:,2) = ((c2(0)*x + c2(1))*dpn(:,1) + c2(2)*dr2*dpn(:,0))
      endif
      if(Debug > +3)then
         call arrayprint("roots_poly: x,dpr,p,p',p'',p'''",x,dr2,pn,dpn(:,:ord-1),fmtx=fmt22)
      endif
      dr2 = dx/dr2
      Select Case(ord)
      case(1)
         pn(:,0) = pn(:,0) + dx*(pn(:,1) + 0.5_float*dr2*dpn(:,0))
         pn(:,1) = pn(:,1) + dr2*dpn(:,0)
      case(2)
         pn(:,0) = pn(:,0) + dx*(pn(:,1) + 0.5_float*dr2*(dpn(:,0) + div3*dr2*dpn(:,1)))
         pn(:,1) = pn(:,1) + dr2*(dpn(:,0) + 0.5_float*dr2*dpn(:,1))
      case(3)
         pn(:,0) = pn(:,0) + dx*(pn(:,1) + 0.5_float*dr2*(dpn(:,0) + div3*dr2*dpn(:,1)))
         pn(:,1) = pn(:,1) + dr2*(dpn(:,0) + 0.5_float*dr2*(dpn(:,1) + div3*dr2*dpn(:,2)))
      end select
   End Subroutine Jacobi_Pn_update_X
   !----------------------------------------------------------------------------
   Function RootEstimate0(n,nx,a,b,Shift)   result(x)
      Integer, intent(in) :: n,nx
      Real(float), intent(in) :: a,b
      logical, optional :: Shift
      Real(float) :: x(nx),zn(nx),abn
      logical :: shft
      Integer :: k
      ! -------------------------------
      ! Simple root estimates using generalized Chebyshev formula
      ! -------------------------------
      shft = Default_shift;  if(Present(Shift))Shft = Shift
      zn = (/(Real(2*k,float),k=n-nx+1,n)/)
      abn = Real(2*n+1,float) + a + b
      zn = pi*(zn + b - 0.5_float)/abn
      x = -cos(zn)
      if(shft)x = (x + 1.0_float)*0.5_float
   End Function RootEstimate0
   !----------------------------------------------------------------------------
   Function RootEstimate(n,nx,a,b,nneg)  result(x)
      Integer, intent(in) :: n,nx
      Real(float), intent(in) :: a,b
      integer, optional :: nneg
      Real(float) :: x(nx),an(n),ae,xe,zn
      integer :: nn,ib,bmeth,i1,i2,ix,k
      Logical :: Legendre,Lobatto,Ultraspherical,Shortcut
      ! roots (theta) for n = 2 & 3:
      !  ------------------------------------------------
      !  Estimate roots of Jacobi polynomials using up to 4 different
      !  approximations (no more than 2 for a given set of roots).
      !  It is optimized for Gauss, Lobatto & Radau base points.
      !  Cutoff values are used to chose between boundary & interior methods.
      !  Boundary correlation used on both ends if not Ultraspherical
      !  for any given set of points.
      !   n - is polynomial degree
      !   nn = n for full, 2n+ for shortcut, i.e. corresponding full
      !   nx - roots determined, either n or n/2 (if ultraspherical & not shortcut))
      !   i1:i2 - range of points using interior correlation
      !   xe - a cutoff value between boundary and interior correlation
      !   ae - cutoff angle arccos(xe)
      !   bmeth - boundary method to use
      !  Some logicals used are:
      !   Shortcut = T when b = +/-0.5
      !   Legendre = T if Legendre (including Shortcut for Legendre)
      !   Ultraspherical - T if Ultraspherical (including Shortcut)
      !  ------------------------------------------------
      i1 = 1;  bmeth = 2;
      ib = nint(b*2.0)  ! -1,0,1,2 ?
      Shortcut = abs(ib) == 1 .and. abs(2.0*b-Real(ib)) < 1.0e-6
      Ultraspherical = (a == b .or. Shortcut)
      Legendre = (a == 0.0_float .and. Ultraspherical);
      Lobatto  = (a == 1.0_float .and. Ultraspherical);
      nn = n
      if(Shortcut)then
         i2 = n
         nn = 2*n + (ib+1)/2     ! nn = n of corresponding full polynomial
      elseif(Ultraspherical)then
         i2 = int(n/2)
      endif
      zn = real(nn,float)
      !write(iout,'(a,3i5,2f5.1,4L2)')'# RootEstimate:',n,nx,nn,a,b,Ultraspherical,Legendre,Lobatto,Shortcut
      x = 0.0_float;
      an = (/(Real(2*k,float),k=n-1,0,-1)/)
      an = pi*(an+a+1.5_float)/(Real(2*n+1,float)+a+b) ! rough estimate for cutoffs
      if(Legendre)then     ! Legendre a=b=0
         if(nn < 19)then   ! linear
            xe = 0.50_float + zn*(0.019_float)
            bmeth = 1
         else     ! large n, bmeth = 2
            xe = 0.57_float - zn*(0.001_float)
            xe = max(xe,0.480_float)
         endif
         !write(iout,'(a,g12.4,i3)')'# Estimates Right end :',xe,bmeth
      elseif(Lobatto)then ! Ultraspherical - Lobatto
         xe = 0.813_float - zn*(0.00075_float)
         xe = max(xe,0.780_float)
         if(nn < 5)xe = 0.60_float  ! small nn
         if(nn < 20)bmeth = 1       ! change boundary methods
      elseif(UltraSpherical)then ! Ultraspherical - not Legendre or Lobatto
         x = an(n-nx+1:)         ! use crude estimate
         return
      else  ! Radau - boundary method througout (easily modified)
      !  xe = -0.5;   ae = acos(xe);  ix = FirstLT(a0,ae)-1;  !uncomment to use
         xe = (0.10_float)*(b-a);   ae = acos(xe);  ix = FirstLT(an,ae)-1
         x(n:n-ix+1:-1) = pi - RootEstBoundary(n,ix,b,a,bmeth);
         !write(iout,'(a,2i6,2g12.4,2i5)')'# Estimates Left end :',n-ix+1,n,xe,ae,ix
         xe = -0.8; ! positive value for interior method in center
         i2 = n-ix;
      endif
      if(Shortcut)xe = 2.0*xe*xe-1.0   ! convert to shortcut values
      ae = acos(xe); i1 = min(i2+1,n-FirstLT(an,ae)+2);  ! last boundary point
      if(nn > 3)then
         x(i1:i2)  = RootEstInterior(n,i1,i2,a,b);
         x(1:i1-1) = RootEstBoundary(n,i1-1,a,b,bmeth)
         if(Debug > 1)then
            write(iout,'(a,2i5,2g12.4,i3)')'# Estimates Right end :',i1,i2,xe,ae,bmeth
            call vectorprint('Interior:',x(i1:i2))
            call vectorprint('Boundary:',x(1:i1-1))
         endif
      else
         x = RootStored(nn,nx)   ! exact roots n < 4
      endif
      if(Ultraspherical .and. nx == n .and. .not.Shortcut)then
         ! copy right half to left half
         x(i2+1) = pi*0.5_float
         x(n:n-i2+1:-1) = x(1:i2)
         x(1:i2) = pi - x(1:i2) ! for theta
      else    ! reverse values, so theta goes large to small (x small to large)
         an(1:nx) = x;    x = an(nx:1:-1)
      endif
      !if(present(nneg) .and. .not.Ultraspherical)then
      if(present(nneg))then
         nneg = FirstLT(x,pi2) - 1
         x(:nneg) = pi - x(:nneg)
      endif
      contains
      ! -------------------------------------
      Function RootStored(nn,nx)  result(x)
         integer :: nn,nx,i1
         real(float) :: x(nx)
         ! The first 3 roots (theta) for Gauss, Lobatto and Radau
         real(float), parameter :: Groots(2) = (/0.955316618124509278D0,0.684719203002282914D0/)
         real(float), parameter :: Lroots(2) = (/1.107148717794090503D0,0.857071947850130988D0/)
         real(float), parameter :: Rroots(6) = (/  &
                          1.910633236249018556D0,1.276676121181329395D0,2.332144396859698049D0, &
                          0.957802316375337495D0,1.752866861904567592D0,2.537158998077295037D0/)
         if(Legendre)then
            if(nn > 1)x(1) = Groots(nn-1)
         elseif(Lobatto)then
            if(nn > 1)x(1) = Lroots(nn-1)
         else
            i1 = max(1,2 + (nn-2)*2)
            if(a > b)then
               x = Rroots(i1:i1+nx-1)
            else
               x = pi - Rroots(i1+nx-1:i1:-1)
            endif
         endif
         if(Shortcut)x = (2.0_float)*x  ! theta for shortcut
      End Function RootStored
   End Function RootEstimate
   !----------------------------------------------------------------------------
   Function RootEstInterior(n,i1,i2,a,b,method)  result(x)
      Integer, intent(in) :: n,i1,i2
      Real(float), intent(in) :: a,b
      Integer, optional :: method
      Integer :: meth,ib,k,nn
      Real(float) :: x(i1:i2),r,a0(i1:i2),d0(i1:i2)
      logical :: Legendre,Shortcut,Ultraspherical
      ! ----------------------------------------------------
      ! Asymptotic correlation for roots accurate away from boundaries
      !  1. Gatteschi and Pittaluga [see Gautschi & Giordano, Hale & Townsend]
      !  2. Tricomi valid only for Legendre case.
      ! Only unique values are calculated for a=b (Ultraspherical case)
      ! Roots are ordered large to small
      !  n - polynomial
      !  nn - corresponding full polynomial, i.e. n if no shortcut, or 2*n+
      ! ----------------------------------------------------
      x = 0.0_float
      if(i1 > i2)return
      meth = 2;  if(present(method))meth = method
      ib = nint(b*2.0_float)  ! -1,0,1,2 ?
      Shortcut = abs(ib) == 1
      Ultraspherical = (a == b .or. Shortcut)
      Legendre = (a == 0.0_float .and. Ultraspherical);
      nn = n
      if(Shortcut)nn = 2*n + (ib+1)/2   ! nn = n of full polynomial
      if(.not.Legendre)meth = 1
      !write(iout,'(a,3i4,2f7.2,i3)')'# RootEstInterior: n,i1,i2,a,b,meth >',n,i1,i2,a,b,meth
      !write(iout,'(a,3L2)')'# RootEstInterior: Ultra, Legendre, Shortcut >',Ultraspherical,Legendre,Shortcut
      Select Case(meth)
      Case(2)
         Call RootEstLeg(nn)
         if(Shortcut)x = (2.0_float)*x  ! theta for shortcut
      Case Default
         r = (1.0_float)/(Real(2*n,float)+a+b+1.0_float);
         x = (/(Real(2*k,float),k=i1,i2)/)
         a0 = (x + a - 0.5_float)*r*pi
         d0 = tan(0.5_float*a0);
         d0 = ((0.25_float-a*a)/d0 - (0.25_float-b*b)*d0)*r*r;
         x = a0 + d0
      End Select
      ! -----------------------------------
      Contains
      ! -----------------------------------
      Subroutine RootEstLeg(n)
         integer, intent(in) :: n
         real(float) :: xn,xn1,xn4,z1
         real(float), parameter :: x39 = (39.0_float)/(384.0_float)
         real(float), parameter :: x28 = (28.0_float)/(384.0_float)
         ! Triconi interior correlation specific for Legendre roots
         x = (/(Real(k,float),k=i1,i2)/)
         xn = Real(n,float);  xn1 = xn - 1.0_float;  xn4 = (1.0_float)/xn**4
         z1 = 1.0_float - ((0.125_float)*xn*xn1+x39)*xn4;
         a0 = pi*(x - 0.25_float)/(xn + 0.5_float);
         d0 = (x28*xn4)/(sin(a0)**2);
         x = (z1 + d0)*cos(a0);
         x = acos(x)
      End Subroutine RootEstLeg
   End Function RootEstInterior
   !----------------------------------------------------------------------------
   Function RootEstBoundary(n,nc,a,b,method)  result(x)
      Integer, intent(in) :: n,nc
      Real(float), intent(in) :: a,b
      Integer, optional :: method
      Real(float) :: x(nc)
      Integer :: meth
      Real(float) :: r,a2,b2,v,xtan(nc),xcot(nc)
      Real(float), parameter :: c45 = 1.0_float/45.0_float, c33 = 1.0_float/3.0_float
      ! ------------------------------------------------
      ! Boundary Estimation methods. Using Bessel zeroes
      !  1. Gatteschi [Gautschi and Giordano (2008)]
      !  2. Olver, more uniform and better for x < 0.80 (appx.)
      !  n - no. roots, nc - no. to estimate (near boundary)
      ! ------------------------------------------------
      if(nc <= 0)return
      meth = 2
      if(present(method))meth = method
      r = Real(2*n,float)+a+b+1.0_float;
      a2 = a*a;  b2 = b*b;
      x = RootEstBessel(nc,a);
      Select Case(meth)
      Case(1)
         v = (1.0_float)/sqrt(r*r + c33 - c33*a2 - b2);
         r = c45*(4.0_float-a2-(15.0_float)*b2)*v**4;
         xcot = (2.0_float)*x*v
         x = (2.0_float)*x*v*(1.0_float - r*((0.5_float)*x**2+a2-1.0_float));
         !call arrayprint('Boundary 1',xcot,x,fmtx=fmt17)
      Case Default
         r = (1.0_float)/r;  b2 = a2-b2
         a2 = (2.0_float)*a2-0.5_float
         x = x*r*2.0_float;
         xtan = tan(x*0.5_float);
         xcot = (1.0_float/x - (0.5_float)*(1.0_float/xtan-xtan))
         xcot = (a2*xcot - b2*xtan)*r*r;
!         call arrayprint('Boundary 2',x,x+xcot,fmtx=fmt17)
         x = x + xcot
      End Select
   End Function RootEstBoundary
   !----------------------------------------------------------------------------
   Function RootEstBessel(n,a)  result(xj)
      Integer, intent(in) :: n
      Real(float), intent(in) :: a
      Real(float) :: xj(n)
      Integer :: j,k,ia
      Logical :: Intorder
      Real(float) :: xmu,c(4),r(n)
      Real(float), parameter :: x75 = (1.0_float)/(0.75_float)
   !  Estimation roots of first n Bessel function of kind "a". 16 are stored for order 0 & 1
   !  Others calculated using first 5 terms of McMahon's approximation, ~2e-16 is the maximum error
      Real(float), Parameter :: cf(2:4) = (/1.0_float/6.0_float,1.0_float/30.0_float,0.125_float/105.0_float/)
      real(float), parameter :: c2(2) = (/-31.0_float,7.0_float/)
      real(float), parameter :: c3(3) = (/3779.0_float,-982.0_float,83.0_float/)
      real(float), parameter :: c4(4) = (/-6277237.0_float,1585743.0_float,-153855.0_float,6949.0_float/)
      Integer, Parameter :: Mxj=16
      Real(float),Parameter :: zj0(Mxj) = &
       (/ 2.40482555769577276862163188_float, 5.52007811028631064959660411_float, 8.65372791291101221695419871_float, &
         11.79153443901428161374304491_float,14.93091770848778594776259400_float,18.07106396791092254314788298_float, &
         21.21163662987925895907839335_float,24.35247153074930273705794476_float,27.49347913204025479587728823_float, &
         30.63460646843197511754957893_float,33.77582021357356868423854635_float,36.91709835366404397976949306_float, &
         40.05842576462823929479930737_float,43.19979171317673035752407273_float,46.34118837166181401868578888_float, &
         49.48260989739781717360276153_float /)
      Real(float),Parameter :: zj1(Mxj) = &
       (/ 3.83170597020751231561443589_float, 7.01558666981561875353704998_float,10.17346813506272207718571178_float, &
         13.32369193631422303239368413_float,16.47063005087763281255246047_float,19.61585851046824202112506588_float, &
         22.76008438059277189805300515_float,25.90367208761838262549585545_float,29.04682853491685506664781988_float, &
         32.18967991097440362662298410_float,35.33230755008386510263447902_float,38.47476623477161511205219756_float, &
         41.61709421281445088586351681_float,44.75931899765282173277935271_float,47.90146088718544712127400872_float, &
         51.04353518357150946873303463_float /)
      Real(float), parameter :: xjz(Mxj,0:1) = Reshape((/zj0,zj1/),(/MxJ,2/))
      !----------------------------------------------------------
      j=1;  ia = max(0,min(1,NINT(a)));  Intorder = abs(Real(ia,float)-a) < 1.0e-6
      !Intorder = .false.  ! for testing McMahon expansion
      if(Intorder)then
         j = min(Mxj,n);  xj(:j) = xjz(:j,ia)
         if(j == n)return
         j = j+1
      endif
      xmu = 4.0*a*a
      c(1) = (xmu - 1.0_float)
      c(2) = cf(2)*c(1)*(c2(1) + xmu*c2(2))
      c(3) = cf(3)*c(1)*(c3(1) + xmu*(c3(2) + xmu*c3(3)))
      c(4) = cf(4)*c(1)*(c4(1) + xmu*(c4(2) + xmu*(c4(3) + xmu*c4(4))))
      xj(j:) = (0.5_float*a - 0.25_float + Real((/(k,k=j,n)/),float))*pi
      r(j:) = 0.125_float/(xj(j:)**2)
      xj(j:) = xj(j:)*(1.0_float - r(j:)*(c(1) + r(j:)*(c(2) + r(j:)*(c(3) + r(j:)*c(4)))))
   End Function RootEstBessel
   !----------------------------------------------------------------------------
   Function Jacobi_Pm(tx,n,nd,alpha,beta,Method,Order,Theta,nneg)  result(p)
      real(float), intent(in) :: tx(:)
      integer, intent(in) :: n,nd
      real(float), optional :: alpha,beta
      integer, optional :: Method,Order,nneg
      logical, optional :: Theta
      real(float) :: p(size(tx),0:nd),x(size(tx)),t1(size(tx)),t2(size(tx)),cs(size(tx)),sn(size(tx))
      integer :: i,meth,nx,kneg,nn,nb,icase
      logical :: thet,Recur
      real(float) :: abx(2),a,b
      ! --------------------------------------
      ! calculates Jacobi Polynomials, Pn  n0 to n together with nd derivatives
      ! results are p(:,n0-n) to p(:,nd), i.e. p(:,0) = Pn & p(:,1) is Pn'
      ! uses either recurrence relations or asymptotic approximation
      !  tx - either x or theta angle, acos(x) when theta = .true.
      !  n  - polynomial number
      !  nd - 0 - no. derivative, 1 1st derivativs calculated
      !  alpha,beta - parameters in Jacobi weight (optional)
      !  Method - method of calculation, code determines if not supplied
      !     =-1 recurrent using shortcut if possible
      !     = 0 recurrent, full calculations
      !     = 1 asymptotic interior method
      !     = 2 asymptotic boundary method
      !     =-2 or not supplied, code decides method(s)
      !  Order - optional order of asymptotic
      !  Theta - .true. for theta supplied (tx) & dPn/dtheta returned
      !  nneg - number negative roots (Radau), when postive are supplied
      !         this option implements folding of negative roots about pi/2
      ! Given tx, theta & nneg, first determine
      !  x - values of x (including any negatives)
      !  tf - folded theta, i.e. 0 to pi/2
      !  t - unfolded theta, i.e. 0 to pi
      ! --------------------------------------
      nx = size(tx)
      abx = abparm(alpha,beta); a = abx(1);  b = abx(2)
      meth = -2;        if(Present(Method))meth = Method
      thet = .false.;   if(Present(Theta))thet = Theta
      kneg = 0;         if(present(nneg))kneg = nneg
      ThetaFac = 1.0_float
      if(Debug > 1)write(iout,'(a,2i6,2i3,L2)') &
         '# Jacobi_Pm : n,nneg,nd,meth,thet =',n,kneg,nd,meth,thet
      p = 0.0_float
      nn = n
      nb = nint(b*2.0_float);   if(abs(nb)==1)nn = 2*n + (nb+1)/2
      Recur = (nn <= Nboundary .and. meth == -2) .or. meth == 0 .or. meth == -1
      if(thet)then
         x = cos(tx)
         t1 = tx
      else
         x = tx
         t1 = acos(tx)
      endif
      t2 = t1
      if(Recur)then  ! only recurrent is used
         x(:kneg) = -x(:kneg)  ! kneg > 0 only if folded x supplied
         p = Jacobi_Recurrent(x,n,nd,alpha,beta)
         if(nd > 0 .and. thet)p(:,1) = -p(:,1)*sqrt(1.0_float-x*x)!  *sin(t)
      else
         ! cases:
         !  0  (0,0) Gauss
         !  11 (1,1) Lobatto
         !  1,10  (1,0) RadauR or (0,1) RadauL
         !  5,6   (0 or 1,+/-0.5) Shortcut Gauss or Lobatto
         ! note also:
         icase = nint(a) + 5*abs(nb)  ! = a + abs(b)*10
         select case(icase)
         case(5,6)  ! b = +/-0.5
            JacobiType = ShortGauss + nint(a)
            ThetaFac = 2.0_float
            t1 = t1*0.5_float ! halve values for calculations
            t2 = t2*0.5_float ! (pi - theta)/2 for negatives
            if(kneg > 0)then
               cs(1:kneg) = sin(t1(1:kneg))  ! cosine for negatives
               cs(kneg+1:)= cos(t1(kneg+1:)) ! cos(0.5*pi-x) = sin(x)
               sn(1:kneg) = cos(t1(1:kneg))  ! sin for negatives
               sn(kneg+1:)= sin(t1(kneg+1:))
               x(:kneg) = -x(:kneg)
               t1(:kneg) = pi2 - t1(:kneg)      ! negative values reversed
               t2(kneg+1:) = pi2 - t2(kneg+1:)  ! for interior asymptotic calculations
            endif
         case(1,10)  ! Radau, must be folded
            JacobiType = RadauR + nint(b)
            if(kneg == 0)then
               kneg = FirstGT(x(:),0.0_float) - 1 ! fold manually
               x(:kneg) = -x(:kneg)
               t1(:kneg) = pi - t1(:kneg)
            endif
            t2 = pi2 - t1  ! for interior asymptotic
            cs = cos(t1)   ! ok, since reversal used for negatives
            sn = sin(t1)
         case default   ! full Gauss and Lobatto must not be folded
            JacobiType = Gauss + nint(a)
            if(kneg > 0)then  ! shouldn't happen for Gauss or Lobatto
               x(:kneg) = -x(:kneg)
               t1(:kneg) = pi - t1(:kneg)
               kneg = 0
            endif
            t2 = pi2 - t1  ! for interior asymptotic
            cs = cos(t1)
            sn = sin(t1)
         end select
         p = Jacobi_Pmx(x,t1,t2,cs,sn,n,nd,a,b,kneg,meth,Order)
         if(nd > 0 .and. .not.thet)p(:,1) = -p(:,1)/sin(t1)   ! convert to x derivative
      endif
      if(Debug > 2)then
         i = max(1,nx-10)
         if(Debug > 3)i = 1
         call ArrayPrint('Jacobi_Pm:',x(i:nx),t1(i:nx),t2(i:nx),p(i:nx,:),fmtx=fmt22)
      endif
   End Function Jacobi_Pm
   !----------------------------------------------------------------------------
   Function Jacobi_Pmx(x,t1,t2,cs,sn,n,nd,a,b,nneg,meth,Order)  result(pn)
      Real(float), intent(inout) :: x(:),t1(:),t2(:),cs(:),sn(:)
      Integer, intent(in) :: n,nd,meth,nneg
      Real(float), intent(in) :: a,b
      Integer, optional :: Order
      Real(float) :: pn(size(t1),0:nd)
      Real(float) :: p(size(t1),0:1),ab(2,2)
      integer :: nt,j,k,i0,i1,icase,iultra(2)
      logical :: rev(2)
      ! ------------------------------------------
      ! Handles reversal of folded negative points - reordered in increasing
      ! absolute values of x.
      !  Pcalc called if all points use same method
      !  Pcomposite called if variety of methods used
      !     x, t - points in x and theta (values are folded)
      !     n - polynomial degree
      !     nd - number of derivatives
      !     a,b - Jacobi alpha,beta weight parameters
      !     nneg - number negative x values
      !     meth - method, > 0 all use same method
      !     Order - number terms used for asymptotic approximations
      !     Ultra_Method = 0 t2 = pi/2 - t1 used, t1 used for interior calcs
      !  t2 is best for left end of shortcut polynomials, otherwise t1 is best
      ! ------------------------------------------
      nt = size(t1)
      ab(:,1) = (/a,b/)
      ab(:,2) = ab(:,1)
      rev = .false.
      iultra = 1        ! ultra0 or ultra2 using t1
      icase = JacobiType
      if(meth > 0)icase = 0
      select case(icase)
      case(RadauR,RadauL)
         ab(:,1) = (/b,a/) ! reversal for negatives
         rev(1) = .true.
      case(ShortGauss,ShortLobatto)
         iultra(1) = 0  ! ultra1 or ultra3 using t2
      end select
      select case(icase)
      case(0)
         ! how to call this now ??? need to fix for negatives
         call Jacobi_Pcalc(x,t1,t2,cs,sn,p,n,ab(:,1),meth,Order)
      case(Gauss,Lobatto)
         call Jacobi_Pcomposite(x,t1,t2,cs,sn,p,n,ab(:,1),reverse=.false.)
      case default   ! Radau and Short, nneg > 0
         i0 = 1; i1=nneg
         do k = 1,2
            Ultra_Method = iultra(k)
            call Jacobi_Pcomposite(x(i0:i1),t1(i0:i1),t2(i0:i1),cs(i0:i1),sn(i0:i1),p(i0:i1,:), &
               n,ab(:,k),reverse=rev(k))
            if(rev(k))then
               j = mod(n+1,2)   ! reverse sign p or p' if n odd or even
               p(i0:i1,j) = -p(i0:i1,j)
            endif
            i0 = i1+1; i1 = nt
         enddo
      end select
      pn = p(:,0:nd)
   End Function Jacobi_Pmx
   !----------------------------------------------------------------------------
   Subroutine Jacobi_Pcalc(x,t1,t2,cs,sn,p,n,ab,meth,Order)
      Integer, intent(in) :: n,meth
      Real(float), intent(in) :: x(:),t1(:),t2(:),cs(:),sn(:)
      Real(float), intent(out) :: p(:,0:)
      Real(float), intent(in) :: ab(2)
      Integer, optional :: Order
      logical :: Ultra
      Integer, parameter :: Mord(2) = (/12,4/)
      ! ------------------------------------------
      !  Calculates polynomial and derivative using different methods
      !     x,t - points
      !     n  - degree
      !     nd - number derivatives (0/1)
      !     a,b - Jacobi weight parameters
      !     meth - method > 0 asymptotic, otherwise recurrent
      !        Jacobi_Recurrent - recurrent method
      !        Jacobi_Asympt_Ultra - asymptotic methods ultraspherical
      !        Jacobi_Asympt_Radau - asymptotic methods Radau
      ! ------------------------------------------
      Ultra = (ab(1)==ab(2)) .or. (abs(nint(ab(2)*2.0_float))==1)
      if(Debug > 2)write(iout,'(a,2i6,i3,L2)')'# Jacobi_Pcalc: n,nx,meth,Ultra =',n,size(t1),meth,Ultra                      !
      if(meth <= 0)then
         p = Jacobi_Recurrent(x,n,1,ab(1),ab(2),Shortcut=.false.)
         p(:,1) = -p(:,1)*sin(ThetaFac*t1)  ! convert to theta derivative
      elseif(Ultra)then
         call Jacobi_Asympt_Ultra(t1,t2,cs,sn,p,n,ab,meth,Order)
      else
         call Jacobi_Asympt_Radau(t1,t2,cs,sn,p,n,ab,meth,Order)
      endif
   End Subroutine Jacobi_Pcalc
   !----------------------------------------------------------------------------
   Subroutine Jacobi_Pcomposite(x,t1,t2,cs,sn,p,n,ab,reverse)
      Real(float), intent(in) :: x(:),t1(:),t2(:),cs(:),sn(:)
      real(float), intent(out) :: p(:,0:)
      integer, intent(in) :: n
      real(float), intent(in) :: ab(2)
      logical :: reverse
      Integer :: nt,k,nn,nb
      Integer, parameter :: NuMM=7
      real(float) :: t0(NuMM-1),xm,zn,tr1,tx0,tx1
      Integer :: i0,i1,jM(0:NuMM),j0(NuMM),j1(NuMM)
      Integer, parameter :: MM(NuMM)  = (/4,5,7,9,12,0,4/)  ! method order (number terms)
      Integer, parameter :: mth(NuMM) = (/1,1,1,1,1,-1,2/)  ! -1 recurrent 1 interior, 2 boundary
      Real(float), parameter :: eps0 = epsilon(1.0_float)   ! epsf is safety factor
      Real(float), parameter :: epsf(NuMM-2) = (/0.1_float,0.1_float,0.5_float,1.0_float,2.0_float/)*eps0
      Real(float), parameter :: bexp(4) = (/1.0_float,1.0_float,0.9_float,1.0_float/)
      Real(float), parameter :: cexp(4) = (/1.507_float,0.789_float,2.527_float,0.900_float/)
      Real(float), parameter :: g0(4) = (/6.945_float,-13.80_float,58.44_float,0.000_float/)*fdegree
      Real(float), parameter :: g1(4) = (/8.671_float, 10.58_float, 2.50_float,9.621_float/)*fdegree
      real(float), Parameter :: bd(0:1) = (/2.26_float,0.11_float/)*fdegree
!     real(float), parameter :: br(0:1) = (/-1.2_float,fdegree*3700.0_float/)
!     tr0 = br(1)*zn**(br(0)) ! recurrent vs interior, close to t0 for M=12
      ! ------------------------------------------
      ! Calculate polynomial and derivative using combination of methods.
      !  x,t - points x & theta ordered increasing x (0 to 1), decreasing t
      !  n - polynomial degree
      !  a,b - Jacobi weight parameters
      ! The most efficient method for the desired accuracy is selected from:
      !  interior M=4,5,7,9,12, recurrent, boundary asymptotic M=4
      ! The method & order selected depends on the accuracy. Error for
      ! interior error correlated - bexp,cexp,g0 & g1, boundary - bd
      ! ------------------------------------------
      nt = size(x)
      nn = n   !;  tfac = 1.0_float
      nb = nint(ab(2)*2.0_float)     ! nb = 1 if b = +/-0.5
      if(abs(nb)==1)nn = 2*n + (nb+1)/2  ! nn is n for corresponding full polynomial
      zn = real(nn,float)
      tr1 = bd(0) + bd(1)*zn     ! boundary method break point
      i0 = 2*nint(ab(1))+1;    i1 = i0+1
      do k=1,NuMM-2
         xm = real(MM(k),float)
         tx0 = (g0(i0)+g1(i0)*xm)/(epsf(k)**(1.0_float/(xm+cexp(i0)))*zn**bexp(i0))
         tx1 = (g0(i1)+g1(i1)*xm)/(epsf(k)**(1.0_float/(xm+cexp(i1)))*zn**bexp(i1))
         t0(k) = max(tx0,tx1)
      enddo
      t0(NuMM-1) = min(tr1,t0(NuMM-2))   ! recurrent to fill in if needed
      if(reverse)then
         jM(0) = nt+1;  jM(NuMM) = 1
         do k=1,NuMM-1
            jM(k) = FirstGT(t1,t0(k))
         enddo
         j0 = jM(1:NuMM);     j1 = jM(0:NuMM-1) - 1
      else
         jM(0) = 1;  jM(NuMM) = nt+1
         do k=1,NuMM-1
            jM(k) = FirstLT(t1,t0(k))
         enddo
         j0 = jM(0:NuMM-1);   j1 = jM(1:NuMM) - 1
      endif
      if(Debug > 1)call PrintDebug()
      if(nn < Nasymptotic)then  ! case when Nasymptotic > Nboundary
         continue !* write code to eliminate interior asymptotic method
         stop 'Jacobi_Pcomposite require modification'
      endif
      do k = 1,NuMM
         i0 = jM(k-1);  i1 = jM(k)-1
         i0 = j0(k);  i1 = j1(k)
         if(i0 > i1)cycle
         call Jacobi_Pcalc(x(i0:i1),t1(i0:i1),t2(i0:i1),cs(i0:i1),sn(i0:i1),p(i0:i1,:), &
            n,ab,mth(k),Order=MM(k))
      enddo
      contains ! ------------
      subroutine PrintDebug()
         write(iout,'(a,L2)')'# Jacobi_Pcomposite: reverse?',reverse
         write(iout,'(a,8x,10f8.3)')'# Composite:',t0 ! cutoff angles
         write(iout,'(a,10i8)')'# 1st values :',j0
         write(iout,'(a,10i8)')'# total calcs:',j1-j0+1
      end subroutine PrintDebug
   End Subroutine Jacobi_Pcomposite
   !----------------------------------------------------------------------------
   Subroutine Jacobi_Asympt_Ultra(t1,t2,cs,sn,pp,n,ab,meth,Order)
      Real(float), intent(in) :: t1(:),t2(:),cs(:),sn(:)
      Real(float), intent(out) :: pp(:,0:)
      Integer, intent(in) :: n,meth
      Real(float), intent(in) :: ab(2)
      Integer, optional :: Order
      Real(float) :: ai(0:1)
      Integer :: nt,nn,nb,ord
      logical ::  short, odd
      Integer, parameter :: Mord(2) = (/12,4/)  ! default order
      ! ------------------------------------------
      !  Implements asymptotic formulas for ultraspherical polynomials
      !  Can be either the full polynomial or its shortcut equivalent
      !     tt - theta angle
      !     n - polynomial degree
      !     a,b - Jacobi Weight parameters
      !     meth, Order - method & order (see Jacobi_Asymptotic)
      ! ------------------------------------------
      nt = size(t1);    nn = n
      nb = nint(ab(2)*2.0_float)
      odd = nb==1;   short = abs(nb)==1
      ord = Mord(meth); if(present(Order))ord = Order
      if(short)then
         nn = 2*n + (nb+1)/2     ! ai = 1/an - inverse of scaling factor
         ai = Jacobi_Ends(n,ab(1),ab(2),Monic=.false.)/Jacobi_Ends(nn,ab(1),ab(1),Monic=.false.)
      endif                                                            !
      if(Debug > 1)write(iout,'(a,3i6,3i3,4L2)') &
         '# _Asympt_Ultra: nt,n,nn,meth,ord,Ultra_method,short,odd =', &
          nt,n,nn,meth,ord,Ultra_method,short,odd
      call Jacobi_Asymptotic(t1,t2,cs,sn,nn,ab(1),meth,ord,pp)
      if(short)then
         pp(:,0) = ai(1)*pp(:,0)
         pp(:,1) = pp(:,1)*(ai(1)*0.5_float)
         if(odd)then
            pp(:,0) = pp(:,0)/cs    !
            pp(:,1) = (pp(:,1) + (0.5_float)*sn*pp(:,0))/cs
         endif
      endif
   End Subroutine Jacobi_Asympt_Ultra
   !----------------------------------------------------------------------------
   Subroutine Jacobi_Asympt_Radau(t1,t2,cs,sn,pp,n,ab,meth,Order)
      Real(float), intent(in) :: t1(:),t2(:),cs(:),sn(:)
      Real(float), intent(out) :: pp(:,0:)
      Integer, intent(in) :: n,meth
      Real(float), intent(in) :: ab(2)
      Integer, optional :: Order
      Real(float) :: zn,z1,zx
      Integer :: nt,ord
      logical :: RadauL
      Integer, parameter :: Mord(2) = (/12,4/)
      real(float), parameter :: a0 = 0.0_float
      ! ------------------------------------------
      !  Implements asymptotic formulas for Radau cases
      !     tt - theta angle
      !     n - polynomial number
      !     a,b - Jacobi Weight parameters
      !     meth, Order - method & order (see Jacobi_Asympt)
      ! ------------------------------------------
      nt = size(t1)
      RadauL = ab(1)==0.0_float .and. ab(2)==1.0_float
      ord = Mord(meth); if(present(Order))ord = Order                  !
      if(Debug > 1)write(iout,'(a,2i6,3i3,4L2)') &
         '# _Asympt_Radau: n,nt,meth,ord,Ultra_method,RadauL =', &
         n,nt,meth,ord,Ultra_method,RadauL
      call Jacobi_Asymptotic(t1,t2,cs,sn,n,a0,meth,ord,pp)
      zn = Real(n,float); z1 = (1.0_float)/Real(n+1,float)
      zx = 2.0*z1*(zn+1.0_float)
      if(RadauL)then  ! convert Legendre to Radau
         pp(:,0) = pp(:,0) + z1*xminus(t1)*pp(:,1)/sn
         pp(:,1) = (zx*pp(:,1) - zn*sn*pp(:,0))/xplus(t1)
      else  ! RadauR
         pp(:,0) = pp(:,0) - z1*xplus(t1)*pp(:,1)/sn
         pp(:,1) = (zx*pp(:,1) + zn*sn*pp(:,0))/xminus(t1)
      endif
      contains ! ---- these should inline -------
      pure elemental function xplus(t) result(xp)
         real(float), intent(in) :: t
         real(float) :: xp
         ! calculate 1 + x in theta coordinates
         xp = (2.0_float)*cos(t*0.5_float)**2
      end function xplus
      ! -----------------------------
      pure elemental function xminus(t) result(xm)
         real(float), intent(in) :: t
         real(float) :: xm
         ! calculate 1 - x in theta coordinates
         xm = (2.0_float)*sin(t*0.5_float)**2
      end function xminus
   End Subroutine Jacobi_Asympt_Radau
   !----------------------------------------------------------------------------
   subroutine Jacobi_Asymptotic(t1,t2,cs,sn,nn,a,meth,ord,p)
      Real(float), intent(in) :: t1(:),t2(:),cs(:),sn(:),a
      Integer, intent(in) :: nn
      Real(float), intent(out) :: p(:,0:)
      Integer :: meth,n,nt,ord,M,ia
      Real(float) :: zn,zn2,a1,an
      ! ------------------------------------------
      ! Evaluates asymptotic formulas for Jacobi polynomials
      !  t1 - theta angle
      !  t2 = pi/2 - t1
      !  cs - cos(t1)
      !  sn = sin(t1)
      !  n - polynomial number
      !  a - Jacobi alpha
      !  meth - method number:
      !     1. Jacobi_AsymUltraX interior routines
      !     2. Jacobi_AsymB - near boundary expansion from Bogaert(2014)
      !  Order - order = 1+
      !  Jacobi_AsymUltra0 - uses t2, no derivative calc
      !  Jacobi_AsymUltra1 - uses t1, no derivative calc
      !  Jacobi_AsymUltra2 - uses t2, w/ derivative calc
      !  Jacobi_AsymUltra3 - uses t1, w/ derivative calc
      !  Deriv_meth controls derivative calculation:
      !  dP_n^(alpha) = 0.5*(n + 2*alpha + 1)P_n-1^(alpha+1) Deriv_meth.ne.0 (not used)
      ! ------------------------------------------
      nt = size(t1)
      if(nt == 0)return
      n = nn;  M = max(1,ord)
      ia = nint(a)
      if(Debug > 0)write(iout,'(a,i7,10i3)')'Jacobi_Asymptotic: n,a,methods =,', &
         nn,ia,meth,Ultra_method
      if(meth == 1 .and. Deriv_meth == 0)then
         if(Ultra_method == 0)then
            p = Jacobi_AsymUltra2(t2,cs,sn,a,n,M)   ! interior, t2 = pi/2 - theta
         else
            p = Jacobi_AsymUltra3(t1,cs,sn,a,n,M)   ! interior, t1 = theta
         endif
      elseif(meth == 1)then   ! Deriv_meth .ne. 0
         an = (-0.5_float)*real(n+1,float) - a
         a1 = a + 1.0_float
         if(Ultra_method == 0)then
            p(:,0) = Jacobi_AsymUltra0(t2,sn,a,n,M)
            p(:,1) = an*sn*Jacobi_AsymUltra0(t2,sn,a1,n-1,M)
         else
            p(:,0) = Jacobi_AsymUltra1(t1,sn,a,n,M)
            p(:,1) = an*sn*Jacobi_AsymUltra1(t1,sn,a1,n-1,M)
         endif
      else
         n = nn + ia   ! use n+1 for a = 1
         M = min(M,4) - 1
         call Jacobi_AsymB(t1,cs,sn,n,M,p)       ! boundary method
         if(ia == 1)then
         ! convert Legendre P_n+1 and dP_n+1 to Jacobi(1,1) P_n and dP_n
            zn  = Real(n,float);
            zn2 = (-2.0_float)/Real(n+1,float)
            call swap(p(:,0),p(:,1))   ! function avoids scratch array
            p(:,0) = zn2*p(:,0)/sn
            p(:,1) = (zn*p(:,1) - cs*p(:,0))*2.0_float/sn
         endif
      endif
      contains ! ----------------
      subroutine swap(a,b)
         real(float) :: a(:),b(:),tmp
         integer :: i
         do i = 1,size(a)
            tmp = a(i)
            a(i) = b(i)
            b(i) = tmp
         enddo
      end subroutine swap
   End Subroutine Jacobi_Asymptotic
   !----------------------------------------------------------------------------
   Function Jacobi_AsymUltra0(t,sn,a,n,M)  result(p)
      Real(float), intent(in) :: t(:),sn(:),a
      Integer, intent(in) :: n,M
      Real(float) :: p(size(t))
      Real(float) :: ap,am,bn,hm(0:M)
      Real(float) :: hn(size(t))
      Integer :: k,ia,ncase
      ! ------------------------------
      ! t = pi/2 - theta (formula simplified accordingly)
      ! sn = sin(theta)
      ! a - alpha
      ! Szego gives ultraspherical formula, p. 94, Eq. (4.9.25)
      ! Also at dlmf.nist.gov/18.15. The two seem to disagree in the constant
      ! term when a .ne. 0. The dlmf version is below
      ! dlmf also gives general Jacobi case, requires double summation, will avoid
      ! Can convert for Radau case using Eq. (2.80) in Sec. 2.4.4
      !  cos(n*pi/2-b) = cos(n*pi/2)cos(b) + sin(n*pi/2)sin(b)
      !  sin(n*pi/2-b) = sin(n*pi/2)cos(b) - cos(n*pi/2)sin(b)
      ! ------------------------------                                 !
      ia = nint(a)
      ap = a + 0.5_float
      am = a - 0.5_float
      bn = ap + real(n,float)
      hm(0) = GammaRatio(a+1.0_float,a+1.5_float,n)*real(2**ia,float)/sqrt(pi2)
      hm(1:M) = (/(Real(k,float),k=1,M)/)         ! 1,2,3,4,...
      hm(1:M) = (0.5_float)*(hm(1:M)+am)*(hm(1:M)-ap)/(hm(1:M)*(hm(1:M)+bn))
      hn = sqrt(sn)/(sn**ia)
      p = 0.0_float
      ncase = mod(n,4)
      do k=0,M
         bn = real(n+k,float) + ap
         hn = hm(k)*hn/sn
         select case(ncase)
         case(0)  ! cos(n*pi/2) =+1, sin(n*pi/2)= 0
            p = p + hn*cos(bn*t)
         case(1)  ! cos(n*pi/2) = 0, sin(n*pi/2)=+1
            p = p + hn*sin(bn*t)
         case(2)  ! cos(n*pi/2) =-1, sin(n*pi/2)= 0
            p = p - hn*cos(bn*t)
         case(3)  ! cos(n*pi/2) = 0, sin(n*pi/2)=-1
            p = p - hn*sin(bn*t)
         end select
         if(Debug > +4)then
            ia = min(size(t),15)
            write(iout,'(a,a1,i2,a1,f5.1)')'#_AsymUltra0: k,bn =',tb,k,tb,bn
            if(k == 0)call vectorprint('hm:',hm,fmtx=fmt22)
            if(n < 179 .and. k < 4) &
            call arrayprint('Jacobi_AsymUltra0: t,hn,p',t(:ia),hn(:ia),p(:ia),fmtx=fmt22)
         endif
      enddo
   End Function Jacobi_AsymUltra0
   !----------------------------------------------------------------------------
   Function Jacobi_AsymUltra1(t,sn,a,n,M)  result(p)
      Real(float), intent(in) :: t(:),sn(:),a
      Integer, intent(in) :: n,M
      Real(float) :: p(size(t))
      Real(float) :: ap,am,bn,bk,hm(0:M)
      Real(float) :: hn(size(t))
      Integer :: k,ia
      character :: chk*3
      ! ------------------------------
      !  t - theta
      !  sn - sin(theta)
      !  a - polynomial alpha parameter
      ! Szego gives ultraspherical formula, p. 94, Eq. (4.9.25)
      ! Also at dlmf.nist.gov/18.15. The two seem to disagree in the constant
      ! term when a .ne. 0. The dlmf version is below
      ! dlmf also gives general Jacobi case, requires double summation, will avoid
      ! Can convert for Radau case using Eq. (2.80) in Sec. 2.4.4
      ! ------------------------------
      ia = nint(a)
      ap = a + 0.5_float
      am = a - 0.5_float
      bn = ap + real(n,float)
      hm(0) = GammaRatio(a+1.0_float,a+1.5_float,n)*real(2**ia,float)/sqrt(pi2)      ! dlmf.nist.gov
      hm(1:M) = (/(Real(k,float),k=1,M)/)         ! 1,2,3,4,...
      hm(1:M) = (0.5_float)*(hm(1:M)+am)*(hm(1:M)-ap)/(hm(1:M)*(hm(1:M)+bn))
      hn = sqrt(sn)/(sn**ia)
      p = 0.0_float
      do k=0,M
         bn = real(n+k,float) + ap
         bk = real(k,float) + ap
         hn = hm(k)*hn/sn
         p = p + hn*cos(bn*t - bk*pi2)
         write(chk,'(i3)')k
         if(Debug > +4)then
            ia = min(size(t),15)
            write(iout,'(a,a1,f5.1,3(a1,g23.16))')'#_AsymUltra1: bk,cos(bk),sin(bk) =', &
               tb,bk,tb,cos(bk*pi2),tb,sin(bk*pi2)
            if(k == 0)call vectorprint('hm:',hm,fmtx=fmt22)
            if(n < 179 .and. k < 4) &
            call arrayprint('Jacobi_AsymUltra1: t,hn,p'//chk,t(:ia),hn(:ia),p(:ia),fmtx=fmt22)
         endif
      enddo
   End Function Jacobi_AsymUltra1
   !----------------------------------------------------------------------------
   Function Jacobi_AsymUltra2(t,cs,sn,a,n,M)  result(p)
      Real(float), intent(in) :: t(:),cs(:),sn(:),a
      Integer, intent(in) :: n,M
      Real(float) :: p(size(t),0:1)
      Real(float) :: ap,am,bn,bk,hm(0:M)
      Real(float) :: hn(size(t)),ca(size(t)),sa(size(t))
      Integer :: k,ia,ncase
      character :: chk*3
      ! ------------------------------
      ! Szego gives ultraspherical formula, p. 94, Eq. (4.9.25)
      ! Also at dlmf.nist.gov/18.15. The two seem to disagree in the constant
      ! term when a .ne. 0. The dlmf version is below.
      !  t = pi/2 - theta
      ! code uses expansion, such that
      !  cos(n*pi/2-b) = cos(n*pi/2)cos(b) + sin(n*pi/2)sin(b)
      !  sin(n*pi/2-b) = sin(n*pi/2)cos(b) - cos(n*pi/2)sin(b)
      ! values depends on mod(n,4), either cos(n*pi/2) or sin(n*pi/2) is zero
      ! ------------------------------                                 !
      ia = nint(a)
      ncase = mod(n,4)
      ap = a + 0.5_float
      am = a - 0.5_float
      bn = ap + real(n,float)
      hm(0) = GammaRatio(a+1.0_float,a+1.5_float,n)*real(2**ia,float)/sqrt(pi2)      ! dlmf.nist.gov
      hm(1:M) = (/(Real(k,float),k=1,M)/)         ! 1,2,3,4,...
      hm(1:M) = (0.5_float)*(hm(1:M)+am)*(hm(1:M)-ap)/(hm(1:M)*(hm(1:M)+bn))
      hn = (1.0_float)/(sqrt(sn)*sn**ia)
      p = 0.0_float
      do k=0,M
         bn = real(n+k,float) + ap
         bk = real(k,float) + ap
         select case(ncase)
         case(0)  ! cos(n*pi/2) =+1, sin(n*pi/2)= 0
            ca = cos(bn*t);  sa = -sin(bn*t)
         case(1)  ! cos(n*pi/2) = 0, sin(n*pi/2)=+1
            ca = sin(bn*t);  sa = cos(bn*t)
         case(2)  ! cos(n*pi/2) =-1, sin(n*pi/2)= 0
            ca = -cos(bn*t); sa = sin(bn*t)
         case(3)  ! cos(n*pi/2) = 0, sin(n*pi/2)=-1
            ca = -sin(bn*t); sa = -cos(bn*t)
         end select
         hn = hm(k)*hn
         p(:,0) = p(:,0) + hn*ca
         p(:,1) = p(:,1) - hn*sa*bn
         hn = hn/sn
         p(:,1) = p(:,1) - hn*ca*cs*bk
         if(Debug > +4)then
            write(iout,'(a,a1,i5,a1,i2,4(a1,f5.1))')'#_AsymUltra2: n,k,a,bk =',tb,n,tb,k,tb,a,tb,bk
            if(k == 0)call vectorprint('hm:',hm,fmtx=fmt22)
            write(chk,'(i3)')k
            ia = min(size(t),15)
            if(n < 179 .and. k < 4)call arrayprint('Jacobi_AsymUltra2: t,hn,p'//chk, &
               t(:ia),hn(:ia),p(:ia,:),fmtx=fmt22)
         endif
      enddo
   End Function Jacobi_AsymUltra2
   !----------------------------------------------------------------------------
   Function Jacobi_AsymUltra3(t,cs,sn,a,n,M)  result(p)
      Real(float), intent(in) :: t(:),cs(:),sn(:),a
      Integer, intent(in) :: n,M
      Real(float) :: p(size(t),0:1)
      Real(float) :: am,ap,bk,bn,hm(0:M)
      Real(float) :: hn(size(t)),ca(size(t)),sa(size(t)) !ak(size(t)),
      real(float), parameter :: x707(0:1) = (/-1.0_float/sqrt(2.0_float),1.0_float/sqrt(2.0_float)/)
      Integer :: k,ia,is,ic
      character :: chk*3
      ! ------------------------------
      ! Szego gives ultraspherical formula, p. 94, Eq. (4.9.25)
      ! Also at dlmf.nist.gov/18.15. The two seem to disagree in the constant
      ! term when a .ne. 0. The dlmf version is below
      ! dlmf also gives general Jacobi case, requires double summation, will avoid
      ! Can convert for Radau case using Eq. (2.80) in Sec. 2.4.4
      !   ca = cos(bn*t)*cos(bk*pi2) + sin(bn*t)*sin(bk*pi2)
      !   sa = sin(bn*t)*cos(bk*pi2) - cos(bn*t)*sin(bk*pi2)
      ! ------------------------------
      ia = nint(a)
      ap = a + 0.5_float
      am = a - 0.5_float
      bn = ap + real(n,float)
      hm(0) = GammaRatio(a+1.0_float,a+1.5_float,n)*real(2**ia,float)/sqrt(pi2)      ! dlmf.nist.gov
      hm(1:M) = (/(Real(k,float),k=1,M)/)         ! 1,2,3,4,...
      hm(1:M) = (0.5_float)*(hm(1:M)+am)*(hm(1:M)-ap)/(hm(1:M)*(hm(1:M)+bn))
      hn = (1.0_float)/(sqrt(sn)*sn**ia)
      p = 0.0_float
      do k=0,M
         bn = real(n+k,float) + ap
         bk = real(k,float) + ap
         hn = hm(k)*hn
         ic = mod(nint((bk+2.0)*0.5),2)
         is = mod(nint((bk+1.0)*0.5),2)
         ca = cos(bn*t)*x707(ic) + sin(bn*t)*x707(is)
         sa = sin(bn*t)*x707(ic) - cos(bn*t)*x707(is)
         p(:,0) = p(:,0) + hn*ca
         hn = hn/sn
         p(:,1) = p(:,1) - hn*(bk*ca*cs + bn*sn*sa)
         if(Debug > +4)then    ! ---------
            write(chk,'(i3)')k
            ia = min(size(t),15)
            write(iout,'(a,a1,f5.1,4(a1,f7.4))')'#_AsymUltra3: bk,cos(bk),sin(bk) =', &
               tb,bk,tb,x707(ic),tb,cos(bk*pi2),tb,x707(is),tb,sin(bk*pi2)
            if(k == 0)call vectorprint('hm:',hm,fmtx=fmt22)
            if(n < 179 .and. k < 4)call arrayprint('Jacobi_AsymUltra3: t,hn,p'//chk, &
               t(:ia),hn(:ia),p(:ia,:),fmtx=fmt22)
         endif
      enddo
   End Function Jacobi_AsymUltra3
   !----------------------------------------------------------------------------
   Subroutine Jacobi_AsymB(t,cs,sn,n,M,pn)
      Real(float), intent(in) :: t(:),cs(:),sn(:)
      Integer, intent(in) :: n,M
      Real(float), intent(out) :: pn(:,0:)
      Real(float) :: vn,vi
      Real(float) :: av(size(t),0:3,0:1),bj(size(t),0:3,0:1)
      integer :: i,j
      ! ------------------------------
      ! Legendre polynomials via asymptotic expansion in Bessell functions from:
      ! ITERATION-FREE COMPUTATION OF GAUSS–LEGENDRE QUADRATURE NODES AND WEIGHTS
      ! Bogaert SIAM J. SCI. COMPUT. (2014) - Eq. (2.4)
      !  t - theta, 0 < t < pi/2, or 0 < x < 1
      !  n - polynomial no.
      !  M - order, M+1 terms included
      !  cs - cos(theta)
      !  sn -sin(theta)
      !  vn - 1/(n+1/2)
      !  av(:,k) - coefficient k: av(0)*{bj(0) + av(1)*bj(1)*vn + av(2)*bj(2)*vn**2 ...}
      !  bj(:,k) - bessell functions, first kind order k
      ! ------------------------------                                 !
      !write(iout,'(a,2i15)')'__AsymB...................loc(t), loc(p) =',loc(t),loc(pn(1,0))
      vi = (0.5_float)*real(2*n+1,float)
      vn = 1.0_float/vi
      av = Jacobi_AsymB_Bfunc(t,cs,sn)    ! calculates all 4
      bj = Jacobi_AsymB_Bessell(t*vi)     ! calculates all 4
      do i=1,M  ! derivative of coefficients
         av(:,i,1) = av(:,0,1)*av(:,i,0) + av(:,0,0)*av(:,i,1)
      enddo
      Select case(M)
      case(0)
         pn(:,0) = av(:,0,0)*bj(:,0,0)
         pn(:,1) = av(:,0,1)*bj(:,0,0) + av(:,0,0)*bj(:,0,1)*vi
      case(1)
         pn(:,0) = av(:,0,0)*(bj(:,0,0)+vn*av(:,1,0)*bj(:,1,0))
         pn(:,1) = av(:,0,1)*bj(:,0,0) +vn*av(:,1,1)*bj(:,1,0) &
                 + av(:,0,0)*(bj(:,0,1)*vi+av(:,1,0)*bj(:,1,1))
      case(2)
         pn(:,0) = av(:,0,0)*(bj(:,0,0)+vn*(av(:,1,0)*bj(:,1,0)+vn*av(:,2,0)*bj(:,2,0)))
         pn(:,1) = av(:,0,1)* bj(:,0,0)+vn*(av(:,1,1)*bj(:,1,0)+vn*av(:,2,1)*bj(:,2,0))  &
                 + av(:,0,0)*(bj(:,0,1)*vi +av(:,1,0)*bj(:,1,1)+vn*av(:,2,0)*bj(:,2,1))
      case default   ! case 3
         pn(:,0) = av(:,0,0)*(bj(:,0,0)+vn*(av(:,1,0)*bj(:,1,0)+vn*(av(:,2,0)*bj(:,2,0)+vn*av(:,3,0)*bj(:,3,0))))
         pn(:,1) = av(:,0,1)* bj(:,0,0)+vn*(av(:,1,1)*bj(:,1,0)+vn*(av(:,2,1)*bj(:,2,0)+vn*av(:,3,1)*bj(:,3,0)))  &
                 + av(:,0,0)*(bj(:,0,1)*vi +av(:,1,0)*bj(:,1,1)+vn*(av(:,2,0)*bj(:,2,1)+vn*av(:,3,0)*bj(:,3,1)))
      End Select
      if(Debug > +4)then
         j = max(1,size(t)-10)
         write(iout,'(a,f10.3)')' vi =',vi
         call arrayprint('Bfunc:',t(j:),av(j:,:,0),av(j:,:,1),fmtx=fmt22)
         call arrayprint('Bessl:',t(j:),bj(j:,:,0),bj(j:,:,1),fmtx=fmt22)
         call ArrayPrint('Jacobi_AsymB',t(j:),pn(j:,:),fmtx=fmt22)
      endif
   End Subroutine Jacobi_AsymB
   ! -------------------------------------
   Function Jacobi_AsymB_Bessell(x) result(bs)
      real(float), intent(in) :: x(:)
      integer :: j
      real(float) :: bs(size(x),0:3,0:1),xi(size(x)),x2(size(x))
      ! --------------------------------
      ! Bessell functions needed for Legendre polynomial approximation
      !  x - values will normally be near zeros of J0 (Gauss) or J1 (Lobatto)
      !  *** Could linearize about the zeroes to speedup ***
      ! Could avoid calls to Bessell functions using a procedure similar to
      ! Bogaert, i.e. combine Mahon's formula for the zeroes with the asymptotic
      ! representation of Bessell functions in terms of sin and cos.
      ! --------------------------------
      xi = (2.0_float)/x
      ! values
      bs(:,0,0) = BESSEL_J0(x)
      bs(:,1,0) = BESSEL_J1(x)
      bs(:,2,0) = bs(:,1,0)*xi - bs(:,0,0)
      x2 = xi*2.0_float
      bs(:,3,0) = bs(:,1,0)*(xi*x2 - 1.0_float) - bs(:,0,0)*x2
      ! derivatives
      bs(:,0,1) = -bs(:,1,0)
      bs(:,1,1) = bs(:,0,0) - bs(:,1,0)*xi*0.5_float
      bs(:,2,1) = bs(:,0,0)*xi + bs(:,1,0)*(1.0_float - xi*xi)
      x2 = xi*xi*3.0_float
      bs(:,3,1) = bs(:,1,0)*xi*(2.5_float-x2) - bs(:,0,0)*(1.0_float-x2)
      j = max(1,size(x)-5)
   End Function Jacobi_AsymB_Bessell
   !----------------------------------------------------------------------------
   Function Jacobi_AsymB_Bfunc(t,cs,sn)  result(b)
      real(float), intent(in) :: t(:),cs(:),sn(:)
      integer :: nt,k
      real(float) :: b(size(t),0:3,0:1),si(size(t)),t2(size(t)),s2(size(t)),c2(size(t))
      real(float), parameter :: tsxx(6) = (/0.09_float,0.17_float,0.25_float,0.07_float,0.17_float,0.26_float/)
      Real(float),parameter :: tsx(3,0:1) = Reshape(tsxx,(/3,2/))*(epsilon(t)/2.2e-16)**(1/7.5)
      integer, parameter :: nbt = 6
      real(float), parameter :: b0p(nbt)=(/ &
         1.0_float,3.0_float/20.0_float,61.0_float/3360.0_float,1261.0_float/604800.0_float, &
         79.0_float/337920.0_float,66643.0_float/2583060480.0_float/)*(1.0_float/6.0_float)
      real(float), parameter :: b1t(nbt)=(/ &
         1.0_float,1.0_float/15.0_float,2.0_float/315.0_float,1.0_float/1575.0_float, &
         2.0_float/31185.0_float,1382.0_float/212837625.0_float/)/(-24.0_float)
      real(float), parameter :: b1p(nbt)=(/ &
         1.0_float,1.0_float/5.0_float,2.0_float/63.0_float,1.0_float/225.0_float, &
         2.0_float/3465.0_float,1382.0_float/19348875.0_float/)/(-24.0_float)
      real(float), parameter :: b2t(nbt)=(/ &
         7.0_float,26.0_float/21.0_float,19.0_float/105.0_float,50.0_float/2079.0_float, &
         42842.0_float/14189175.0_float,148.0_float/405405.0_float/)/(1920.0_float)
      real(float), parameter :: b2p(nbt)=(/ &
         14.0_float,104.0_float/21.0_float,38.0_float/35.0_float,400.0_float/2079.0_float, &
         85684.0_float/2837835.0_float,592.0_float/135135.0_float/)/(1920.0_float)
      real(float), parameter :: b3t(nbt)=(/31.0_float/3.0_float, &
         157.0_float/45.0_float,1093.0_float/1485.0_float,257456.0_float/2027025.0_float, &
         118982.0_float/6081075.0_float,111373.0_float/39760875.0_float/)/(-21504.0_float)
      real(float), parameter :: b3p(nbt)=(/31.0_float,  &
         157.0_float/9.0_float,7651.0_float/1485.0_float,257456.0_float/225225.0_float, &
         118982.0_float/552825.0_float,1447849.0_float/39760875.0_float/)/(-21504.0_float)
      ! ------------------------------------
      ! Functions given by Bogaert Eq. (2.9), derivatives from Maple
      ! most of these are subject to roundoff near t = 0, so 6 term Taylor
      ! series is used (probably not necessary, but done anyway)
      !  b0 = sqrt(t/sin(t)) (this is the multiplier out front)
      !  b1 = (t*cos(t) - sin(t))/(8*t*sin(t))
      !  db0/dt = -4*b1*b0
      !  ... others given in text
      ! ------------------------------------
      nt = size(t);    t2 = t*t
      si = (0.125_float)/(t*sn)
      s2 = sn*sn;      c2 =cs*cs
      k = 1;  b = 0.0_float
      where(t > tsx(1,0))  ! 1st term
         b(:,1,0) = (t*cs-sn)*si
      elsewhere
         b(:,1,0) = (b1t(1)+t2*(b1t(2)+t2*(b1t(3)+t2*(b1t(4)+t2*b1t(5)+t2*b1t(6)))))*t
      endwhere
      b(:,0,0) = sqrt(t/sn) ! 0th term
      b(:,0,1) = (-4.0_float)*b(:,1,0)*b(:,0,0)
      where(t > tsx(1,1))  ! 1st term derivative
         b(:,1,1) = (8.0_float)*(s2 - t2)*si*si
      elsewhere
         b(:,1,1) = (b1p(1)+t2*(b1p(2)+t2*(b1p(3)+t2*(b1p(4)+t2*b1p(5)+t2*b1p(6)))))
      endwhere
      where(t > tsx(2,0))  ! 2nd term
         b(:,2,0) = ((6.0_float)*t*cs*sn+(15.0_float+t2)*cs*cs+(8.0_float)*t2-15.0_float) &
                     *(0.5_float)*si*si
      elsewhere
         b(:,2,0) = (b2t(1)+t2*(b2t(2)+t2*(b2t(3)+t2*(b2t(4)+t2*b2t(5)+t2*b2t(6)))))*t2
      endwhere
      where(t > tsx(2,1))  ! 2nd term derivative
         b(:,2,1) = (cs*(5.0_float*cs*sn-t*cs*cs+3.0_float*t*t2+t) + sn*(t2-5.0_float)) &
                     *(-24.0_float)*si**3
      elsewhere
         b(:,2,1) = (b2p(1)+t2*(b2p(2)+t2*(b2p(3)+t2*(b2p(4)+t2*b2p(5)+t2*b2p(6)))))*t
      endwhere
      where(t > tsx(3,0))  ! 3rd term
         b(:,3,0) = ((16.0_float*t2+21.0_float - (t2+21.0_float)*c2)*t*cs  &
                + ((3.0_float*t2 + 63.0_float)*c2 + 24.0_float*t2 - 63.0_float)*sn)  &
                  *(2.5_float)*(si**3)
      elsewhere
         b(:,3,0) = (b3t(1)+t2*(b3t(2)+t2*(b3t(3)+t2*(b3t(4)+t2*b3t(5)+t2*b3t(6)))))*t*t2
      endwhere
      where(t > tsx(3,1))  ! 3rd term derivative
         b(:,3,1)=(((3.0_float*t2 + 189.0_float)*c2 + (42.0_float-29.0_float*t2)*t2 - 378.0_float)*c2 &
            + (42.0_float*c2 - 54.0_float*t2 - 42.0_float)*t*sn*cs - (45.0_float+16.0_float*t2)*t2 &
            + 189.0_float)*(20.0_float)*(si**4)
      elsewhere
         b(:,3,1) = (b3p(1)+t2*(b3p(2)+t2*(b3p(3)+t2*(b3p(4)+t2*b3p(5)+t2*b3p(6)))))*t2
      endwhere
      k = max(1,nt-50)
   End Function Jacobi_AsymB_Bfunc
   !----------------------------------------------------------------------------
   Function Jacobi_Recurrent(x,n,nd,alpha,beta,Shortcut) result(p)
      real(float), intent(in) :: x(:)
      integer, intent(in) :: n,nd
      real(float), optional :: alpha,beta
      logical, optional :: Shortcut
      real(float) :: p(size(x),n:n+nd),xx(size(x))
      real(float) :: abx(2),a,b,bb,an(0:1),fac4
      integer :: n2
      logical :: shrt,shft
      ! --------------------------------------
      ! uses recurrent relationships to calculate Jacobi Polynomials n0 to n
      ! together with nd derivatives
      !  x - x values
      !  t - theta values must be folded to use boundary asymptotic
      !  n0,n - 1st, last polynomial return
      !  nd - number of derivatives (only 1st calculated here)
      !  alpha,beta - Jacobi weight parameters
      !  meth = 0 forces full recurrent, i.e. no shortcut, no boundary asymptotic
      !  short - .true. for shortcut calcs
      !  xx = 2*x^2 - 1 or x^2 if shifted
      ! --------------------------------------
      abx = abparm(alpha,beta); a = abx(1);  b = abx(2)
      shrt = ShortOK;    if(present(Shortcut))shrt = Shortcut
      shft = ShiftOK;
      shrt = shrt .and. (a == b) .and. (n > 1)
      if(Debug > 0) &
         write(iout,'(a,2i6,2f5.1,4L2)')'# Jac_Recurrent: nx,n,a,b,short,shft=',size(x),n,a,b,shrt,shft
      if(.not.shrt)then
         p = Jacobi_Rec(x,n,n,nd,a,b)
      else  ! shortcut calcs
         n2 = n/2
         bb = real(2*mod(n,2)-1,float)*0.5_float
         fac4 = 2.0_float
         xx = x*x
         if(.not.shft)then
            fac4 = 4.0_float
            xx = (2.0_float)*xx - 1.0_float
         endif
         p = Jacobi_Rec(xx,n2,n2,nd,a,bb,shft)
         an = Jacobi_Ends(n,a,a,Monic=.false.)/Jacobi_Ends(n2,a,bb,Monic=.false.)
         p = an(1)*p
         if(mod(n,2)==1)then  ! odd
            if(nd > 0) &
               p(:,n+1) = p(:,n) + fac4*x*x*p(:,n+1)
            p(:,n) = x*p(:,n)
         elseif(nd > 0)then
            p(:,n+1) = fac4*x*p(:,n+1)
         endif
      endif
   End Function Jacobi_Recurrent
   !----------------------------------------------------------------------------
   Function Jacobi_Rec(x,n0,n,nd,a,b,shift) result(p)
      real(float), intent(in) :: x(:),a,b
      integer, intent(in) :: n0,n,nd
      logical, optional :: shift
      real(float) :: p(size(x),n0:n+nd),pr(size(x),3)
      real(float) :: ar(0:n,4),c(3)
      logical :: shft
      real(float), Parameter :: epsx = epsilon(x)*1.0_float
      ! --------------------------------------
      ! uses recurrence to calculate Jacobi Polynomials n0 to n and
      ! 1st derivative when nd > 0
      ! --------------------------------------
      shft = Default_shift;  if(present(shift))shft = shift
      ar = x_jacobi(n+1,a,b,Monic=.false.,Shift=shift)   ! recurrence coefficients
      pr = Jacobi_RecP3(x,n,ar,Monic=.false.)
      p(:,n0:n) = pr(:,3-n+n0:3)  ! save those desired (only valid for n-2 to n
      if(nd > 0)then    ! calculate 1st derivative
         c = Jacobi_Deriv1(n,a,b,Monic=.false.,Shift=shift)
         p(:,n+1) = (c(1)*x+c(2))*pr(:,3) + c(3)*pr(:,2)
         if(.not.shft)then
            p(:,n+1) = p(:,n+1)/max(epsx,(1.0_float-x*x))
         else
            p(:,n+1) = p(:,n+1)/max(epsx,(x*(1.0_float-x)))
         endif
      endif
   End Function Jacobi_Rec
   !----------------------------------------------------------------------------
   Function Jacobi_RecP3(x,n,ab,Monic)  result(p)
      real(float), intent(in) :: x(:),ab(0:,:)
      integer, intent(in) :: n
      Logical, optional :: Monic
      real(float) :: p(size(x),3)
      integer :: k,k0,k1,k2,k3
      integer, parameter :: kx(0:2,0:2) = reshape((/3,1,2,2,3,1,1,2,3/),(/3,3/))
      logical :: monc
      ! ------------------------------------------
      ! Uses recurrent relations to calculate Jacobi polynomials n-2 to n
      !  x  - x where polynomials are calculated
      !  n  - highest order polynomial calculated
      !  ab - recurrence coef., etc.
      !  Monic - .true. if monic polynomial
      ! Monic form gives
      !  Pn+1 = (x - ab(n,1))*Pn - ab(n,2)*Pn-1
      !  ab(n,3) - leading coefficient of conventional form
      ! Conventional form (Monic = false) give coefficients:
      !  Pn+1 = (x*ab(n,3) - ab(n,1))*Pn - ab(n,2)*Pn-1
      ! In either case:
      !  ab(n,4) - Integral(Pn^2*w(x))
      ! ------------------------------------------
      monc = Default_monic;    if(present(Monic))monc = Monic
      k = mod(n,3)
      k0 = kx(0,k); k1 = kx(1,k); k2 = kx(2,k)
      p(:,k0) = 1.0_float
      if(.not.monc)then
         p(:,k1) = x*ab(0,3) - ab(0,1)
         do k = 1,n-1
            p(:,k2) = (x*ab(k,3) - ab(k,1))*p(:,k1) - ab(k,2)*p(:,k0)
            k3 = k0; k0 = k1; k1 = k2; k2 = k3
         enddo
      else
         p(:,k1) = x - ab(0,1)
         do k = 1,n-1
            p(:,k2) = (x - ab(k,1))*p(:,k1) - ab(k,2)*p(:,k0)
            k3 = k0; k0 = k1; k1 = k2; k2 = k3
         enddo
      endif
   End Function Jacobi_RecP3
   !----------------------------------------------------------------------------
   Function Jacobi_All(x,n,ab,Monic)  result(p)
      real(float), intent(in) :: x(:),ab(0:,:)
      integer, intent(in) :: n
      Logical, optional :: Monic
      real(float) :: p(size(x),0:n)
      integer :: k
      logical :: monc
      ! ------------------------------------------
      ! Uses recurrent relations to calculate All Jacobi polynomials from 0 to n
      !  x  - x where polynomials are calculated for 0 to n
      !  n  - highest order polynomial calculated
      !  ab - recurrence coef., etc.
      !  Monic - .true. if monic polynomial
      ! Monic form gives
      !  Pn+1 = (x - ab(n,1))*Pn - ab(n,2)*Pn-1
      !  ab(n,3) - leading coefficient of conventional form
      ! Conventional form (Monic = false) give coefficients:
      !  Pn+1 = (x*ab(n,3) - ab(n,1))*Pn - ab(n,2)*Pn-1
      ! In either case:
      !  ab(n,4) - Integral(Pn^2*w(x))
      ! ------------------------------------------
      !write(*,'(a,2i4)')'Jacobi_All:',n,size(ab,1)
      monc = Default_monic;    if(present(Monic))monc = Monic
      p(:,0) = 1.0_float
      if(.not.monc)then
         p(:,1) = x*ab(0,3) - ab(0,1)
         do k = 1,n-1
            p(:,k+1) = (x*ab(k,3) - ab(k,1))*p(:,k) - ab(k,2)*p(:,k-1)
         enddo
      else
         p(:,1) = x - ab(0,1)
         do k = 1,n-1
            p(:,k+1) = (x - ab(k,1))*p(:,k) - ab(k,2)*p(:,k-1)
         enddo
      endif
   End Function Jacobi_All
   !----------------------------------------------------------------------------
   Function Jacobi_Series(an,ar,x,Monic)  result(u)
      real(float), intent(in) :: an(0:),ar(0:,:),x(:)
      Logical, optional :: Monic
      Real(float) :: u(size(x)),p(size(x),0:size(an)-1)
      integer :: n
      ! ------------------------------------------
      ! Function calculates the values of function which is a
      ! finite series of Jacobi polynomials
      !  an - coefficients of the polynomials
      !  x  - values where calculated
      !  alpha,beta - Jacobi weights
      !  ap - recurrence coefficients returned if argument present
      !  Monic - .true. if monic polynomial
      !  Shift - .true. if shifted to 0,1
      ! ------------------------------------------
      n = size(an) - 1
      p = Jacobi_All(x,n,ar,Monic)
      u = MatMul(p,an)
   End Function Jacobi_Series
   !----------------------------------------------------------------------------
   Function Legendre(x,n,order)  result(pn)
      Real(float) :: x(:)
      integer :: n, order,j,k
      Real(float) :: pn(size(x),0:n,0:order),ar(0:n-1,4),dx
      integer :: d2(0:n)
      ! ------------------------------------------
      ! This routine calculates values and derivatives of Legendre polynomial
      ! x - points where calculated
      ! n - highest degree
      ! order - highest derivative 0 to 2
      ! pn(i,j,k) - value at x(i) of jth polynomials, kth derivative
      ! ------------------------------------------
      pn = 0.0_float
      ar = x_jacobi(n,Monic=.false.)        ! recursion coefficients
      pn(:,:,0) = Jacobi_All(x,n,ar,Monic=.false.)
      if(order > 0)then
         dx = 1.0_float    ! 1st derivative
         pn(:,1,1) = pn(:,0,0)
         do k = 2,n
            dx = dx + 2.0_float     ! coefficients calculated on the fly
            pn(:,k,1) = dx*pn(:,k-1,0) + pn(:,k-2,1)
         enddo
      endif
      if(order > 1)then
         do j = 2,n
            d2(0:j) = Legendre_Deriv(j,2) ! calculates coefficients
            do k = 0,j-2
               if(d2(k) == 0)cycle
               pn(:,j,2) = pn(:,j,2) + Real(d2(k),float)*pn(:,k,0)
            enddo
         enddo
      endif
   End Function Legendre
   !----------------------------------------------------------------------------
   Function Legendre_Deriv(n,order)  result(b)
      integer, intent(in) :: n,order
      integer :: k,k1,b(0:n)
      ! ------------------------------------------
      ! This routine calculates coefficients of derivatives of Legendre polynomial
      ! in terms of the undifferentiated polynomials
      ! b(0:n) - coefficient values
      ! dp/dx = sum(b(i)*p(i)) for order 1
      ! order = 1 1st derivative
      ! order = 2 2nd derivative
      ! order = 3 difference of 2nd derivative polynomials
      ! ------------------------------------------
      b = 0
      Select Case(order)
      case(1)              ! coefficients of P'n
         k1 = 1 - mod(n,2) ! start at 1 if even, 0 if odd
         do k = k1,n-1,2
            b(k) = 2*k+1
         enddo
      case(2)              ! coefficients of P"n
         k1 = mod(n,2)     ! start at 0 if even, 1 if odd
         do k = k1,n-2,2
            b(k) = (2*k+1)*(n*(n+1) - k*(k+1))/2
         enddo
      case(3)
         k1 = mod(n,2)     ! start at 0 if even, 1 if odd
         do k = k1,n,2     ! coefficients of P"n+2 - P"n
            b(k) = (2*n+3)*(2*k+1)
         enddo
      Case Default
      End Select
   End Function Legendre_Deriv
!----------------------------------------------------------
   Function Jacobi_Deriv1(n,alpha,beta,Monic,Shift) result(c)
      integer, intent(in) :: n
      real(float), optional :: alpha,beta
      logical, optional :: Monic,Shift
      real(float) :: c(3),abx(2),a,b,xn,t
      logical :: monc,shft
      ! -----------------------------------------
      ! Produces the following coefficients for first derivative calcs.
      ! (1-x^2)Pn' = [c(1)*x + c(2)]Pn + c(3)*Pn-1
      ! Valid for convention on [-1,1] or Shifted [0,1]
      ! For monic must modify with leading terms from x_jacobi:
      !  c(3) = c1(3)*ar(n-1,3)/ar(n,3)
      !  or from two calls to Jacobi_Lead
      ! -----------------------------------------
      monc = Default_monic;   if(present(Monic))monc = Monic
      shft = Default_shift;   if(present(Shift))shft = Shift
      abx = abparm(alpha,beta); a = abx(1);  b = abx(2)
      xn = Real(n,float)
      t = (1.0_float)/(xn*2.0_float+a+b)
      c(1) = -xn
      c(2) = xn*(a-b)*t
      c(3) = (2.0_float)*(xn+a)*(xn+b)*t
      if(shft)then
         c(2) = (c(2) - c(1))*0.5_float
         c(3) = c(3)*0.5_float
      endif
      if(monc)then  ! convert to monic
        c(3) = c(3)*Jacobi_Lead(n-1,a,b,shft)/Jacobi_Lead(n,a,b,shft)
      endif
   End Function Jacobi_Deriv1
!----------------------------------------------------------
   Function Jacobi_Deriv2(n,alpha,beta,nd,Shift) result(c)
      integer, intent(in) :: n
      real(float), optional :: alpha,beta
      integer, optional :: nd
      logical, optional :: Shift
      real(float) :: c(3),abx(2),a,b,xn
      integer :: j
      logical :: shft
      ! -----------------------------------------
      ! Produces the following coefficients for 2nd derivative
      ! (1-x^2)Pn" = [c(1)*x + c(2)]Pn' + c(3)*Pn
      ! Coefficients valid for monic or conventional Pn
      ! -----------------------------------------
      shft = Default_shift;   if(present(Shift))shft = Shift
      abx = abparm(alpha,beta); a = abx(1);  b = abx(2)
      xn = Real(n,float)
      c(1) = 2.0_float+a+b
      c(2) = a-b
      c(3) = -xn*(xn+a+b+1.0_float)
      if(present(nd))then
         do j=3,nd         ! coefficients for 3rd and higher
            c(3) = c(3) + c(1)
            c(1) = c(1) + 2.0_float
         enddo
      endif
      if(shft)then
         c(2) = (c(2) - c(1))*0.5_float
      endif
   End Function Jacobi_Deriv2
!----------------------------------------------------------
   Function Jacobi_Deriv(n,alpha,beta,ab,Monic,Shift) result(c)
      integer, intent(in) :: n
      real(float), optional :: alpha,beta,ab(0:,:)
      logical, optional :: Monic,Shift
      logical :: monc,shft
      integer :: k
      real(float) :: c(0:n-1,0:n-1),a,b,abx(2),apb,ab2
      real(float) :: ad(0:n,3),zn(0:n+1)
      real(float) :: ar(0:n,4)
      logical :: ultra
      real(float), parameter :: eps = 1.0e-6_float
      ! -----------------------------------------
      ! This function calculates coefficients to express Jacobi polynomial
      ! 1st derivative as a sum of lower order Jacobi polynomial values
      ! dPn/dx = sum(d(n-1,0:n-1)*Pn(0:n-1))
      ! This code should work for all combinations of Monic & Shift
      ! -----------------------------------------
      if(n < 2)return   ! not set up for n = 0 or 1
      abx = abparm(alpha,beta); a = abx(1);  b = abx(2)  ! alpha & beta
      monc = Default_monic;   if(present(Monic))monc = Monic
      shft = Default_shift;   if(present(Shift))shft = Shift
      ultra = abs(a-b) < eps
      ! basic coefficients to Pk+1' = ad_k3*Pk + ad_k1*Pk' + ad_k2*Pk-1'
      apb = a + b;  ab2 = 0.0_float
      if(abs(apb) > eps)ab2 = (2.0_float)/apb  ! 2/(a+b)
      zn = (/(Real(k,float),k=0,n+1)/)    ! 0,1,2,...,n+1
      ad = 0.0_float
      if(.not.ultra)then
         ad(:,1) = zn(1:n+1)*ab2             ! 2(k+1)/(a+b) for bk or Pk' coeff
      endif
      ad(0,2) = ab2*0.5_float             ! 1/(a+b)
      ad(2:,2)= zn(3:n+1)/(zn(2:n)+apb)   ! (k+1)/(k+a+b) ck or Pk-1' coeff
      ad(:,3) = zn(1:n+1)                 ! k+1 for ak or Pk coeff.
      if(present(ab))then
         ar = ab  ! values passed in must not be shifted!
      else
         ar = x_jacobi(n+1,alpha,beta,Monic,Shift=.false.)   ! recursion coefficients
      endif
      if(shft .and. monc)then
         ar(:,1) = ar(:,1)- 0.5_float
      elseif(shft)then
         ar(:,1) = ar(:,1) - ar(:,3)*0.5_float
      endif
      ad(:,1:2) = ad(:,1:2)*ar(:,1:2)
      if(.not.monc)ad(:,3) = ad(:,3)*ar(:,3)
      if(Debug > 3)then
         write(iout,'(a,i2,a,2f8.2)')'Jacobi_Deriv: n =',n,', a,b =',a,b
!         Call ArrayPrint('Jacobi_Deriv - ad',transpose(ad))
      endif
      ! build coefficients recursively
      c = 0.0_float;    c(0,0) = ad(0,3)
      c(1,1) = ad(1,3); c(1,0) = ad(1,1)*c(0,0)
      do k = 2,n-1
         c(k,k) = ad(k,3)
         c(k,0:k-1) = ad(k,1)*c(k-1,0:k-1)
         c(k,0:k-2) = c(k,0:k-2) + ad(k,2)*c(k-2,0:k-2)
      enddo
   End Function Jacobi_Deriv
   !----------------------------------------------------------------------------
   Function Jacobi_Ends(n,alpha,beta,Monic,Shift)  result(p)
      Integer, intent(in) :: n
      Real(float), optional :: alpha,beta
      Logical, optional :: Monic,Shift
      Real(float) :: p(2),abx(2),an
      logical monc,shft
      !----------------------------------------------------
      ! Calculates polynomial values at left and right ends
      ! Could easily be modified to calculate end derivatives
      ! Direct calculation used, because Gamma(>170) fails
      !----------------------------------------------------
      abx = abparm(alpha,beta)
      monc = Default_monic;  if(present(Monic))monc = Monic
      shft = Default_shift;  if(present(Shift))shft = Shift
      p = GammaR(n,abx(2))
      if(abx(2) /= abx(1)) &
         p(2) = GammaR(n,abx(1))
      if(monc)then
         an = Jacobi_Lead(n,abx(1),abx(2),shft)
         p = p/an
      endif
      p(1) = p(1)*Real(1-2*mod(n,2),float) ! sign
!      write(iout,'(a,i5,2f7.2,10es24.15)')'Jacobi_Ends:',n,abx,p
      ! ----------------------------------
      Contains
      Function GammaR(n,c) result(p)
         integer :: n
         real(float) :: c,p
         ! calculate endpoint value: c - integer value not required
         !  p = Gamma(n+c+1)/[(n!)Gamma(c+1)
         if(c == 0.0_float)then
            p = 1.0_float
         elseif(c == 1.0_float)then
            p = Real(n+1,float)
         else
            p = GammaRatio(c+1.0_float,1.0_float,n)/Gamma(c+1.0_float)
         endif
      End Function
   End Function Jacobi_Ends
   !----------------------------------------------------------------------------
   Function Jacobi_Lead(n,alpha,beta,Shift)  result(an)
      Integer, intent(in) :: n
      Real(float), optional :: alpha,beta
      Logical, optional :: Shift
      Real(float) :: an,abx(2),xn,c,ab   !,zn(0:n-1)
      logical shft
      !----------------------------------------------------
      ! Calculate leading coefficient:
      !   n - number or degree of polynomial
      !   alpha, beta -  optional Jacobi weight parameters
      !   Shift -> shft optional true for shifted
      ! The formula calculated (non shifted) is:
      !   an = Cn*Gamma(n+c+1)*Gamma(n+c+0.5)/[Gamma(n+1)*Gamma(n+2c+1)]
      ! where:
      !   c = (a+b)/2
      !   Cn = [2^(n+a+b)]/sqrt(pi)
      ! Gamma function fails for Gamma(n), n >~ 170 (Real*8)
      ! The normal expression has arguments O(2n), which can
      ! be reduced to two formulas of O(n). In this form
      ! we end up with four gamma functions, 2 divided by 2
      ! all with different arguments in general. Can be reduced in
      ! special cases. We use the GammaRatio function which is called twice.
      ! A direct calculation for small n and Stirling approximation for
      ! large n.
      !----------------------------------------------------
      abx = abparm(alpha,beta)
      shft = Default_shift;  if(present(Shift))shft = Shift
      xn = Real(n,float)
      ab = abx(1)+abx(2);  c = ab*0.5_float
      an = GammaRatio(c+0.5_float,ab+1.0_float,n)
      an = GammaRatio(c+1.0_float,1.0_float,n)*an/sqrt(pi)
      if(shft)then
         an = an*((2.0_float)**(xn*2.0_float+ab))
      else
         an = an*((2.0_float)**(xn+ab))
      endif
      if(Debug > 0)write(iout,'(a,i5,2f5.1,g14.7,2(a,2f7.2))')'Jacobi_Lead: n,a,b,an =',n,abx,an, &
         ', Gamma:',c+0.5_float,ab+1.0_float,', Gamma:',c+1.0_float
   End Function Jacobi_Lead
   !----------------------------------------------------------------------------
   Function GammaRatio(a,b,n) result(g)
      integer :: n,k
      real(float) :: a,b,zn(n),g,z1,z2,zr
      Real(floatq) :: z1q,z2q ! poor accuracy real exponentiation
      real(float), parameter :: eps0 = min(epsilon(a),1.0e-9_float)*0.5_float
      real(float), parameter :: zswitch =max(10.0_float,(7.0e-5_float/eps0)**(1.0_float/6.0_float))
      real(float), parameter :: eps = 1.0e-6_float
      ! --------------------------------------
      ! Calculate gamma(n+a)/gamma(n+b)
      ! --------------------------------------
      z1 = real(n,float) + a;   z2 = real(n,float) + b
      if(abs(a-b) < eps)then
         g = 1.0_float
      elseif(min(z1,z2) < zswitch)then    ! based on next term in expansion (see Sfunc)
         zn = (/(Real(k,float),k=0,n-1)/) ! 0,1,2,...,n-1
         g = product((zn+a)/(zn+b))
         g = g*gamma(a)/gamma(b)
         !write(iout,'(a,i7,2f6.2,es24.17,L2)')'GammaRatio: ',n,a,b,g,min(z1,z2) < zswitch
      else   ! Stirling approximation good to 1e-16 for n > ~50 to 60
       ! Note: zr term must be in higher precision, not sure why
         z1q = real(n,floatq) + real(a,floatq);   z2q = real(n,floatq) + real(b,floatq)
         z1q = (z1q/z2q)**(z1q-0.5_floatq);       zr = real(z1q,float)
         g = exp(b-a)*(z2**(a-b))*zr*Sfunc(z1)/Sfunc(z2)
         !write(iout,'(a,i7,2f6.2,es24.17,L2)')'GammaRatio Stirling: ',n,a,b,g,min(z1,z2) < zswitch
      endif
      contains
      ! --------------------------
      Function Sfunc(x)  result(s)
         real(float) :: x,s,z
         real(float), parameter :: a(0:5) = (/1.0_float,1.0_float/12.0_float,1.0_float/288.0_float, &
         -139.0_float/51840.0_float,-571.0_float/2488320.0_float,163879.0_float/209018880.0_float/)
         z = 1.0_float/x
         s = a(0) + z*(a(1) + z*(a(2) + z*(a(3) + z*(a(4) + z*a(5)))))
      End Function Sfunc
   End Function GammaRatio
!----------------------------------------------------------
   function Jacobi_Transform(n,t,g,s,shft,func)  result(a)
      integer, intent(in) :: n,t,g,s
      logical, intent(in) :: shft
      real(float) :: a(0:n)
      real(float) :: ab(2),wbfac
      real(float) :: ar(0:n,4),xw(s:n+2,4),fx(n+1),pn(n+1,0:n),x(n+1),xs(n+1)
      integer :: k,n1
      Interface
         Function func(x)  result(f)
            include 'defs.fi'
            Real(float), intent(in) :: x(:)
            Real(float) :: f(size(x))
         End Function func
      End Interface
      ! --------------------------------
      ! calculate Jacobi transform coefficients for function "func"
      ! n - highest polynomial degree
      ! t - determine Jacobi alpha,beta as follows:
      !   0,0 - Gauss; 1,1 - Lobatto; 1,0 - RadauR; 0,1 - RadauL
      ! g - 0,1,2,3 for nonsymmetric, planar, cylindrical & spherical
      ! s - 0,1 for nonsymmetic & symmetric
      ! func- f(x) for which trasform is calculated
      ! Integrals are approximated by quadrature, so this is actually interpolation
      ! using Gauss-Jacobi requires n+1 points to get accuracy through 2n+1
      ! Is is possible to adjust last factor to do with n quadrature points??
      ! --------------------------------
      n1 = n + 1
      ab = Polyset(t,g)  ! alpha & beta
      if(Debug > 0)write(iout,'(a,i6,3i3,2f5.1,L2)')'Jacobi_Transform:',n,t,g,s,ab,shft
      xw = Jacobi_Quadrature(n1,ab,g,wbfac)
      ar = x_jacobi(n1,ab(1),ab(2),Shift=shft)! recursion coefficients
      x = xw(1:n1,1)
      xs = x
      if(g > Nonsymmetric)then
         x = x**2
         select case(t)
         case(Gauss)! modifications, extra factor of two
            xw(1:n1,4) = xw(1:n1,4)*2.0_float
         case(Lobatto)
            xw(1:n1,4) = xw(1:n1,4)*(1.0_float - x)*2.0_float
         end select
      elseif(shft)then
         xs = (1.0_float + x)*0.5_float
         x = (1.0_float + x)*0.5_float
      endif
      fx = func(xs)
      pn = Jacobi_All(x,n,ar,Monic=.false.)     ! polynomials
      select case(t+10*g)
      case(Lobatto)  ! convert Lobatto & Radau quadrature to Jacobi-Gauss
         xw(1:n1,4) = xw(1:n1,4)*(1.0_float - xw(1:n1,1)**2)
      case(RadauL)
         xw(1:n1,4) = xw(1:n1,4)*(1.0_float + xw(1:n1,1))
      case(RadauR)
         xw(1:n1,4) = xw(1:n1,4)*(1.0_float - xw(1:n1,1))
      case default
         continue
      end select
      a = Matmul(fx*xw(1:n1,4),pn)/ar(:,4)
      ! not sure whether to convert these are not (make it an option?)
      select case(g + (t-Lobatto)*10) ! convert Lobatto to full polynomials
      case(Planar)      ! coefficient of Pn(1,1) n = 0,2,4,...
         a = a*real((/(k,k=1,n+1)/),float)/real((/(2*k+1,k=0,n)/),float)
      case(Spherical)   ! coefficients of Pn(1,1)/x n =1,3,5,...
         a = a*0.5_float
      end select
      if(Debug > 2)then
         xw(1:n1,2) = Matmul(pn,a) ! to check (wrong for Lobatto symmetric)
         Call arrayprint('an, Recursion:',a,ar,fmtx=fmt22)
         call Arrayprint('Jacobi_Trans: x,f(x),x,f(x),wb,w,Pn',x,fx,xw(1:n1,:),pn,fmtx=fmt22)
      endif
   end function Jacobi_Transform
!----------------------------------------------------------
   Function Legendre_Transform(x,w,ityp,shift) result(Q)
      real(float), intent(in) :: x(0:),w(0:)
      integer, intent(in) :: ityp
      logical, optional :: shift
      logical :: shft
      real(float) :: Q(0:size(x)-1,0:size(x)-1)
      real(float) :: xc(0:size(x)-1),wc(0:size(x)-1)
      ! -------------------------------------
      ! This routine calculates Legendre Transforms
      ! x,w - nodal points and quadrature weights
      ! nc - interior nodes
      ! ityp - Gauss or Lobatto type
      ! pn - Legendre polynomials at points
      ! Q - Legendre transform matrix
      ! ------------------------------------------------
      shft = Default_shift;  if(present(Shift))shft = Shift
      if(shft)then
         wc = w*2.0_float;  xc = x*2.0_float - 1.0_float
      else
         wc = w;  xc = x
      endif
      select case(ityp)
      case(Gauss)
         Q = Legendre_Jacobi_Trans(xc,wc) ! Legendre Transform via Jacobi Transform
      case(Lobatto)
         Q = Legendre_Trans(x,w,shft)        ! straight Legendre Transform
      end select
   End Function Legendre_Transform
!----------------------------------------------------------
   Function Legendre_Trans(x,w,shft) result(Q)
      real(float), intent(in) :: x(0:),w(0:)
      logical, intent(in) :: shft
      real(float) :: Q(0:size(x)-1,0:size(x)-1)
      real(float) :: ar(0:size(x)-1,4),pn(0:size(x)-1,0:size(x)-1)
      integer :: k,n,nc
      ! -------------------------------------
      ! This routine calculates Legendre Transforms
      ! x,w - Lobatto quadrature points and weights
      ! nc - interior nodes
      ! n - polynomial degree
      ! pn - Legendre polynomials at points
      ! Q - Legendre transform matrix
      ! ar(k,4) = 2/(2k+1) for k < n
      ! ar(n,4) = 2/n
      ! ------------------------------------------------
      nc = size(x)-2;  n = nc + 1
      write(iout,'(a,2i4,L2,2i4)')'Legendre_Tran - Lobatto pts:',nc,n,shft,shape(ar)
      ar = x_jacobi(n+1,Monic=.false.,Shift=shft)        ! recursion coefficients
      pn = Jacobi_All(x,n,ar,Monic=.false.)
      ar(n,4) = ar(n,4)*Real(2*n+1,float)/Real(n,float)  ! modify for last point
      do k = 0,nc+1
         Q(k,:) = w*pn(:,k)/ar(k,4)
      enddo
      ! ------------------------------------------------
      if(Debug > 3)then
         Call ArrayPrint('Legendre_Trans x,w,ar4,pn:',x,w,ar(:,4),pn)
      !   Call ArrayPrint('Legendre_Trans Matrix Q:',Q)
      endif
   End Function Legendre_Trans
!----------------------------------------------------------
   Function Legendre_Jacobi_Trans(x,w) result(Q)
      real(float), intent(in) :: x(0:),w(0:)
      real(float) :: Q(0:size(x)-1,0:size(x)-1)
      real(float) :: ar(0:size(x)-3,4),xn(0:size(x)-3),pj(1:size(x)-2,0:size(x)-3)
      real(float) :: Qj(0:size(x)-3,1:size(x)-2),Qs(0:size(x)-1,0:size(x)-1)
      real(float), parameter :: alpha = 1.0_float
      logical, parameter :: shft = .false.
      integer :: n,k,nc,nc1
      ! -------------------------------------
      ! This code calculates a Legendre Transform matrix for Gauss points
      ! Converts a Jacobi(1,1) transform matrix
      !  x,w - Gauss points and quadrature weights on [0,1] basis
      !  an(:,0) - the integrals of Pn(1,1)^2
      !  Qj - Jacobi transform
      !  Q  - Legendre transform
      !  Qs - Q scratch array (intermediate transform)
      !  n  - highest order of the Jacobi polynomial
      !  nc - interior pts.
      !  ityp - Gauss or Lobatto
      !  ar(k,4) = 8(k+1)/[(2k+3)(k+2)] k < n
      !  ar(k,4) = 8(k+1)/(k+2)^2       k = n
      ! -------------------------------------
      nc = size(x)-2; n = nc-1;  nc1 = nc+1
      write(iout,'(a,10i4)')'Legendre_Jacobi_Trans - Gauss pts:',nc,n,shape(ar)
      xn = (/(Real(2*k,float),k=1,n+1)/)
      xn = -xn/(xn + 1.0_float) ! conversion from Jacobi to Legendre
      ar = x_jacobi(n+1,alpha,Monic=.false.,Shift=shft)  ! recursion coefficients
      pj = Jacobi_All(x(1:nc),n,ar,Monic=.false.)        ! Jacobi (1,1)
      ar(n,4) = ar(n,4)*Real(2*n+3,float)/Real(n+2,float)
      do k = 0,n
         Qj(k,:) = w(1:nc)*pj(:,k)/ar(k,4)
      enddo
      ! calculate Legendre Transform, from Jacobi transform
      do k = 0,n
         Qj(k,:) = xn(k)*Qj(k,:)  ! *(-1)(2k+2)/(2k+3)
         Qs(k,0) = (-0.5_float)*sum(Qj(k,:)*(1.0_float - x(1:nc)))
         Qs(k,1:nc) = Qj(k,:)
         Qs(k,nc1) = (-0.5_float)*sum(Qj(k,:)*(1.0_float + x(1:nc)))
      enddo
      Qs(n+1:n+2,:) = 0.0_float
      Q(0,0) =  0.5_float - Qs(0,0);  Q(0,nc1) = 0.5_float - Qs(0,nc1)
      Q(1,0) = -0.5_float - Qs(1,0);  Q(1,nc1) = 0.5_float - Qs(1,nc1)
      Q(0:1,1:nc) = -Qs(0:1,1:nc)
      do k = 2,n+2
         Q(k,:) = Qs(k-2,:) - Qs(k,:)
      enddo
      ! ----------------------------------------
      if(Debug > 3)then
         Call ArrayPrint('Recurrence:',ar)
         Call ArrayPrint('Legendre_Jacobi_Trans x,w,ar4,pn:',x(1:nc),w(1:nc),ar(:,4),pj)
         Call ArrayPrint('Jacobi Transform Matrix:',Qj)
         Call ArrayPrint('Intermediate Transform Matrix:',Qs(0:n,:))
      !   Call ArrayPrint('Legendre Transform Matrix:',Q)
      endif
   End Function Legendre_Jacobi_Trans
   !----------------------------------------------------------------------------
   Function x_Jacobi(n,alpha,beta,Monic,Shift)  result(ab)
      Integer, intent(in) :: n
      Real(float), optional :: alpha,beta
      Logical, optional :: Monic,Shift
      Real(float) :: ab(0:n-1,4)
      Real(float) :: abx(2),a,b,am(0:n-1,2)
      Real(float) :: b2a2,apb,a1,b1,ab1,gab1
      Real(float) :: zn(0:n),z2nab(0:n),r(0:n),fac2
      integer :: k,n1
      logical monc,shft
      !-------------------------------------------------------------------------
      ! Function to calculate the coefficients of a Jacobi polynomial
      ! Monic = true/false for monic or conventional form      (Default = .true.)
      ! Shift = true/false shifted (0,1) vs not shifted (-1,1) (Default = .false.)
      ! Monic form gives
      !  Pn+1 = (x - ab(n,1))*Pn - ab(n,2)*Pn-1
      !  ab(n,3) - leading coefficient of conventional form
      ! Conventional form (Monic = false) give coefficients:
      !  Pn+1 = (x*ab(n,3) - ab(n,1))*Pn - ab(n,2)*Pn-1
      ! In either case:
      !  ab(n,4) - Integral(Pn^2*w(x))
      !-------------------------------------------------------------------------
      abx = abparm(alpha,beta); a = abx(1);  b = abx(2)
      monc = Default_monic;   if(present(Monic))monc = Monic
      shft = Default_shift;   if(present(Shift))shft = Shift
      fac2 = 2.0_float
      if(shft)fac2 = 1.0_float
      ab = 0.0_float
      n1 = n-1
      a1 = a + 1.0_float;  b1 = b + 1.0_float;  ab1 = a1 + b;
      b2a2 = b*b - a*a;  apb = a + b;  gab1 = gamma(a1)*gamma(b1)
      zn = (/(Real(k,float),k=0,n)/);   z2nab = (2.0_float)*zn + apb
      r(0) = (a1+b1)/fac2
      r(1:) = (2.0_float*zn(1:)+a1+b1)*(2.0_float*zn(1:)+ab1)/((zn(1:)+1.0_float)*(zn(1:)+ab1)*fac2)
      ab(0,1) = (b-a)/(apb+2.0_float);
      ab(0,2) = ((2.0_float)**(apb+1.0_float))*gab1/gamma(a1+b1);
      ab(0,3) = 1.0_float
      ab(0,4) = gab1*(fac2**ab1)/gamma(ab1+1.0_float)
      if(n > 1)then
         ab(1:,1) = b2a2/(z2nab(1:n1)*(z2nab(1:n1) + 2.0_float))
         ab(1,2) = (4.0_float)*(a1)*(b1)/((apb+3.0_float)*(a1+b1)**2);
         ab(2:,2) = (4.0_float)*zn(2:n1)*(zn(2:n1)+a)*(zn(2:n1)+b)*(zn(2:n1)+apb) &
               /((z2nab(2:n1)+1.0_float)*(z2nab(2:n1)-1.0_float)*(z2nab(2:n1)**2));
      endif
      Call Shift_Poly
      do k=1,n-1
         ab(k,3) = ab(k-1,3)*r(k-1)
         ab(k,4) = ab(k-1,4)*ab(k,2)*r(k-1)**2
      enddo
      Call Convert_Monic
      ! ------------------------------------
      contains
      ! ------------------------------------
      Subroutine Convert_Monic
         if(.not.monc)then
            am(:,1) = ab(:,1)*r(0:n-1)
            am(1:,2) = ab(1:,2)*r(1:n-1)*r(0:n-2)
            am(0,2) = a*b*(a1+b1)  ! default?
            if((ab1*apb)>0.0_float)am(0,2) = am(0,2)/(ab1*apb)
            ab(:,3) = r(0:n-1)
            ab(:,1:2) = am
         endif
      End Subroutine Convert_Monic
      Subroutine Shift_Poly
         if(shft)then
            ab(:,1) = (1.0_float + ab(:,1))*0.5_float
            ab(:,2) = ab(:,2)*0.25_float
         endif
      End Subroutine Shift_Poly
   End Function x_Jacobi
   !----------------------------------------------------------------------------
   Function d_Jacobi(n,alpha,beta,Shift)  result(dn)
      Integer, intent(in) :: n
      Real(float), optional :: alpha,beta
      Real(float) :: dn(0:n)
      Logical, optional :: Shift
      Real(float) :: abx(2),a,b,ab1,xk,xk2ab
      Logical :: shft
      Integer :: k,k1
      ! --------------------------------
      ! this function calculates the integrals of Pn^2
      ! on interval (-1,1) or (0,1) if Shift = .true.
      ! Gamma(ab1+1.0_float) = ab1*Gamma(ab1) is an identity
      ! used in dn(0) to avoid divide by zero for a = b = -0.5
      ! --------------------------------
      shft = Default_shift;   if(present(Shift))shft = Shift
      abx = abparm(alpha,beta)  ! handles defaults
      a = abx(1);  b = abx(2); ab1 = a + b + 1.0_float
      dn(0) = Gamma(a+1.0_float)*Gamma(b+1.0_float)/Gamma(ab1+1.0_float)
      if(abs(ab1) < epsilon(a)*10.0_float)then
         dn(1) = dn(0)  ! from L'Hopitals rule
         k1 = 2
      else
         k1 = 1
         if(.not.shft)dn(0) = dn(0)*(2.0_float)**(a+b+1.0_float)
      endif
      do k = k1,n
        xk = Real(k,float);  xk2ab = (2.0_float)*xk+ab1
        dn(k) = dn(k-1)*(xk+a)*(xk+b)*(xk2ab-2.0_float)/(xk*(xk+a+b)*(xk2ab))
      enddo
   End Function d_Jacobi
   !----------------------------------------------------------------------------
   Function Polyset(t,g) result(ab)
      integer, intent(in) :: t,g
      real(float) :: ab(2)             !Gauss,     Lobatto,  Chebyshev, RadauR,   RadauL,  Chebyshev1
      real(float), parameter :: a(6) = (/0.0_float,1.0_float,0.5_float,1.0_float,0.0_float,-0.5_float/)
      real(float), parameter :: b(6) = (/0.0_float,1.0_float,0.5_float,0.0_float,1.0_float,-0.5_float/)
      ! -----------------------------------------------------
      ! alpha, beta & n for calculation
      ! t - type - 1-6 => Gauss,Lobatto,Chebyshev,RadauR,RadauL,Chebyshev1
      ! g - geometry, 0/1/2/3 => nonsymmetric/full/planar/cyl/spherical
      ! -----------------------------------------------------
      ab(1) = a(t);  ab(2) = b(t)
      if(t <= Lobatto .and. g >= Planar)then
         ab(2) = Real(g-2,float)*(0.5_float)
      endif
   End Function Polyset
   !----------------------------------------------------------------------------
   Function abparm(alpha,beta)  result(ab)
      Real(float), optional :: alpha,beta
      Real(float) :: ab(2)
      ! this function handles defaulting of a and b
      ! if neither specified, use Legendre alpha = beta = 0
      ! if only one is specified, set other one to it, alpha = beta
      if(present(alpha) .and. present(beta))then
         ab(1) = alpha
         ab(2) = beta
      elseif(present(alpha))then
         ab(1:2) = alpha
      elseif(present(beta))then
         ab(1:2) = beta
      else
         ab = 0.0_float
      endif
   End Function abparm
   ! -----------------------------------------
   Function FirstLT(x,x0) result(k)
      real(float) :: x(:),x0
      integer :: k
      ! function to find the first value of k where x(k) < x0
      do k = 1,size(x)
         if(x(k) < x0)exit
      enddo
   End Function FirstLT
   ! -----------------------------------------
   Function IFirstLT(ix,i0) result(k)
      integer :: ix(:),i0
      integer :: k
      ! function to find the first value of k where x(k) < x0
      do k = 1,size(ix)
         if(ix(k) < i0)exit
      enddo
   End Function IFirstLT
   ! -----------------------------------------
   Function FirstGT(x,x0,Back) result(k)
      real(float) :: x(:),x0
      integer :: k,nx
      logical, optional :: Back  ! option to search backward
      logical :: bk
      ! function to find the first value of k where x(k) > x0
      bk = .false.;  if(present(Back))bk = Back
      nx = size(x)
      if(bk)then
         do k = nx,1,-1
            if(x(k) > x0)exit
         enddo
         k = nx - k + 1
      else
         do k = 1,size(x)
            if(x(k) > x0)exit
         enddo
      endif
   End Function FirstGT
   !----------------------------------------------------------------------------
End Module Orth_Poly

