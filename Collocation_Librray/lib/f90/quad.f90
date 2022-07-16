Subroutine Quadrature(n,ab,geom,w,wbfac,nneg)
   Use Orth_Poly, only : Jacobi_Quadrature
   include 'defs.fi'
   integer, intent(in) :: n,geom
   real(float), intent(in) :: ab(2)
   real(float), intent(out) ::  w(max(0,min(1,geom)):n+1,0:3),wbfac
   integer, intent(out) :: nneg
   ! quadrature points and weight returns n+1 or n+2 x 4 array:
   ! w(:,i) - abs(x), theta, wbary, wquad
   ! wbfac - scaling factor
   ! nneg - negative x's
   w = Jacobi_Quadrature(n,ab,geom,wbfac,Nneg)
end subroutine Quadrature
