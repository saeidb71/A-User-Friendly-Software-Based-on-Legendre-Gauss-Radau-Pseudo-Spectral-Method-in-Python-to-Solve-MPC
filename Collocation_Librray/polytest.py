def main():
   import numpy as np
   import arrayprint as ap
   import jacobi as jp
   ap.fopen('polytest.dat')
   pi_const = jp.pi_const
   n = 6          # manually set values
   nd = 2         # derivative
   n0 = 1         # Pn to P_n-n0, how many lower?
   abeta = [1,0]  # alpha & beta
   Monic = False
   Shift = False
   Ultraspherical = (abeta[0] == abeta[1])
   print('n, ab, n0, nd =',n,abeta,n0,nd,file=ap.file())
   # recursion coefficients
   ar = jp.jac_reccoef(n,abeta,Monic,Shift)
   ap.arrayprint('\nRecursiont coefficients:',ar)
   # derivative relations
   c = np.empty((3,3))
   c[0,:] = jp.jac_deriv1(n,abeta,Monic,Shift)  # 1st derivative
   c[1,:] = jp.jac_deriv2(n,abeta,Shift)        # 2nd
   c[2,:] = jp.jac_deriv2(n,abeta,Shift,nd=3)   # 3rd
   # end points
   abeta1 = [abeta[0]+1,abeta[1]+1]
   p0 = jp.jac_ends(n,abeta,Monic,Shift)
   pd0 = jp.jac_ends(n-1,abeta1,Monic,Shift)
   pd0 = pd0*(n + abeta[0] + abeta[1] + 1)/2 #/#
   print('\nend points:',p0,pd0,file=ap.file())
   ap.arrayprint('\nDerivatives relationships:',c)
   # calculate polynomial and derivatives
   x0 = -1.0
   x1 = 1.0
   ncalc = 21
   x = np.linspace(x0,x1,num=ncalc)
   p = jp.jac_poly_recurrent(x,n,abeta,nd,n0,Monic,Shift)
   if nd > 0: # correct endpoint for 1st derivative
      if x1 == 1.0:  p[n0+nd-1,ncalc-1] = pd0[1]
      if x0 == -1.0: p[n0+nd-1,0] = pd0[0]
   ap.arrayprint("\nx,Pn,Pn':",x,np.transpose(p))
   # calculate nodes and weights
   w,wbfac = jp.jac_quad(n,abeta)
   ap.arrayprint("\nx,theta,Wb,W:, wb factor = "+str(wbfac),w.transpose())
   ap.fclose()

main()
   

