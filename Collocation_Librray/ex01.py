import numpy as np
import numpy.linalg as linalg
import math
import arrayprint as ap
import occ
from occ import OrthColloc
def Example01(n,Thiele,Var0,Bi):
   ''' Solve Example BVP, first order, variable coefficients '''
   # implements catalyst pellet problem, linear source, variable coefficient
   # problem is solved twice: (1) collocation & (2) Galerkin (Lobatto) or Moments (Gauss)
   global ptyp,geom,sym,shift,Var1,Var2
   n1 = n + 1
   nt = n + 2
   nextra = 2  # extra quadrature points in mass matrix calc.
   shift = True   # [-1,1] or [0,1], i.e. shifted
   sym = min(1,geom)
   Th2 = (2.0*Thiele)**2
   Var1 = Var0*Th2
   Var2 = 4.0*(1.0 - Var0)*Th2
   Bi2 = 2.0*Bi
   ConstCoef = (Var2 < 1.0e-8)
   Dirichlet = Bi2.max() < 1.0e-8
   if Dirichlet:
      i0 = 1
      i1 = n1
   else:
      i0 = 0
      i1 = nt
   oc = OrthColloc(n,ptyp,geom,shift)
   x = oc.Xpoints()
   w = oc.WeightQ()
   Ac = oc.Deriv1()
   Bc = oc.Deriv2()
   Cc = oc.Stiffness()
   Dm = oc.MassMatrix()
   rx = VarCoef(x)
   Df = oc.MassMatrix(VarCoef,nextra)
   ReportBasic(Thiele,Bi,x,w,rx,Ac,Bc,Cc,Dm,Df,nextra)
   y = np.zeros((2,nt))
   flux = np.zeros((2,3))
   A = np.empty((nt,nt))
   Aflux = np.empty((2,nt))
   # 0. Solve with collocation - symmetric matrix formulation
   rw = rx*w
   A[:,:] = Cc    # note: A = Cc makes A equivalent to Cc
   for i in range(nt):
      A[i,i] += rw[i]
   Aflux[0,:] = A[i,i]
   Aflux[1,:] = A[n1,i]
   A[0,0] += Bi2[0]
   A[n1,n1] += Bi2[1]
   y[0,i0:i1] = linalg.solve(A[i0:i1,i0:i1],rw[i0:i1])
   flux[0,0] = (rw[0]-A[0,:].dot(y[0,:]))/Th2   #/#
   flux[0,1] = (rw[n1]-A[n1,:].dot(y[0,:]))/Th2  #/#
   flux[1,0] = Ac[0,:].dot(y[0,:])/Th2    #/#
   flux[1,1] = -Ac[n1,:].dot(y[0,:])/Th2   #/#
   flux[:,2] = flux[:,0] + flux[:,1] # total
   ReportFinal(A,rw,y[0,:],flux,'Symmetric Collocation Formulation')
   # 1. Solve with Galerkin method or method of moments
   A = Cc + Df
   rw = Df.sum(axis=1)  # = w*rx except for small n
   Aflux[0,:] = A[0,:]
   Aflux[1,:] = A[n1,:]
   A[0,0] += Bi2[0]
   A[n1,n1] += Bi2[1]
   y[1,i0:i1] = linalg.solve(A[i0:i1,i0:i1],rw[i0:i1])
   flux[0,0] = (rw[0]-Aflux[0,:].dot(y[1,:]))/Th2   #/#
   flux[0,1] = (rw[n1]-Aflux[1,:].dot(y[1,:]))/Th2  #/#
   flux[1,0] = Ac[0,:].dot(y[1,:])/Th2    #/#
   flux[1,1] = -Ac[n1,:].dot(y[1,:])/Th2  #/#
   flux[:,2] = flux[:,0] + flux[:,1] # total
   ReportFinal(A,rw,y[1,:],flux,'Galerkin/moments method')
   ReportTabular(oc,x,y,Bc)

def ReportTabular(oc,x,y,Bc):
   ''' Gemerate tabular report, calculated on finer x. '''
   nshapey = y.shape
   ny = nshapey[0]
   nx = nshapey[1]
   nv = 22
   xv = np.arange(nv-1,-1,-1)/(nv-1)   #/#
   xv = (np.cos(xv*occ.pi_const) + 1.0)*0.50
   rv = VarCoef(xv)
   Lc = oc.Lcoeff(xv)
   v = np.empty((nv,3))
   for i in range(ny):
      yv = Lc.dot(y[i,:])
      d2c = Bc.dot(y[i,:])
      v[:,0] = Lc.dot(d2c)
      v[:,1] = rv*(1.0 - yv)
      v[:,2] = v[:,0] + v[:,1]
      ap.arrayprint('x, y, d2y, r(x,y), Resid',xv,yv,v,fmtf='%14.7g')

def ReportFinal(A,b,y,flux,title):
   ''' Report final results. '''
   sfmt = '%14.7g'
   print('\nFinal Results -',title,file=ap.file())
   ap.arrayprint('Matrix,rhs,solution:',A,b,y,fmtf=sfmt)
   ap.arrayprint('Left, right & total flux (correct & by derivative):',flux,fmtf=sfmt)

def ReportBasic(Thiele,Bi,x,w,rx,A,B,C,D,Df,nextra):
   ''' Report basic data for problem. '''
   global ptyp,geom,sym,shift,Var1,Var2
   sfmt = '%14.7g'
   txtt = str(occ.txt_type[ptyp]).strip("[']")
   txtg = str(occ.txt_geom[geom]).strip("[']")
   n = x.size - 2
   print('\nSolution of BVP with',txtt,'points and n =',n,file=ap.file())
   print('Thiele parameter =',Thiele,' and Bi =',Bi,file=ap.file())
   ap.arrayprint("\nPoints, Quad. Wts., Rate coefficient:",x,w,rx,fmtf='%18.15f')
   ap.arrayprint("\nFirst Derivative: ",A,fmtf=sfmt)
   ap.arrayprint("\nSecond Derivative: ",B,fmtf=sfmt)
   ap.arrayprint("\nStiffness Matrix: ",C,fmtf=sfmt)
   ap.arrayprint("\nMass Matrix, f(x) = 1: ",D,fmtf=sfmt)
   ap.arrayprint("\nMass Matrix, with cubic f(x) and "+str(nextra)+" extra pts.:",Df,fmtf=sfmt)

def VarCoef(x):
   ''' Calculate variable coefficients at values of x. '''
   global ptyp,geom,sym,shift,Var1,Var2
   f = Var1 + Var2*(1.5 - x)*x*x  # BVP with cubic coef., 0.2 to 1.8
   return(f)

def main():
   ''' Main program to read data and launch solution routine. '''
   global ptyp,geom,sym,shift
   geom = occ.Geom.Nonsymmetric.value
   sym = min(1,geom)
   c = input('Enter Thiele parameter >')
   Thiele = float(c)
   c = input('Variable coefficients y/n? >')
   Var0 = 1.0
   if (c == 'y'):
      Var0 = 0.2
   Bi = np.empty(2)
   Bi[:] = 0.0    # change manually if desired
   c = input('Enter g/l/r/x for Gauss/Lobatto/Radau(R/L) >')
   ptyp = ('glcrx').find(c) + 1
   if(c == 'c'):
      ptyp = occ.Typ.Lobatto.value
      print('** Chebyshev not currently supported **')
   c = input('Enter n values separated by whitespace: first last delta >')
   nums = [int(i) for i in c.split(' ')]
   c = input('Enter root of file name >')
   ap.fopen(c+'.dat')
   for n in range(nums[0],nums[1]+1,nums[2]):
      Example01(n,Thiele,Var0,Bi)

main()

