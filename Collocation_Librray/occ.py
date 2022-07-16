import numpy as np
import arrayprint as ap
import jacobi as jp
from enum import Enum
Debug = 0
pi_const = 3.141592653589793238462643383279
txt_type = (['Generic'],['Gauss'],['Lobatto'],['Chebyshev'],['Radau Right'],['Radau Left'],['Chebyshev 1'])
txt_geom = (['No symmetry'],['Planar'],['Cylindrical'],['Spherical'])
txt_symm = (['Nonsymmetric'],['Symmetric'])
class Typ(Enum):
   ''' e.g. occ.Typ.Gauss.value = numerical value '''
   Generic = 0
   Gauss = 1
   Lobatto = 2
   Chebyshev = 3
   RadauR = 4
   RadauL = 5
   Cheb1 = 6

class Geom(Enum):
   ''' e.g. occ.Geom.Planar.value = numerical value '''
   Nonsymmetric = 0
   Planar=1
   Cylindrical=2
   Spherical=3

def _Wshift(n,w,wbfac):
   ''' Convert nonsymmetric weights and points on [-1,1] to [0,1].'''
   wbfac = (2.0**n)*wbfac
   w[0,:] = (w[0,:]+1.0)*0.5  # X
   w[2,:] = 2.0*w[2,:]        # Wb
   w[3,:] = 0.5*w[3,:]        # Wq
   return(w,wbfac)

def _Deriv1(sym,xx,wb):
   ''' Calculate first derivative matrix. '''
   nt = xx.size
   if Debug > 0:
      print('Deriv1:',nt,sym)
   if sym:
      x = xx**2
   else:
      x = xx
   A = np.empty((nt,nt))
   for i in range(nt):
      px = x[i] - x
      px[i] = 1.0
      A[i,:] = wb/(px*wb[i])  #/#
      A[i,i] = 1.0 - np.sum(A[i,:])
   if sym:
      for i in range(nt):
         A[i,:] = 2.0*xx[i]*A[i,:]
   return(A)
   
def _Deriv1odd(geom,sym,x,wb):
   ''' Calculate first derivative matrix for an odd function. '''
   nt = x.size
   An = np.empty((nt,nt))
   A = _Deriv1(sym,x,wb)
   if sym:
      gx = float(geom)
      for i in range(nt):
         An[i,:] = x[i]*A[i,:]/x[:]
         An[i,i] = An[i,i] + gx/x[i]
   else:
      An = A   # in case of call for nonsymmetric case
   return(An)

def _Deriv2(geom,sym,x,wb):
   ''' Calculate 2nd derivative or Laplacian matrix. '''
   A = _Deriv1(sym,x,wb)
   if Debug > 0:
      print('Deriv2:',x.size)
   if sym:
      B = _Deriv2sym(geom,x,A)
   else:
      B = _Deriv2non(x,A)
   return(B)

def _Deriv2non(x,A):
   ''' Calculate 2nd derivative for nonsymmetric case. '''
   nt = x.size
   B = np.empty((nt,nt))
   for i in range(nt):
      B[i,:] = x[i] - x
      B[i,i] = 1.0
      B[i,:] = 2.0*A[i,:]*(A[i,i] - 1.0/B[i,:]) #/#
      B[i,i] = 0.0
      B[i,i] = -np.sum(B[i,:])
   return(B)

def _Deriv2sym(geom,x,A):
   ''' Calculate 2nd derivative for symmetric case. '''
   nt = x.size
   B = np.empty((nt,nt))
   gx = 0.5*float(geom)
   xx = x*x
   for i in range(nt):
      B[i,:] = xx[i] - xx
      B[i,i] = 1.0
      B[i,:] = 2.0*A[i,:]*(A[i,i] + gx/x[i] - 2.0*x[i]/B[i,:])
      B[i,i] = 0.0
      B[i,i] = -np.sum(B[i,:])
   return(B)
   
def _Stiffness(geom,sym,x,wb,wq):
   ''' Calculate Stifffness matrix. '''
   nt = x.size
   if Debug > 0:
      print('Stiffness:',nt)
   A = _Deriv1(sym,x,wb)
   B = _Deriv2(geom,sym,x,wb)
   C = np.empty((nt,nt))
   for i in range(nt):
      C[i,:] = -wq[i]*B[i,:]
   C[nt-1,:] = C[nt-1,:] + A[nt-1,:]
   if not sym:
      C[0,:] = C[0,:] - A[0,:]
   return(C)

def _Lcoeff(x,xc,wb=None,sym=0):
   ''' Interpolating polynomials thru xc evaluated at x '''
   if sym:
      L = _Lcoefx(x*x,xc*xc,wb)
   else:
      L = _Lcoefx(x,xc,wb)
   return (L)

def _Lcoefx(x,xc,wb=None):
   ''' Interpolating polynomials thru xc evaluated at x '''
   #  Calculate values of the Lagrange interpolating polynomials:
   #     x  - points where values are desired (evaluation points)
   #     xc - the n interpolation points (collocation points)
   #     wb - barycentric weights
   #     L(i,j) - value jth interpolating polynomial at xx(i)
   #  example usage: given values yc at the n interpolation points
   #     y(x) = MatMul(L,yc) is the interpolant at x(:)
   nx = x.size
   nc = xc.size
   L = np.empty((nx,nc))
   if wb is None:
      wb = _WBcalc(xc)
   for k in range(nx):
      if np.min(abs(x[k]-xc[:])) > 0.0001:
         L[k,:] = np.prod(x[k]-xc)*wb/(x[k]-xc) #/#
      else:  # long calculation
         px = np.empty(nc)
         pp = np.empty(nc)
         for i in range(nc):
            px[:] = x[k] - xc
            px[i] = 1.0
            pp[i] = px.prod()
         L[k,:] = pp*wb
   return(L)

def _WBcalc(x):
   ''' Calculate barycentric weight from product formula. '''
   n = x.size
   wb = np.empty(n)
   px = np.empty(n)
   for i in range(n):
      px = x[i] - x
      px[i] = 1.0
      wb[i] = 1.0/px.prod()   #/#
   return(wb)
   
def _MassMatrix(oc,func=None,nextra=1):
   nt = oc.nt
   if func is None:
      if not oc.sym and oc.ptyp == Typ.Lobatto.value:
         D = _MassMatL(oc.x,oc.wb,oc.shft)
      elif not oc.sym and oc.ptyp == Typ.Gauss.value:
         D = _MassMatG(oc.x,oc.wb,oc.shft)
      else:
         D = np.zeros((nt,nt))
         for i in range(nt):
            D[i,i] = oc.wq[i]
   elif oc.ptyp <= Typ.Lobatto.value:
      if Debug > 0:
         print('_MassMatFn:',nt,nextra,oc.sym,oc.ptyp)
      D = _MassMatFn(oc,func,nextra)
   else:
      if Debug > 0:
         print('_MassMatrix: other',nt,oc.sym,oc.ptyp)
      D = np.zeros((nt,nt))
      fx = func(oc.x)
      for i in range(nt):
         D[i,i] = oc.wq[i]*fx[i]
   return(D)

def _MassMatFn(oc,func,nextra=1):
   ''' Calculate Galerkin or Moments mass matrix with f(x). '''
   # approximated using quadrature with nextra extra points
   sym = oc.sym;   nt = oc.nt;   n = oc.n
   shift = oc.shft
   Gauss = (oc.ptyp == Typ.Gauss.value)
   Lobatto = (oc.ptyp == Typ.Lobatto.value)
   nq = n + nextra   # interior points for quadrature
   oq = OrthColloc(nq,oc.ptyp,oc.geom,oc.shft)
   i1 = 1 - sym
   i2 = nq + i1
   Lq = oc.Lcoeff(oq.x[i1:i2])
   wq = func(oq.x)
   if Gauss and (not sym) and shift:
      wq[i1:i2] = wq[i1:i2]*(oq.wq[i1:i2])/(oq.x[i1:i2]*(1.0 - oq.x[i1:i2])) #/# moments shifted
   elif Gauss:
      wq[i1:i2] = wq[i1:i2]*(oq.wq[i1:i2])/(1.0 - (oq.x[i1:i2])**2)   #/#
   else:
      wq = wq*(oq.wq)     # Galerkin
   D = np.empty((nt,nt))
   for i in range(nt):
      for j in range(i,nt):
         D[i,j] = np.sum(wq[i1:i2]*Lq[:,i]*Lq[:,j])
         D[j,i] = D[i,j]
   if Gauss and (not sym) and shift:
      for i in range(nt):
         D[i,:] = D[i,:]*(oc.x[i])*(1.0 - oc.x[i])
   elif Gauss:
      for i in range(nt):
         D[i,:] = D[i,:]*(1.0 - (oc.x[i])**2)
   elif Lobatto and sym:
      D[n,n] = D[n,n] + wq[nq]    # end pt. weights
   elif Lobatto:
      D[0,0] = D[0,0] + wq[0]
      D[n+1,n+1] = D[n+1,n+1] + wq[nq+1]
   return(D)

def _MassMatG(x,wb,shift):
   ''' Calculate moments-Gauss mass matrix f(x) = 1. '''
   n = x.size - 2
   if Debug > 0:
      print('_MassMatG:',n)
   D = np.zeros((n+2,n+2))
   Cn = -1.0/(2.0*n + 1.0)  #/#
   Cd = -2.0*n
   if shift:
      wx = Cn*x*(1.0 - x)*wb
   else:
      wx = 2.0*Cn*(1.0 - x*x)*wb
   for i in range(n+1):
      D[i,:] = wx[i]*wb
      D[i,i] = Cd*D[i,i]
   return(D)

def _MassMatL(x,wb,shift):
   ''' Calculate Galerkin-Lobatto mass matrix for f(x) = 1. '''
   n = x.size - 2
   D = np.empty((n+2,n+2))
   Cn = float(-(n+1))/float((2*n+3)*(n+2)) #/#
   Cn = Cn if shift else Cn*8.0
   Cd = float(-2*(n+1))
   for i in range(n+2):
      D[i,:] = Cn*wb[i]*wb
      D[i,i] = Cd*D[i,i]
   if Debug > 0:
      print('_MassMatL:',n,Cn,Cd,float(-(n+1))/float((2*n+3)*(n+2)),file=ap.file()) #/#
      ap.vectorprint('MassMatL: wb',wb)
   return(D)

# ------------- OrthColloc Class for calculating basic data -----------
class OrthColloc:
   ' Class for Orthogonal Collocation & MWR basic data '
   def __init__(self,n,TType=Typ.Lobatto.value,Geometry=0,Shift=True):
      self.n = n
      self.ptyp = TType
      self.geom = Geometry
      self.sym  = max(0,min(Geometry,1))
      if self.sym:
         self.shft = True  # symmetric always shifted
      else:
         self.shft = Shift
      self.nt = self.n + 2 - self.sym
      self.abeta = jp.jac_abeta(TType,Geometry)
      w,wbfac = jp.jac_quadrature(n,self.abeta,Geometry)
      if self.shft and self.geom==0:
         w,wbfac = _Wshift(n,w,wbfac)
      self.x  = w[0,:]
      self.th = w[1,:]
      self.wb = w[2,:]
      self.wq = w[3,:]
      self.wbfac = wbfac
      if Debug > 2:
         print('n,type,geometry,shift:',n,TType,Geometry,Shift,file=ap.file())
         ap.arrayprint("x,theta,wb,wq "+str(wbfac)+":",np.transpose(w),fmtf='%.15f')

   def Xpoints(self):
      ''' Return points or roots of polynomial and end points. '''
      return(self.x)

   def WeightQ(self):
      ''' Return quadrature weights. '''
      return(self.wq)

   def WeightB(self):
      ''' Return barycentric weights. '''
      return(self.wb)
      
   def WBfac(self):
      ''' Return scaling factor for barycentric weights. '''
      return(self.wbfac)
      
   def Deriv1(self):
      ''' Calculate first derivative matrix. '''
      A = _Deriv1(self.sym,self.x,self.wb)
      return(A)

   def Deriv1odd(self):
      ''' Calculate first derivative matrix for odd (antisymmetric) function. '''
      An = _Deriv1odd(self.geom,self.sym,self.x,self.wb)
      return(An)

   def Deriv2(self):
      ''' Calculate second derivative or Laplacian matrix. '''
      B = _Deriv2(self.geom,self.sym,self.x,self.wb)
      return(B)

   def Stiffness(self):
      ''' Calculate second derivative or Laplacian matrix. '''
      C = _Stiffness(self.geom,self.sym,self.x,self.wb,self.wq)
      return(C)

   def Lcoeff(self,x):
      L =_Lcoeff(x,self.x,(self.wb)*(self.wbfac),self.sym)
      return(L)

   def MassMatrix(self,func=None,nextra=1):
      D = _MassMatrix(self,func,nextra)
      return(D)

   def info(self):
      txtt = str(txt_type[self.ptyp]).strip("[']")
      txtg = str(txt_geom[self.geom]).strip("[']")
      print('Calculations with',txtt,'points and',txtg,'geometry, n =',self.n)


