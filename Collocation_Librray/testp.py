import numpy as np
import math
import arrayprint as ap
import occ
from occ import OrthColloc
def test(n,testi):
   global ptyp,geom,sym,shift
   shift = True   # [-1,1] or [0,1], i.e. shifted
   sym = min(1,geom)
   oc = OrthColloc(n,ptyp,geom,shift)
   oc.info()
   x = oc.Xpoints()
   w = oc.WeightQ()
   wb = oc.WeightB()
   A = oc.Deriv1()
   if geom > occ.Geom.Nonsymmetric.value:
      An = oc.Deriv1odd()
   B = oc.Deriv2()
   C = oc.Stiffness()
   D = oc.MassMatrix()
   fmass = MassFunc(x)
   Df = oc.MassMatrix(MassFunc,2)
   sfmt = '%14.7g'
   txtt = str(occ.txt_type[ptyp]).strip("[']")
   txtg = str(occ.txt_geom[geom]).strip("[']")
   print('\nCollocation Parameters for',txtg,'and',txtt,'points for n =',n,file=ap.file())
   ap.arrayprint("\nPoints, Quad. Wts., Bar. Wts, factor = "+str(oc.WBfac())+":",x,w,wb,fmtf='%18.15f')
   ap.arrayprint("\nFirst Derivative: ",A,fmtf=sfmt)
   if geom > occ.Geom.Nonsymmetric.value:
      ap.arrayprint("\nFirst Derivative of odd function: ",An,fmtf=sfmt)
   ap.arrayprint("\nSecond Derivative: ",B,fmtf=sfmt)
   ap.arrayprint("\nStiffness Matrix: ",C,fmtf=sfmt)
   if D is not None:
      ap.arrayprint("\nMass Matrix, f(x) = 1: ",D,fmtf=sfmt)
      ap.vectorprint("\nMass Matrix f(x) Values:",fmass)
      ap.arrayprint("\nMass Matrix, with f(x): ",Df,fmtf=sfmt)
   if testi:
      InterpCalcs(oc,x,w,A)

def InterpCalcs(oc,x,w,A):
   ''' Demonstrate approximate integration, differentiation & interpolation. '''
   global ptyp,geom,sym,shift
   nt = x.size
   n = nt - 2 + sym
   ii = 0
   ni = 21
   if sym:
      ii = ni//2
   xi = np.arange(ii,ni,dtype=float)
   xi = -np.cos(occ.pi_const*xi/(ni-1.0))  #/#
   if shift and not sym:
      xi = (xi + 1.0)*0.5
   ni = xi.size
   L = oc.Lcoeff(xi)
   df = np.empty((3,nt))
   fi = np.empty((3,ni))
   use_exp = True #<< change manually for alternate function
   if use_exp:
      f,df[0,:],fint = ExpFunc(x)      # f(x) = exp(-5x^2)
      fi[0,:],dfi,fint = ExpFunc(xi)
   else:
      f,df[0,:],fint = Ex01Func(x)     # f(x) = cosh(5x)/cosh(5)
      fi[0,:],dfi,fint = Ex01Func(xi)  # integral invalid cylindrical or spherical
   fintq = np.sum(w*f)  # approximate integration
   df[1,:] = A.dot(f)   # approximate derivative at coll. pts
   fi[1,:] = L.dot(f)  # interpolated function
   finterr = 1.0 - fintq/fint  #/#
   df[2,:] = np.abs(df[1,:] - df[0,:])
   fi[2,:] = np.abs(fi[1,:] - fi[0,:])
   nlog = n*math.log10(n)
   ap.arrayprint('\nLagrange Interpolating Polynomials:',xi,L,fmtf='%12.8f')
   print('exp(-5*x^2) n (n)log(n) integral exact, approx, error:',file=ap.file())
   print('Integral n =',n,nlog,fint,fintq,finterr,file=ap.file())
   ap.arrayprint('x, df/dx, df/dx num. error, f:',x,df.transpose(),f,fmtf='%12.8f')
   ap.arrayprint('x, f, f-interp., error, deriv.:',xi,fi.transpose(),dfi,fmtf='%12.8f')

def MassFunc(xx):
   global ptyp,geom,sym,shift
   if sym:
      f = 1.0 - xx*xx            # Graetz problem
   else:
      x = xx
      if not shift:
         x = (xx + 1.0)*0.5
      a0 = 0.20
      a1 = (1.0 - a0)*4.0
      f = a0 + a1*(1.5 - x)*x*x  # BVP with cubic coef., 0.2 to 1.8
   return(f)

def ExpFunc(x):
   ''' exp(-5x^2) test function, value, derivative and integral '''
   global ptyp,geom,sym,shift
   import scipy.special as sc
   import math
   aexp = -5.0
   #if not shift:
   #   x = (xx + 1.0)*0.5
   #else:
   #   x = xx
   f = np.exp(aexp*x*x)    # function values
   df = (2.0)*aexp*x*f     # derivative values
   if geom == occ.Geom.Cylindrical.value: # integrals with dx, xdx & x^2dx
      fint = (-0.5/aexp)*(1.0 - np.exp(aexp))   #/#
   elif geom == occ.Geom.Spherical.value:
      fint = (0.25)*math.sqrt(-occ.pi_const/aexp**3)*sc.erf(math.sqrt(-aexp)) \
      + (0.5/aexp)*np.exp(aexp)
   else:
      fint = sc.erf(math.sqrt(-aexp))*math.sqrt(-occ.pi_const/aexp)*(0.5)   #/#
   #if not shift:
   #   fint = fint*2.0
   #   df = df*0.5
   return(f,df,fint)

def Ex01Func(xx):
   ''' cosh(5x)/cosh(5) - value, derivative & integral '''
   # integral not correct for cylindrical & spherical
   global ptyp,geom,sym,shift
   parm = 5.0
   if shift:
      x = xx*2.0 - 1.0
      dfac = 2.0
   else:
      x = xx
      dfac = 1.0
   div = 1.0/math.cosh(parm)              #/#
   f = np.cosh(parm*x)*div                #   function values
   df = np.sinh(parm*x)*div*parm*dfac     #   derivative values
   fint = 2.0*math.tanh(parm)/(parm*dfac) #/# integral
   return(f,df,fint)

def main():
   global ptyp,geom,sym,shift
   c = input('Enter: 0/1/2/3 for Nonsymmetric, Planar, Cylindrical or Spherical geometry >')
   geom = int(c)
   geom = max(0,min(3,geom))
   sym = min(1,geom)
   if sym:
      Shift = True   # symmetric always shifted
      c = input('Enter: g/l for Gauss/Lobatto >')
   else:
      c = input('Enter g/l/r/x for Gauss/Lobatto/Radau(R/L) >')
   ptyp = ('gl#rx').find(c) + 1
   c = input('Enter n values separated by whitespace: first last delta >')
   nums = [int(i) for i in c.split(' ')]
   c = input('Enter root of file name >')
   print('file',c+'.dat')
   ap.fopen(c+'.dat')
   c = input('Demonstrate integration, differentiation & interpolation y/n? >')
   testi = c == 'y'
   for n in range(nums[0],nums[1]+1,nums[2]):
      test(n,testi)
main()

