import numpy as np
import numpy.linalg as linalg
import math
import arrayprint as ap
import occ
from occ import OrthColloc
sym = 1
shift = True
def Example02(n,Thx,Bi):
   ''' Solve Example BVP, nonlinear source, planar, cylinder or spherical'''
   global ptyp,geom,korder,Kads
   maxiter = 8
   n1 = n + 1
   Generalized = True   # to supply value of generalized parameter
   Thiele = Thx
   if Generalized:
      Thiele = Thx*geom/GenThiele()    #/#
   Th2 = (Thiele)**2
   Dirichlet = Bi <= 0 or Bi > 1.0e+8
   oc = OrthColloc(n,ptyp,geom)
   x = oc.Xpoints()
   w = oc.WeightQ()
   Ac = oc.Deriv1()
   Cc = oc.Stiffness()
   ReportBasic(Thiele,Bi,x,w,Ac,Cc)
   wt = Th2*w
   y = np.empty(n1)
   dy = np.zeros(n1)
   flux = np.zeros(3)
   A = np.empty((n1,n1))
   Aflux = np.empty((2,n1))
   # Solve with collocation - symmetric matrix formulation
   y = Guess(x)
   for iter in range(maxiter):
      r,dr = RxnRate(y)
      A[:,:] = Cc
      for i in range(n1):
         A[i,i] -= wt[i]*dr[i]
      b = wt*r - Cc.dot(y)
      if Dirichlet:
         dy[:n] = linalg.solve(A[:n,:n],b[:n])
      else:
         A[n,n] += Bi
         b[n] -= Bi*y[n]
         dy = linalg.solve(A,b)
      y = y + dy
   r,dr = RxnRate(y)
   flux[0] = np.sum(w*r)*geom
   flux[1] = (w[n]*r[n] - Cc[n,:].dot(y)/Th2)*(geom) #/# -
   flux[2] = -Ac[n,:].dot(y)*(geom/Th2)  #/#
   ReportFinal(oc,x,y,r,flux,Th2)

def ReportFinal(oc,x,y,r,flux,Th2):
   ''' Report final results. '''
   sfmt = '%16.9g'
   print('\nFinal Results:',file=ap.file())
   ap.arrayprint('x, y, r:',x,y,r,fmtf=sfmt)
   ap.vectorprint('Normalized flux - integral, stiffness, derivative:',flux,fmtf=sfmt)
   nv = 22
   pif = 0.5*occ.pi_const/(nv-1)  #/#
   xv = np.arange(nv-1,-1,-1)
   xv = np.cos(xv*pif)
   Bc = oc.Deriv2()
   Lc = oc.Lcoeff(xv)
   yv = Lc.dot(y)
   rv,drv = RxnRate(yv)
   v = np.empty((nv,3))
   yv = Lc.dot(y)
   d2c = Bc.dot(y)
   v[:,0] = Lc.dot(d2c)       # second derivative at xv
   v[:,1] = Th2*rv            # rate at xv
   v[:,2] = v[:,0] + v[:,1]   # residual
   ap.arrayprint('x, y, d2y, r(y), Resid',xv,yv,v,fmtf=sfmt)

def ReportBasic(Thiele,Bi,x,w,A,C):
   ''' Report basic data for problem. '''
   global ptyp,geom,korder,Kads
   sfmt = '%16.9g'
   txtt = str(occ.txt_type[ptyp]).strip("[']")
   txtg = str(occ.txt_geom[geom]).strip("[']")
   n = x.size - 1
   Thgen = GenThiele()*Thiele/geom  #/#
   print('\nSolution of BVP with',txtt,'points and n =',n,file=ap.file())
   print('Symmetric problem in '+txtg+' geometry',file=ap.file())
   print('Thiele parameter, nominal, generalized =',Thiele,Thgen,' and Bi =',Bi,file=ap.file())
   print('Reaction order =',korder,', denominator (adsorption) term Ka =',Kads,file=ap.file())
   ap.arrayprint("\nPoints, Quad. Wts.:",x,w,fmtf='%18.15f')
   ap.arrayprint("\nFirst Derivative: ",A,fmtf=sfmt)
   ap.arrayprint("\nStiffness Matrix: ",C,fmtf=sfmt)

def RxnRate(y):
   ''' Calculate reaction rate given conversion y. '''
   global ptyp,geom,korder,Kads,yguess
   nt = y.size
   r = np.empty(nt)
   dr = np.empty(nt)
   drext = 0.0
   if korder == 1:
      drext = -1.0/(1.0 - Kads)**2  #/#
   for i in range(nt):
      if y[i] < 1.0:
         d = 1.0/(1.0 - Kads*y[i])  #/#
         r[i] = (1.0 - y[i])**korder
         dr[i]= -korder*(1.0 - y[i])**(korder-1) + 2.0*r[i]*Kads*d
         r[i] = r[i]*(d*d)
         dr[i] = dr[i]*(d*d)
      else:
         dr[i] = drext
         r[i] = dr[i]*(y[i] - 1.0)
   return(r,dr)

def GenThiele():
   ''' Calculate T*(g+1)/T General Thiele parameter ratio to nominal value. '''
   # where T* is the generalized value, g = 1,2,3 for planar,cylindrical,spherical
   global ptyp,geom,korder,Kads,yguess
   if Kads > 0 and korder == 1:        # 1st order, Kads > 0
      Thfac = math.sqrt(-2.0*(Kads + math.log(1.0 - Kads)))
      Thfac = Kads/Thfac   #/#
   elif Kads == 0:                     # kth order Kads = 0
      Thfac = 0.5*math.sqrt((korder+1.0)*2.0)
   else:
      Thfac = 1.0                      # others?
   return(Thfac)

def Guess(x):
   global ptyp,geom,korder,Kads,yguess
   xg = 0.50
   fac = 1.0/(1.0 - xg*xg) #/#
   y = yguess*(1.0 - x*x)*fac
   y = np.minimum(y,0.999)
   return(y)

def main():
   ''' Main program to read data and launch solution routine. '''
   global ptyp,geom,korder,Kads,yguess
   c = input('Enter g/l for Gauss/Lobatto points >')
   ptyp = ('gl').find(c) + 1
   ptyp = max(1,min(ptyp,2))
   c = input('Enter 1, 2 or 3 for Planar, Cylindrical or Spherical Geometry > ')
   geom = int(c)
   c = input('Enter n values separated by whitespace: first last delta >')
   nums = [int(i) for i in c.split(' ')]
   c = input('Enter 0 for Dirichlet or Bi for flux b.c. > ')
   Bi = float(c)
   c = input('Reaction order and denominator Ka term >')
   parms = [float(x) for x in c.split(' ')]
   korder = parms[0]
   Kads = parms[1]
   c = input('Thiele parameter: first,last,delta >')
   Thp = [float(x) for x in c.split(' ')]
   c = input('Enter initial guess > ')
   yguess = float(c)
   c = input('Enter root of file name >')
   ap.fopen(c+'.dat')
   for th in np.arange(Thp[0],Thp[1]+Thp[2],Thp[2]):
      for n in range(nums[0],nums[1]+1,nums[2]):
         Example02(n,th,Bi)

main()

