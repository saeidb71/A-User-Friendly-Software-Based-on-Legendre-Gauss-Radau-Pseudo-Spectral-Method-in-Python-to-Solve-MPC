import numpy as np
import numpy.linalg as linalg
import scipy.linalg as la
import math
import arrayprint as ap
import occ
import pdb
Bcolloc = False   # use boundary collocation?
sfmt = '%14.8f'
sfma = '%20.12f'
def FallingFilm(ptyp,zmeth,zmax,nstep,nmax):
   ''' Solve falling film problem, continuous time '''
   # ----------------------------------------------------------------
   # This code solves mass transfer from falling film problem as described
   # in B,S&L p. 538 and solved by Villadsen and Michelsen [1978, Ch. 4]
   # and Finlayson [1972, 2014, pp 41, 58]
   # The code is set up for solution by collocation with
   #    Gauss, Lobatto, Chebyshev, Radau (L or R) points
   # It can also treat the no flux boundary condition either naturally or
   # by boundary collocation (set BColloc = .true. above)
   # The system of linear algebraic equations is solved using various
   # stepping methods
   # xc,wc,Ac,Cc - points, weights, 1st Deriv, Stiffness
   # A - matrix approximation
   # wx = wx*(1-x^2) - diagonal mass matrix
   # ys[:,i] - 0 yavg, 1 flux, 2 Sh at each step
   # yc[:,i] - y at point i for each step
   # r,y,dy - rhs, y & delt y at interior nodes, current step
   # ----------------------------------------------------------------
   global n,ys,yc,C0,Cn1,wx,txt_method
   n1 = n + 1
   oc = occ.OrthColloc(n,ptyp)
   xc = oc.Xpoints()    # points
   wc = oc.WeightQ()    # weights
   Cc = oc.Stiffness()  # stiffness
   wx = wc[1:n1]*(1.0 - xc[1:n1]**2)
   dz = zmax/float(nmax-1)
   dzi = 1.0/dz;
   wdz = -dz/wx   #/# inverse value
   ys = np.zeros((nmax,3),dtype=float)
   yc = np.zeros((nmax,n1),dtype=float)
   y = np.ones((n),dtype=float)
   dy = np.empty((n),dtype=float)
   dy1 = np.empty((n),dtype=float)
   r = np.empty((n),dtype=float)
   z = 0.0
   A = np.copy(Cc[1:n+1,1:n+1])
   if(Bcolloc):
      Ac = oc.Deriv1()
      Cc[n1,:] = Ac[n1,:]           # to implement boundary collocation
   Cn1 = -Cc[n1,1:n1]/Cc[n1,n1]     #/# y[n+1] coefficients
   C0 = Cc[0,1:n1] + Cc[0,n1]*Cn1   #/# flux coefficients/
   for j in range(0,n): # calculate half & copy
      A[j:,j] = A[j:,j] + Cc[j+1:n+1,n+1]*Cn1[j]
      A[j,j+1:] = A[j+1:,j]
   txt_m = str(txt_method[zmeth-1]).strip("[']")
   txt_t = str(occ.txt_type[ptyp]).strip("[']")
   print('Falling Film Problem, '+txt_t+' points, '+txt_m+', n ={:3d}, dz ={:8.6f}'.format(n,dz))
   print('\nFalling Film Problem, '+txt_t+' points, '+txt_m+', n ={:3d}, dz ={:8.6f}'.format(n,dz),file=ap.file())
   ap.arrayprint('x,w,Stiffness',xc,wc,Cc,fmtf=sfmt)
   ap.arrayprint('wx,A',wx,A,fmtf=sfmt)
   yx = CalcFlux(y,0)
   if(zmeth == 1):      # forward Euler
      for k in range(1,nmax):
         dy = A.dot(y)*wdz
         y = y + dy
         yx = CalcFlux(y,k)
         if yx[1] <= 0.0: break   # unstable
   elif(zmeth == 2):    # 2nd Order R-K, Heun, improved Euler
      for k in range(1,nmax):
         dy = A.dot(y)*wdz
         dy1 = A.dot(y+dy)*wdz
         y = y + (dy + dy1)*0.5
         yx = CalcFlux(y,k)
         if yx[1] < 0.0 or yx[0] > 1.1: break   # unstable
   elif(zmeth == 3):    # trapezoidal rule
      Alu = BuildF(A,wx*(2.0*dzi)) # factor matrix
      for k in range(1,nmax):
         r = (-2.0)*A.dot(y)
         dy = la.lu_solve(Alu,r) # forward & back solve on rhs
         y = y + dy
         yx = CalcFlux(y,k)
   elif(zmeth >= 4 and zmeth <= 6): # Backward difference methods
      nbk = zmeth - 4
      if nbk > 0:
         mdrk = np.array([1,3])[nbk-1] #Dirk method for starting
         Dirk(A,mdrk,dz,nbk+1)         # starting
      BackDiff(A,dz,nbk,nmax)
   elif(zmeth >= 7):    # Diagonally Implicit R-K
      Dirk(A,zmeth-6,dz,nmax)
   ReportFinal(nmax,dz)

def ReportFinal(nmax,dz):
   ''' Report final results '''
   global n,ys,yc,C0,Cn1,wx
   tk = dz*np.arange(nmax)
   print('   step \t z \t yavg \t flux \t Sh ',file=ap.file())
   for k in range(nmax):
      if ys[k,1] <= 0.0: break   # unstable
      print('{:6d}\t{:8.6f}'.format(k,k*dz),('\t{:12.10f}'*3).format(ys[k,0],ys[k,1],ys[k,2]),file=ap.file())
      if ys[k,1] < 0.0 or ys[k,0] > 1.1: break   # unstable
   ap.arrayprint('\t z \ty(1:nc+1):',tk,yc,fmtf=sfmt)

def BuildF(A,w):
   Aw = np.diag(w)
   Aw = Aw + A
   Alu = la.lu_factor(Aw)
   return(Alu)

def CalcFlux(y,k):
   global n,ys,yc,wx,C0,Cn1
   yx = np.empty((3),dtype=float)
   yc[k,:n] = y                  # save values
   yc[k,n] = sum(Cn1*y)        # y(n+1)
   yx[0] = 1.5*sum(wx*y)         # average y
   yx[1] = -sum(C0*y)            # flux
   yx[2] = yx[1]/yx[0]           #/# Sh no.
   ys[k,:] = yx
   return(yx)

def BackDiff(A,dz,nbk,nmax):
   '''solution with Backward difference methods'''
   global n,ys,yc,wx,C0,Cn1
   bkx = np.array([[1.0,-1.0,0.0,0.0],
                   [1.5,-2.0,0.5,0.0],
                   [11/6,-3.0,1.5,-1/3]])
   ak = bkx[nbk,:nbk+2]/dz
   ak[2:] = ak[2:]/ak[1]
   wi = -wx*ak[1]
   #print('bk',bkx[nbk,:nbk+2],file=ap.file())
   #print('ak',ak,file=ap.file())
   #ap.arrayprint('Start:',yc[:nbk+1,:])
   Alu = BuildF(A,wx*ak[0]) # factor matrix
   for k in range(nbk+1,nmax):
      r = yc[k-1,:n]
      for i in range(2,nbk+2):
         r = r + ak[i]*yc[k-i,:n]
      r = r*wi
      y = la.lu_solve(Alu,r) # forward & back solve on rhs
      yx = CalcFlux(y,k)

def Dirk(A,method,dz,nmax):   #
   ''' solution with Diagonally Implicit Runge-Kutta'''
   global n,ys,yc,wx,C0,Cn1
   ar,br,cr = DirkParm(method-1)
   ns = cr.size
   explicit1 = np.amax(ar[0,:]) == 0.0
   #print('method, ns, Exp =',method,ns,explicit1,file=ap.file())
   #ap.arrayprint('Butcher c,a:',cr,ar,fmtf=sfmt)
   #ap.vectorprint('Butcher b:',br,fmtf=sfmt)
   ar = ar*dz; br = br*dz; cr = cr*dz
   ad = ar[ns-1,ns-1]; adi = 1.0/ad    #/#
   wa = wx*adi;  wi = (-1.0)/wx        #/#
   y = np.ones((n),dtype=float)
   f = np.zeros((n,ns),dtype=float)
   yx = CalcFlux(y,0)
   Alu = BuildF(A,wa) # factor matrix
   for k in range(1,nmax):
      i = 0
      y0 = y
      if explicit1:  # explicit first stage
         f[:,i] = A.dot(y)*wi
         i = 1
      for i in range(i,ns): # loop over implicit stages
         y = y0 + f[:,:i].dot(ar[i,:i])
         f[:,i] = (la.lu_solve(Alu,wa*y) - y)*adi   # forward & back solve on rhs
      y = y0 + f.dot(br)
      yx = CalcFlux(y,k)

def DirkParm(method):
   '''returns Butcher tableu'''
   # General framework for DIRK methods
   #   DirkParm - returns the Butcher tableu
   #   methods:
   #     1. SDIRK - 2nd order, 2 stage, 2 implicit stages
   #     2. EDIRK - 2nd order, 3 stage, 2 implicit stages
   #     3. SDIRK - 3rd order, 3 stage, 3 implicit stages
   #     4. EDIRK - 3rd order, 4 stage, 3 implicit stages
   #     5. EDIRK - 3rd order, 5 stage, 4 implicit stages
   NDIRK = 5   # number methods stored
   sqr2 = np.sqrt(2.0)
   E3g = 1.0 - (1.0)/sqr2  #/#
   # Parameters for ESDIRK3(2)4L[2]SA (3 implicit stages)
   E4g = 0.43586652150845899941601945
   E4c3 = 0.6
   E4xx = E4c3-2.0*E4g
   E4a32 = E4c3*E4xx/(4.0*E4g)   #/#
   E4b2 = (-2.0 + 3.0*E4c3 + 6.0*E4g*(1.0 - E4c3))/(12.0*E4g*E4xx) #/#
   E4b3 =(1.0 - 6.0*E4g + 6.0*E4g**2)/(3.0*E4c3*E4xx) #/#
   # Parameters for ESDIRK3(2)5L[2]SA (4 implicit stages)
   E5g = 9.0/40.0 #/#
   E5c = np.array([0.0,2.0*E5g,(2.0 + sqr2)*E5g,3.0/5.0,1.0])
   E5a43 = (-7.0)/((1.0+sqr2)*40.0)
   E5a54 = (5827.0)/(7560.0)
   E5a53 =(-2374.0)*(1.0+2.0*sqr2)/(2835.0*(5.0+3.0*sqr2))
   Diag = np.array([E3g,E3g,E4g,E4g,E5g])
   first = np.array([0,1,0,1,1])
   nstage = np.array([2,3,3,4,5])
   # -------------------------------
   mth = max(0,min(method,NDIRK-1))
   ns = nstage[mth]
   a = np.zeros((ns,ns),dtype=float)
   b = np.zeros((ns),dtype=float)
   c = np.zeros((ns),dtype=float)
   np.fill_diagonal(a[first[mth]:,first[mth]:],Diag[mth])
   if mth == 0:
      a[1,0] = 1.0 - E3g
   elif mth == 1:
      a[2,0] = (0.25 - 0.50*E3g)/E3g   #/#
      a[1,0] = E3g
      a[2,1] = a[2,0]
   elif mth == 2:
      a[2,1] = 1.25 + E4g*((1.5)*E4g - 5.0)
      a[2,0] = 1.0 - a[2,1] - a[2,2]
      a[1,0] = (1.0 - E4g)*(0.5)
   elif mth == 3:
      a[1,0] = E4g
      a[2,1] = E4a32
      a[2,0] = E4c3 - np.sum(a[2,1:3])
      a[3,2] = E4b3
      a[3,1] = E4b2
      a[3,0] = 1.0 - np.sum(a[3,1:4])
   elif mth == 4:
      a[1,0] = E5g
      a[3,2] = E5a43
      a[4,3] = E5a54
      a[4,2] = E5a53
      for i in range(2,5):
         a[i,0:2] = (E5c[i] - np.sum(a[i,2:i+1]))*0.5
   b = a[ns-1,:]  # stiffly accurate assumed (make exception if not true)
   for i in range(ns):
      c[i] = np.sum(a[i,:])   # c = row sum(a)

   return(a,b,c)

def main():
   ''' Main program to read data, then solve falling film problem '''
   global n,txt_method
   txt_method = (['Forward Euler'],['Improved Euler'],['Trapezoid Rule'],['Backward Euler'],  \
      ['2nd Order Bkwd'],['3rd Order Bkwd'],['SDIRK 2nd Ord.'],['EDIRK 2nd Ord.'],['SDIRK3 3 stage'], \
      ['EDIRK3 4 stage'],['EDIRK3 5 stage'])
   c = input('Enter g/l/r/x for Gauss/Lobatto/Radau(R/L) >')
   ptyp = ('gl!rx').find(c) + 1
   n = int(input('Enter value of n >'))
   print('Enter the Time Integration Method from:')
   i = 0
   for item in txt_method:
      i += 1
      txt_m = str(item).strip("[']")
      print('{:3d} {}'.format(i,txt_m))
   zmeth = int(input('>'))
   c = input('Enter time step sizes <= 10 >')
   dzsteps = [float(i) for i in c.split(' ')]
   c = input('Enter root of file name >')
   ap.fopen(c+'.dat')
   zmax = 0.4
   z01 =  0.1      # hit this z exactly (comparison point)
   for dz in dzsteps:
      nstep = int((z01 + dz*0.10)/dz)
      dz = z01/float(nstep)
      FallingFilm(ptyp,zmeth,zmax,nstep,int(zmax/z01)*nstep+1) #/#
main()

