import numpy as np
import numpy.linalg as linalg
from scipy.linalg import eigh
import math
import arrayprint as ap
import occ
sfmt = '%14.8f'
sfma = '%14.8f'
def FallingFilm(n,ptyp):
   ''' Solve falling film problem, continuous time '''
   # ----------------------------------------------------------------
   # This code solves mass transfer from falling film problem as described
   # in B,S&L p. 538 and solved by Villadsen and Michelsen [1978, Ch. 4]
   # and Finlayson [1972, 2014, pp 41 & 58]
   # The code is set up for solution by collocation with
   #    Gauss, Lobatto, Radau (L or R) points
   # Chebyshev points are not implemented in this release.
   # It can also treat the no flux boundary condition either naturally or
   # by boundary collocation (set BColloc = .true. above)
   # with minor modifications it can be used with the Galerkin method.
   # The system of linear algebraic equations is solved analytically
   # ----------------------------------------------------------------
   oc = occ.OrthColloc(n,ptyp)
   xc = oc.Xpoints()    # points
   wc = oc.WeightQ()    # weights
   Cc = oc.Stiffness()  # stiffness
   wx = wc[1:n+1]*(1.0 - xc[1:n+1]**2)
   A = np.copy(Cc[1:n+1,1:n+1])
   yk = np.empty((n+1,n),dtype=float)
   for j in range(0,n): # calculate half & copy
      A[j:,j] = A[j:,j] - Cc[j+1:n+1,n+1]*Cc[n+1,j+1]/Cc[n+1,n+1] #/#
      A[j,j+1:] = A[j+1:,j]
   B = np.diag(wx)
   xlx,v = eigh(A,B)    # generalized eigenvalue problem symmetric (Hermitian)
   c = linalg.solve(v,np.ones(n,dtype=float))   # initial conditions
   yk[:n,:] = v*c             # y interior
   yk[n,:] = -Cc[n+1,1:n+1].dot(yk[:n,:])/Cc[n+1,n+1] #/# boundary at x = 1
   s = 1.5*wx.dot(yk[:n,:])   # mixing cup average
   sf = -Cc[0,1:].dot(yk)     # boundary flux, x = 0
   txtt = str(occ.txt_type[ptyp]).strip("[']")
   print('\nFalling Film Problem, '+txtt+' points, n =',n,file=ap.file())
   ap.arrayprint('x,w,Stiffness',xc,wc,Cc,fmtf=sfma)
   ap.arrayprint('w,A',wx,A,fmtf=sfma)
   ap.arrayprint('xlx,v,c - eigenvalues, vectors, i.c.',xlx,v,c,fmtf=sfma)
   ap.arrayprint('yk - coefficients of y',yk,fmtf=sfma)
   ap.arrayprint('s,sf - yavg & flux coefficients',s,sf,xlx*s/1.5,fmtf=sfma)   #/#
   ReportFinal(oc,xc,xlx,yk,s,sf)

def ReportFinal(oc,xc,xlx,yk,s,sf):
   ''' Calculate solution and report results '''
   # z - y(x,z) value for profile calculations
   zx = np.array([0.0001,0.001,0.002,0.004,0.01,0.02,0.04,0.1,0.2,0.4])
   nz = zx.size
   nv = 41
   zv = 10**(-4.0 + 0.1*np.arange(nv))
   yc = np.zeros((oc.n+2,nz),dtype=float)
   yv = np.empty((nv,3),dtype=float)
   for i in range(0,nz):
      expk = np.exp(-xlx*zx[i])
      yc[1:,i] = yk.dot(expk)
   for i in range(0,nv):
      expk = np.exp(-xlx*zv[i])  # exponential
      yv[i,0] = np.sum(s*expk)   # yavg
      yv[i,1] = np.sum(sf*expk)  # flux
   yv[:,2] = yv[:,1]/yv[:,0]     #/# Sh
   print('\nFinal Results:',file=ap.file())
   ap.vectorprint('values at (columns)',zx,fmtf=sfmt)
   ap.arrayprint('x, y(zx)',xc,yc,fmtf=sfmt)
   ap.arrayprint('z, yavg, flux, Sh:',zv,yv,fmtf=sfmt)

def main():
   ''' Main program to read data, then solve falling film problem '''
   c = input('Enter g/l/r/x for Gauss/Lobatto/Radau(R/L) >')
   ptyp = ('gl?rx').find(c) + 1  # omit Chebyshev
   c = input('Enter n values separated by whitespace: first last delta >')
   nums = [int(i) for i in c.split(' ')]
   c = input('Enter root of file name >')
   ap.fopen(c+'.dat')
   for n in range(nums[0],nums[1]+1,nums[2]):
      FallingFilm(n,ptyp)
main()

