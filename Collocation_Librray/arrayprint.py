global fn
fmtdefault = '%.6e '
import numpy as np
def fopen(fname):
   global fn
   fn = open(fname,'w')

def fclose():
   fn.close()

def file():
   global fn
   return(fn)

def vectorprint(title,x,nl=20,fmtf=None):
   global fn
   n = np.size(x)
   if fmtf is None:
      fmtf = fmtdefault
   print(title,file=fn)
   for i0 in range(0,n,nl):
      i1 = min(i0+nl,n)
      np.savetxt(fn,x[i0:i1],fmtf,delimiter='',newline='\t') # normally a column
      print('',file=fn)

def arrayprint(title,a0,a1=None,a2=None,a3=None,fmtf=None):
   global fn
   if fmtf is None:
      fmtf = fmtdefault
   ns = np.zeros((4,2),dtype=int)   # shapes
   nc = np.zeros(5,dtype=int)       # columns
   a0vec = False
   a1vec = False
   a2vec = False
   a3vec = False
   nshape = a0.shape
   #print('shape a0:',nshape,file=fn)
   # count columns
   if len(nshape) == 1:
      a0vec = True
      a0.shape = (nshape[0],1)
   ns[0,:] = np.shape(a0)
   if a1 is not None:
      nshape = np.shape(a1)
      #print('shape a1:',nshape,file=fn)
      if len(nshape) == 1:
         a1vec = True
         a1.shape = (nshape[0],1)
      ns[1,:] = np.shape(a1)
   if a2 is not None:
      nshape = np.shape(a2)
      #print('shape a2:',nshape,file=fn)
      if len(nshape) == 1:
         a2vec = True
         a2.shape = (nshape[0],1)
      ns[2,:] = np.shape(a2)
   if a3 is not None:
      nshape = np.shape(a3)
      #print('shape a3:',nshape,file=fn)
      if len(nshape) == 1:
         a3vec = True
         a3.shape = (nshape[0],1)
      ns[3,:] = np.shape(a3)
   n = np.amax(ns[:,0],axis=0)
   nc[1] = ns[0,1]
   for i in range(1,4):
      nc[i+1] = nc[i] + ns[i,1]
   a = np.empty(nc[4])
   #print('shapes:',ns,file=fn)
   #print('rows, cols:',n,nc,file=fn)
   print(title,file=fn)
   for i in range(n):
      a[nc[0]:nc[1]] = a0[i,:]
      if ns[1,0] > 0:
         a[nc[1]:nc[2]] = a1[i,:]
      if ns[2,0] > 0:
         a[nc[2]:nc[3]] = a2[i,:]
      if ns[3,0] > 0:
         a[nc[3]:nc[4]] = a3[i,:]
      np.savetxt(fn,a,fmtf,delimiter='',newline='\t')
      print('',file=fn)
   if a0vec:
      a0.shape = (ns[0,0],)
   if a1vec:
      a1.shape = (ns[1,0],)
   if a2vec:
      a2.shape = (ns[2,0],)
   if a3vec:
      a3.shape = (ns[3,0],)

#def main():
#   import numpy as np
#   #import arrayprint as ap
#   a = np.arange(10).reshape(5,2)
#   b = np.ones(5)
#   fopen('atest.dat')
#   print('shape',a.shape,b.shape)
#   arrayprint(" a,b",a,b)
#   print('shape',a.shape,b.shape)
#   fclose()
#
#main()
   

