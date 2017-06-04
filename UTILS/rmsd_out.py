# -*- coding: utf-8 -*-
# @Author: zyc
# @Date:   2017-05-27 10:25:53
# @Last Modified by:   zyc
# @Last Modified time: 2017-05-27 19:38:33
import base
import readers
try:
   import numpy as np
except:
   import mynumpy as np
import os.path
import sys

if len(sys.argv) < 3:
   base.Logger.log("Usage is %s trajectory topology [output]" %sys.argv[0],base.Logger.CRITICAL)
   sys.exit()

l=readers.LorenzoReader (sys.argv[1],sys.argv[2])
stmp=l.get_system()
nconf=0


msd=np.zeros(nconf, dtype=np.float64)
nmsd=np.zeros(nconf, dtype=np.int64)

while stmp:
   s=[]
   nconf+=1
   print >> sys.stderr, " --- Working on conf %d" % (nconf)

   for i in range(len(s)):
      for j in range(i):
         strands1=s[i]._strands
         strands2=s[j]._strands

         dt= i - j
         for k in range (len(strands1)):
            r1 = strands1[k].cm_pos
            r2 = strands2[k].cm_pos
            dr=r1-r2
            msd[dt]+=np.dot(dr, dr)
            nmsd[dt] += 1

   s.append(stmp)
   stmp = l.get_system()

out = open ("rmsd.dat", "w")
for i in range(1, 200):
   print >> out, "%7i %10.5g" % (i, msd[i])
   #print >> out, "%7i %10.5g" % (i, msd[i] / nmsd[i])
out.close()

print >> sys.stderr, "All DONE"
