# -*- coding: utf-8 -*-
# @Author: zyc
# @Date:   2017-05-22 03:41:22
# @Last Modified by:   zyc
# @Last Modified time: 2017-05-22 13:14:16
# This script is to generate the pdb file at specific timestep

import base
import readers
try:
    import numpy as np
except:
    import mynumpy as np
import os.path
import sys

if len(sys.argv)<3:
    base.Logger.log("Usage is %s configuration topology [output]" %sys.argv[0], base.Logger.CRITICAL)
    sys.exit()

if len(sys.argv) > 3:
    output = sys.argv[3]
else: output =sys.argv[1]+".pdb"

l=readers.LorenzoReader (sys.argv[1], sys.argv[2])
s=l.get_system ()

append=False

while s:
    pass

