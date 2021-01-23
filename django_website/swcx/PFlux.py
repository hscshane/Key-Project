#!/usr/local/bin/python

from __future__ import print_function
import numpy as np
#import math
from math import pi, sin,cos,atan2,log,exp
from pyhdf.V import *
from pyhdf.VS import *
from pyhdf.HDF import *
from pyhdf.SD import *
#from datetime import date
import datetime
#import sys

external_dir = 'D:\\OneDrive - University of Miami\\work\\COSMOS\\swcx\\Sicong\\swcx\\'

def PFlux(year1,month1,day1):
  filename = external_dir + "swepam_data_1day.hdf"
  hdf = HDF(filename)

#Initialize the V interface on the HDF file.
  v = hdf.vgstart()
  vs = hdf.vstart()
#Scan all vgroups in the file
  ref = -1
  refnum = v.getid(ref)
  vg = v.attach(refnum)
#print "----------------"
#print "name:", vg._name, "class:",vg._class, "tag,ref:",
#print vg._tag, vg._refnum
# Show the number of members of each main object type.
#print "members: ", vg._nmembers,
#print "datasets:", vg.nrefs(HC.DFTAG_NDG),
#print "vdataset:  ", vg.nrefs(HC.DFTAG_VH),
#print "vgroups: ", vg.nrefs(HC.DFTAG_VG)

# Read the contents of the vgroup.
  members = vg.tagrefs()
# Display info about each member.
  index = -1
  for tag, ref in members:
    index += 1
#    print "member index", index
  # Vdata tag
    if tag == HC.DFTAG_VH:
      vd = vs.attach(ref)
      nrecs, intmode, fields, size, name = vd.inquire()
#      print "  vdata:",name, "tag,ref:",tag, ref
#      print "    fields:",fields
#      print "    nrecs:",nrecs
      vd.detach()
  # SDS tag
    elif tag == HC.DFTAG_NDG:
      sds = sd.select(sd.reftoindex(ref))
      name, rank, dims, type, nattrs = sds.info()
#      print "  dataset:",name, "tag,ref:", tag, ref
#      print "    dims:",dims
#      print "    type:",type
      sds.endaccess()
 
  # VS tag
    elif tag == HC.DFTAG_VG:
      vg0 = v.attach(ref)
#      print "  vgroup:", vg0._name, "tag,ref:", tag, ref
      vg0.detach()
 
  # Unhandled tag
    else:
      print("unhandled tag,ref",tag,ref)

  # Close vgroup

  members = vg.tagrefs()
  (tag, ref) in members
  vd = vs.attach(ref)
  nrecs, intmode, fields, size, name = vd.inquire()
  alldata = vd.read(nrecs)
  vd.detach()
  vg.detach()
  v.end()
  hdf.close()

  data = np.array(alldata)

# input 
#  (y,m,d) = (1998,3,1)
  year = data[:,0]
  day = data[:,1]
  hr = data[:,2]
  min = data[:,3]
  sec  = data[:,4]
  fp_year = data[:,5]
  fp_doy = data[:,6]
  pdensity = data[:,8]
  speed  = data[:,11]

  start = datetime.date(int(year[0]),1,1) + datetime.timedelta(int(day[0])-1)
  end = datetime.date(int(year[-1]),1,1) + datetime.timedelta(int(day[-1])-1)
#  print "star Date and End Date:", start, end

  delta1 = datetime.date(year1,month1,day1) - start
  index = delta1.days
#  print "index ....", index
  if index >= (nrecs-1):
    print("ERROR: the input time is too new")
#  break
  elif index < 0:
    print("ERROR: the input time is too old")
#  break
  else:
    avePF = 0
    num = 0
    for i in range(0,90):
      j = index - i
      if j<0:
        print("ERROR: the index is out of the array")
        break
      else:
        if pdensity[j] > 0 and speed[j] >0:
          avePF = avePF + pdensity[j]*speed[j]
          num = num + 1
    avePF = avePF/num
#    print "the 90 days average proton flux is:",avePF, num
  return avePF 
