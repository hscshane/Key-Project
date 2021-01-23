#!/usr/local/bin/python

from __future__ import print_function
import numpy as np
#import math
from math import pi, sin,cos,atan2,log,exp

external_dir = 'D:\\OneDrive - University of Miami\\work\\COSMOS\\swcx\\Sicong\\swcx\\'

def ReadDensGrid(year):
  global HR, FH, HER, FHE, HLONG, HLAT, HeLONG, HeLAT

  NHlong = 361
  NHlat = 181
  NHR = 43
  NHelong = 73
  NHelat = 37
  NHeR = 59
  HLONG = np.empty((NHlong))
  HLAT = np.empty((NHlat))
  HeLONG = np.empty((NHelong))
  HeLAT = np.empty((NHelat))
  HR = np.empty((NHlong, NHlat, NHR),dtype=np.float32)
  FH = np.empty((NHlong, NHlat, NHR),dtype=np.float32)
  HER = np.empty((NHelong, NHelat, NHeR),dtype=np.float32)
  FHE = np.empty((NHelong, NHelat, NHeR),dtype=np.float32)
  nyear = str(year)
  INH = external_dir + "NeutralDensity\\"+'HDens'+nyear+'.txt'
  INHE = external_dir + "NeutralDensity\\"+'HeDens'+nyear+'eit.txt'
#  print(INH)
#  INH = "HDens2000.txt"
#  INHE = "HeDens2000eit.txt"
#  print(INH,INHE)
  fh = open(INH, 'r')
  fhe = open(INHE, 'r')
  hlines = fh.readlines()
  helines = fhe.readlines()
  fh.close()
  fhe.close()
#  print(helines[30],hlines[27])
  Hstart = 21
  Hestart = 21
  n=Hstart
  for i in range(0,NHlong):
    for j in range(0,NHlat):
      for k in range(0,NHR):
        string = hlines[n].split()
        HLONG[i] = string[0]
        HLAT[j] = string[1]
        HR[i,j,k] = string[2]
        FH[i,j,k] = float(string[3]) * 0.1
#        print(string[0], string[1], string[2], string[3])
        n=n+1
#  print "n=",n
  n=Hestart
  for i in range(0,NHelong):
    for j in range(0,NHelat):
      for k in range(0,NHeR):
        string = helines[n].split()
        HeLONG[i] = string[0]
        HeLAT[j] = string[1]
        HER[i,j,k] = string[3]
        FHE[i,j,k] = float(string[4]) * 0.015
#        print(string[0], string[1], string[2], string[3])
        n=n+1
#  print "n=",n

def Emissivity():
  global HR, FH, HER, FHE, HLONG, HLAT, HeLONG, HeLAT
  global OBS, DLOS, SX, SY, SZ, SD
  global RP, DSS, EH, EHE, SELON, SELAT, SHLON, SHLAT

  RP = np.empty((200))
  DSS = np.empty((200))
  EH = np.empty((200))
  EHE = np.empty((200))
  SELON = np.copy(RP)
  SELAT = np.copy(RP)
  SHLON = np.copy(RP)
  SHLAT = np.copy(RP) 

  DOBS = np.empty((3),dtype=np.float32)
  SMAX = 95.0
  I = 0
  S = 0.
  DS = 0.01
  while (S + DS * 0.5) < SMAX:
    S = S + DS * 0.5
  
    for i in range(0,3):
      DOBS[i] = OBS[i] + DLOS[i] * S
    RP[I] = (DOBS[0]**2.+ DOBS[1]**2.+ DOBS[2]**2.)**0.5
    XD = RP[I]
    LON = atan2(DOBS[1], DOBS[0]) * 180/pi
    if LON < 0.:
      LON = LON + 360
    LAT = atan2(DOBS[2], (DOBS[0]**2.+ DOBS[1]**2.)**0.5) * 180/pi

# Helio-Ecliptic COORDINATES of LOS STEP
    SELON[I] = LON
    SELAT[I] = LAT

# SDEFINITION of Helio-Ecliptic RADIAL VECTOR of LOS STEP
    SNLOS = np.empty((3),dtype=np.float32)
    SNLOS[0] = cos(SELON[I]*pi/180) * cos(SELAT[I]*pi/180)
    SNLOS[1] = sin(SELON[I]*pi/180) * cos(SELAT[I]*pi/180)
    SNLOS[2] = sin(SELAT[I]*pi/180)

# TRANSFORMATION to Heliographic COORDINATES SYSTEM
# ROTATION AROUND THE SUN REFERENCE AXES
    SD = np.empty((3),dtype=np.float32)
    SD[0] = SNLOS[0] * SX[0] + SNLOS[1] * SX[1] + SNLOS[2] * SX[2]
    SD[1] = SNLOS[0] * SY[0] + SNLOS[1] * SY[1] + SNLOS[2] * SY[2]
    SD[2] = SNLOS[0] * SZ[0] + SNLOS[1] * SZ[1] + SNLOS[2] * SZ[2]

# SHLON, SHLAT: HELIOGRAPHIC COORDINATES OF LOS STEP
    SHLON[I] = atan2(SD[1],SD[0]) * 180/pi
    SHLAT[I] = atan2(SD[2], (SD[0]**2. + SD[1]**2.)**0.5) * 180/pi

    ILON = int(LON) + 1
    ILAT = int(LAT+90.) + 1
#    print XD, LON, LAT
    NC_INT_HE = InterPolateHe(XD, LON, LAT)
    EHE[I] = NC_INT_HE

    if RP[I] < 0.35 and DS < 0.04:
      EH[I] = 0.
    else:
      NC_INT_H = InterPolateH(XD, LON, LAT)
      EH[I] = NC_INT_H

    DSS[I] = DS
    S = S + DS * 0.5
    DS = 0.1 + (RP[I] - 1.) * 0.1
    I = I + 1
  
  NRP = I - 1
  return NRP


def InterPolateH(D,LON, LAT):
  global HLONG, HLAT, FH, HR
  NHLONG = 361; NHLAT = 181; NHR = 43
  NI = 0; NJ = 0; NK = 0
  for i in range(1,NHLONG):
    if HLONG[i] > LON:
      NI = i
      break
  for j in range(1, NHLAT):
    if HLAT[j] > LAT:
      NJ = j
      break
  if D < HR[NI,NJ,0]:
    NC_INT_H = 0.0
  else:
    for k in range(1,NHR):
      if HR[NI,NJ,k] > D:
        NK = k
        break
      else:
        NK = NHR - 1

  DLON = (LON - HLONG[NI-1])/(HLONG[NI] - HLONG[NI-1])
  DLAT = (LAT - HLAT[NJ-1])/(HLAT[NJ] - HLAT[NJ-1])
  PR = (D - HR[NI,NJ,NK-1])/(HR[NI,NJ,NK] - HR[NI,NJ,NK-1])

  NC_INT_H = FH[NI-1,NJ-1,NK-1]*(1.0 - DLON)*(1.0 - DLAT)*(1.0 - PR) \
           + FH[NI,NJ-1,NK-1] * DLON * (1.0 - DLAT) * (1.0 - PR) \
           + FH[NI-1,NJ,NK-1] * (1.0 - DLON) * DLAT * (1.0 - PR) \
           + FH[NI-1,NJ-1,NK] * (1.0 - DLON) * (1.0 - DLAT) * PR \
           + FH[NI,NJ,NK-1] * DLON * DLAT * (1.0 - PR) \
           + FH[NI,NJ-1,NK] * DLON * (1.0 - DLAT) * PR \
           + FH[NI-1,NJ,NK] * (1.0 - DLON) * DLAT * PR \
           + FH[NI,NJ,NK] * DLON * DLAT * PR
  return NC_INT_H


def InterPolateHe(D,LON, LAT):
  global FHE, HER, HeLONG, HeLAT 
#  print "***",D,LON,LAT
  NHeLONG = 73; NHeLAT = 37; NHeR = 59
  NI = 0.
  NJ = 0.
  NK = 0.
  for i in range(1,NHeLONG):
    if HeLONG[i] > LON:
      NI = i
      break
  for j in range(1, NHeLAT):
    if HeLAT[j] > LAT:
      NJ = j
      break
  if D < 1.5:
    for k in range(1,NHeR):
      if HER[NI,NJ,k] > D:
        NK = k
        break

    DLON = (LON - HeLONG[NI-1])/(HeLONG[NI] - HeLONG[NI-1])
    DLAT = (LAT - HeLAT[NJ-1])/(HeLAT[NJ] - HeLAT[NJ-1])
    PR = HER[NI,NJ,NK] - HER[NI,NJ,NK-1]
    X = D /HER[NI,NJ,NK-1]
    PA = log(HER[NI,NJ,NK]/HER[NI,NJ,NK-1])

    FHE1 = FHE[NI-1,NJ-1,NK-1] * (1.0 - DLON) * (1.0 - DLAT) \
         + FHE[NI-1,NJ,NK-1] * (1.0 - DLON) * DLAT \
         + FHE[NI,NJ-1,NK-1] * DLON * (1.0 - DLAT) \
         + FHE[NI,NJ,NK-1] * DLON * DLAT
    FHE2 = FHE[NI-1,NJ-1,NK] * (1.0 - DLON) * (1.0 - DLAT) \
         + FHE[NI-1,NJ,NK] * (1.0 - DLON) * DLAT \
         + FHE[NI,NJ-1,NK] * DLON * (1.0 - DLAT) \
         + FHE[NI,NJ,NK] * DLON * DLAT
    if FHE1 != 0.0:
      NC_INT_HE = FHE1 * exp(log(FHE2/FHE1) * log(X) / PA)
    else:
      NC_INT_HE = 0.0
  else:
    for k in range(1,NHeR):
      if HER[NI,NJ,k] > D:
        NK = k
#        print "k=",k
        break
      else:
        NK = NHeR - 1
#    print(NI,NJ,NK,D,HER[NI,NJ,NK],k)
    DLON = (LON - HeLONG[NI-1])/(HeLONG[NI] - HeLONG[NI-1])
    DLAT = (LAT - HeLAT[NJ-1])/(HeLAT[NJ] - HeLAT[NJ-1])
#    print NI,NJ,NK,D,NHeR
    PR = (D - HER[NI,NJ,NK-1])/(HER[NI,NJ,NK] - HER[NI,NJ,NK-1])

    NC_INT_HE = FHE[NI-1,NJ-1,NK-1] * (1.0 - DLON) * (1.0 - DLAT) * (1.0 - PR) \
              + FHE[NI,NJ-1,NK-1] * DLON * (1.0 - DLAT) * (1.0 - PR) \
              + FHE[NI-1,NJ,NK-1] * (1.0 - DLON) * DLAT * (1.0 - PR) \
              + FHE[NI-1,NJ-1,NK] * (1.0 - DLON) * (1.0 - DLAT) * PR \
              + FHE[NI,NJ,NK-1] * DLON * DLAT * (1.0 - PR) \
              + FHE[NI,NJ-1,NK] * DLON * (1.0 - DLAT) * PR \
              + FHE[NI-1,NJ,NK] * (1.0 - DLON) * DLAT * PR \
              + FHE[NI,NJ,NK] * DLON * DLAT * PR
#  print PR,NC_INT_HE
  return NC_INT_HE


def NeuInt(EarthLong, elong, elat,year):
  global OBS, DLOS, SX, SY, SZ, SD
  UA = 1.495985e15 # AU in cm

#--- DEFINITION OF HELIOGRAPHIC REFERENCE SYSTEM SX, SY, SZ
#--- OMEGA: LONGITUDE OF THE SOLAR EQUATOR ASCENDING NODE 
#--- BETA: INCLINATION OF SOLAR AXE WITH RESPECT TO THE ECLIPTIC AXE
#--- NOT NECESSARY FOR DENSITY INTEGRALS THROUGH INTERPOLATION
  OMEGA = 73.67 
  BETA = 7.25
  SX = np.empty((3),dtype=np.float32)
  SY = np.empty((3),dtype=np.float32)
  SZ = np.empty((3),dtype=np.float32)

  SX[0] = cos(OMEGA * pi/180)
  SX[1] = sin(OMEGA * pi/180)
  SX[2] = 0.0

  SY[0] = -sin(OMEGA * pi/180) * cos(BETA * pi/180)
  SY[1] = cos(OMEGA * pi/180) * cos(BETA * pi/180)
  SY[2] = sin(BETA * pi/180)

  SZ[0] = sin(OMEGA * pi/180) * sin(BETA * pi/180)
  SZ[1] = -cos(OMEGA * pi/180) * sin(BETA * pi/180)
  SZ[2] = cos(BETA * pi/180)

#  print("SX SY SZ\n",SX,SY,SZ)
  ReadDensGrid(year)

#--- DEFINITION HELIUM WIND AXIS LAMDA = 74.7, BETA = -5.3
#--- A varier is on veut faire des tests sur les parametres IBEX
  HE = np.empty((3),dtype=np.float32)
  HELON = 74.7
  HEDEC = -5.3
  HE[0] = cos(HELON * pi /180.0) * cos(HEDEC * pi /180.0)
  HE[1] = sin(HELON * pi /180.0) * cos(HEDEC * pi /180.0)
  HE[2] = sin(HEDEC * pi /180.0)
#  print("HELIUM AXIS",HE)

# read LOS file
#  fp = open('LOS.txt', 'r')
#  fout = open('out.txt', 'w')
#  lines = fp.readlines()
#  print(lines,lines[0])
#  fp.close()
#  for i in range(0,1):
#    list = [float(f) for f in lines[i].split()]
  LOBS = EarthLong
  ELON = elong
  ELAT = elat
#  print("LOS.txt",LOBS, ELON, ELAT)

  ROBS = 1.0
  BOBS = 0.0

#---DEFINITION OF EARTH POSITION VECTOR
  OBS = np.empty((3),dtype=np.float32)
  OBS[0] = ROBS * cos(LOBS * pi/180.0) * cos(BOBS * pi/180.0)
  OBS[1] = ROBS * sin(LOBS * pi/180.0) * cos(BOBS * pi/180.0)
  OBS[2] = ROBS * sin(BOBS * pi/180.0)
#  print("OBS", OBS)

#--- DEFINITION OF LOS VECTOR
  DLOS = np.empty((3),dtype=np.float32)
  DLOS[0] = cos(ELON *pi/180.0) * cos(ELAT * pi/180.0)
  DLOS[1] = sin(ELON *pi/180.0) * cos(ELAT * pi/180.0)
  DLOS[2] = sin(ELAT *pi/180.0)
#  print("DLOS", DLOS)

  NRP = Emissivity()
#  print(NRP)
  EMIH = 0.0
  EMIHE = 0.0
  for j in range(0, NRP+1):
    EMIH = EMIH + EH[j] * DSS[j]/RP[j]**2.0
    EMIHE = EMIHE + EHE[j] * DSS[j]/RP[j]**2.0

#  print(LOBS, ELON, ELAT, EMIH, EMIHE)
  return (EMIH, EMIHE)
