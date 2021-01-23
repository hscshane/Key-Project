from __future__ import print_function
import os
import sys
import cgi
import ephem
import datetime
import math
import time
import numpy as np
#from pyslalib import slalib
import astropy.units as u
from astropy.coordinates import BarycentricTrueEcliptic
from astropy.coordinates import FK5
from astropy.coordinates import Galactic
import NeuInt
import PFlux

def CheckDate(date):
  if "/" in date:
    month = date.split("/")[0]
    day = date.split("/")[1]
    year = date.split("/")[2]
  else:
    date = date.strip()
    month = date.split()[0]
    day = date.split()[1]
    year = date.split()[2]
  month = int(month)
  day = int(day)
  year = int(year)
  return (month, day, year)

def CheckCoord(Coord, CoordType):
  degtorad = math.pi/180
  Coord = Coord.strip()
  if "," in Coord:
    lon = Coord.split(",")[0].strip()
    lat = Coord.split(",")[1].strip()
    if len(lon.split()) == 3:
      lon_hr = float (lon.split()[0])
      lon_min = float (lon.split()[1])
      lon_sec = float (lon.split()[2])
      lat_deg = float (lat.split()[0])
      lat_min = float (lat.split()[1])
      lat_sec = float (lat.split()[2])
      lon = lon_hr*15. + lon_min/60. + lon_sec/3600.
      lat = lat_deg + lat_min/60. + lat_sec/3600.
    else:
      lon = float(lon)
      lat = float(lat)
  else:
    lon = Coord.split()[0]
    lat = Coord.split()[1]
    lon = float (lon)
    lat = float (lat)
  if CoordType == "J2000":
#mjd for 2000-01-01
    #(elon,elat) = slalib.sla_eqecl(lon*degtorad, lat*degtorad, 51544.0)
    gc = FK5(ra=lon*u.degree, dec=lat*u.degree)
    hec = gc.transform_to(BarycentricTrueEcliptic)
    (elon,elat) = (hec.lon.value*degtorad, hec.lat.value*degtorad)
  elif CoordType == "Galactic":
    #(ra,dec) = slalib.sla_galeq(lon*degtorad, lat*degtorad)
    #(elon,elat) = slalib.sla_eqecl(ra,dec,51544.0)
    gc = Galactic(l=lon*u.degree,b=lat*u.degree)
    hec = gc.transform_to(BarycentricTrueEcliptic)
    (elon,elat) = (hec.lon.value*degtorad, hec.lat.value*degtorad)
  else:
    elon = lon*degtorad
    elat = lat*degtorad
#    return (elon, elat)
  return (elon, elat)

def EarthEclip(date):
  sun = ephem.Sun(date)
#heliocentric ecliptic coordinate for Earth
  ElonEarth = sun.hlon
  ElatEarth = sun.hlat
  return (ElonEarth.real, ElatEarth.real)

# form = cgi.FieldStorage()
# newdate1 = form.getvalue("DateVal")
# coord = form.getvalue("CoordVal")
# coordtype = form.getvalue("CoordType")
# HeHRatio = form.getvalue("CSRatio")
coord = sys.argv[1]
coordtype = sys.argv[2]
newdate1 = sys.argv[3]
#HeHRatio = sys.argv[4]

(month, day, year) = CheckDate(newdate1)
date = datetime.date(year,month,day)
(elonEarth, elatEarth) = EarthEclip(date)

(losElon, losElat) = CheckCoord(coord, coordtype)
degtorad = math.pi/180
#sys.stdout = os.fdopen(sys.stdout.fileno(),'w',0)

avePF = PFlux.PFlux(year,month,day)
print("Date:",str(year)+"-"+str(month)+"-"+str(day),'|')
coord = coord.strip()
coord1 = coord.split(',')[0]
coord2 = coord.split(',')[1]
if len(coord1.split()) == 3:
  x = [coord1.split()[i] for i in range(3)]
  y = [coord2.split()[i] for i in range(3)]
  print("Coordinates: (",x[0],"H",x[1],"M",x[2],"S",",",y[0],"D",y[1],"M",y[2],"S",") | ",sep='',)
else:
  print("Coordinates: (",coord1,"degs,",coord2,"degs)"," | ",sep='')
print("System:",coordtype, '| ')

neutralint = NeuInt.NeuInt(elonEarth/degtorad, losElon/degtorad, losElat/degtorad,year)

norm = 8.89E-7
index = 2.254
HeH = 2.
#HeHRatio = float(HeHRatio)
O7 = norm*avePF**index * (neutralint[0]+HeH*neutralint[1])
#O7 = 1.0
print("He&H Ratio:",HeH,"|")

O7 = "{0:.4g}".format(O7)

print ("OVII =", O7, "LU") 