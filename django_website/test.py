# -*- coding: utf-8 -*-
"""
Created on Fri Jun  5 00:36:48 2020

@author: hscsh
"""

import sys
import datetime

time=datetime.datetime.now()

output=f"Hi {sys.argv[1]} current time is {time}"

print(output)