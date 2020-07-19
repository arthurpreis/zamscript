#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#  file_filter.py
#
#  Copyright 2019 Arthur Reis <arthurreis@INPE>
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with this program; if not, write to the Free Software
#  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
#  MA 02110-1301, USA.
#
#


import os

def check(t):
    '''returns true to files that do contain a solution'''
    ti = open(t,'r')
    for line in ti:
        if "FINAL MODEL" in line:
            return True 
    return False

#check was integrated in autom.
#deprecated:
files = []

for f in os.listdir("./"):
    if f.endswith(".dat") and int(f.split('.')[0])>i:
        files.append(os.path.join(f))

good_files = []
for f in files:
    print(f)
    if not check(f):
        os.remove(f)