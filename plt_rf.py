#!/usr/bin/env python3

import os
import obspy
import numpy as np
from matplotlib import pyplot as plt

def get_file_list(basis_dir="./", begin="", end=""):
  path_list = os.listdir(basis_dir)
  list_final = []
  for partial in path_list:
    if begin and end:
      if partial[:len(begin)] == begin and partial[-len(end):] == end:
        list_final.append(partial)
    elif end:
      if partial[-len(end):] == end:
        list_final.append(partial)
    elif begin:
      if partial[:len(begin)] == begin:
        list_final.append(partial)
    else:
      list_final.append(partial)
  return list_final

staname='LBDL'
netname='7A'
scale_num=0.005
file_list = get_file_list("./data",begin=staname+'.'+netname,end=".saci")
plt.figure()
for ifile in file_list:
  print('File: %s'%(ifile))
  st=obspy.read("./data/"+ifile)
  print(st[0].stats.sac.get('user0'))
  plt.plot(np.arange(st[0].stats.sac.get('b'),st[0].stats.sac.get('e'),\
          st[0].stats.sac.get('delta')),st[0].stats.sac.get('user0')+\
          scale_num*st[0].data/max(st[0].data),linewidth='1',color='black')
  plt.fill_between(np.arange(st[0].stats.sac.get('b'),st[0].stats.sac.get('e'),\
          st[0].stats.sac.get('delta')),st[0].stats.sac.get('user0'),\
          st[0].stats.sac.get('user0')+scale_num*st[0].data/max(st[0].data),\
          where=scale_num*st[0].data>0,facecolor='black')
plt.xlim(-3,20)
plt.xlabel('Time (s)')
plt.ylabel('p (s/km)')
plt.show()
