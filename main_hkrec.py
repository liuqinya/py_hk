#!/usr/bin/env python3

import pickle,  sys
import obspy

from hkrec import recfun
# to find the station active time range: https://ds.iris.edu/SeismiQuery/station.htm (select time range)
# first test the example from HK package example directory KUL.BH?.01   2003/6/16 (167) 22:08:02.140
staname='KUL'
netname='XH'
#rec=recfun(sta='KUL',net='XH',chan='BH',loc_code='',data_starttime='2002-12-22',data_endtime='2003-12-31')
#rec=recfun(sta='KUL',net='XH',chan='BH',loc_code='',data_starttime='2003-06-15',data_endtime='2003-07-16')
#rec=recfun(sta='WINE',net='7A',chan='BH',loc_code='',data_starttime='2013-10-01',data_endtime='2016-10-01')
rec=recfun(sta=staname,net=netname,chan='BH',loc_code='',data_starttime='2003-06-15',data_endtime='2003-07-16')
rec.event_par.magnitude_range.min_mag=6.5
rec.get_station_info()
rec.get_events()
rec.get_and_process_event_data()
# more on data processing: http://eqseis.geosc.psu.edu/cammon/HTML/RftnDocs/prep01.html
# TODO: more to be added here
print('selected_event ids:', rec.accepted_event_list)
pickle.dump(rec,open('rec.pkl','wb'))

# output bash script to run hk
with open(staname+'.'+netname+'.sh', 'w') as fpp:
    fpp.write('#!/bin/bash\n')
    fpp.write('f1=3\n')
    fpp.write('f2=5\n')
    for ievn in rec.accepted_event_list:
        fpp.write('./hk/iter_decon -C-2/-10/80 -F1/$f1/-$f2 -N100 -T0.1 %s %s\n'%('./data/'+staname+'.'+netname+\
                '..BHZ.'+str(ievn)+'.p.sac','./data/'+staname+'.'+netname+'..BHR.'+str(ievn)+'.p.sac'))
    fpp.write('\n')
    fpp.write('./hk/k_stack -R30/60/1.6/2.0 -I0.5/0.01 -G%s %s\n'%('./'+staname+'.'+netname+\
            '.grd','./data/'+staname+'.'+netname+'..BHR.*.p.saci'))
    fpp.write('./hk/grdmin -D %s\n'%('./'+staname+'.'+netname+'.grd'))
    fpp.write('./hk/grdmin -D %s | ./hk/hk_plt.pl > %s'%('./'+staname+'.'+netname+'.grd','./'+staname+'.'+netname+'.ps'))
