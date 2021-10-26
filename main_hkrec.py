#!/usr/bin/env python3

import pickle,  sys
import obspy

from hkrec import recfun
# to find the station active time range: https://ds.iris.edu/SeismiQuery/station.htm (select time range)
# first test the example from HK package example directory KUL.BH?.01   2003/6/16 (167) 22:08:02.140
#rec=recfun(sta='KUL',net='XH',chan='BH',loc_code='',data_starttime='2002-12-22',data_endtime='2003-12-31')
rec=recfun(sta='KUL',net='XH',chan='BH',loc_code='',data_starttime='2003-06-15',data_endtime='2003-07-16')
#rec=recfun(sta='WINE',net='7A',chan='BH',loc_code='',data_starttime='2013-10-01',data_endtime='2016-10-01')
rec.event_par.magnitude_range.min_mag=6.5
rec.get_station_info()
rec.get_events()
rec.get_and_process_event_data()
# more on data processing: http://eqseis.geosc.psu.edu/cammon/HTML/RftnDocs/prep01.html
# TODO: more to be added here
print('selected_event ids:', rec.accepted_event_list)
pickle.dump(rec,open('rec.pkl','wb'))