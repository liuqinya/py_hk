#!/usr/bin/env python

import pickle,  sys
import recfun
import obspy

download_event_and_station_file, download_waveform_and_write_to_sac, write_sac_cut_file= False,  False,  True
# test event range
if download_event_and_station_file:
    stime=obspy.UTCDateTime('2013-03-01')
    etime=obspy.UTCDateTime.now()
    #stime=obspy.UTCDateTime('2003-06-16T10:00:00')
    #etime=obspy.UTCDateTime('2003-06-17')
    rec=recfun.recdata()
# max_dep=200, search_center=(38,28), min_search_radius=30., max_search_radius=90
    rec.get_events(stime=stime, etime=etime, min_mag=6.5, event_list_txt='events_test')
    rec.event_catalog_plot()
# min_lat=20, max_lat=50, min_lon=20, max_lon=40
    rec.get_stations(network='TD',station='TD008',station_list_txt='stations_test')
    rec.station_inventory_plot()
# pickle rec object
    with open('rec.pkl', 'wb') as fp:
        pickle.dump(rec, fp, -1)    
if download_waveform_and_write_to_sac:
# read back pickle object
    with open('rec.pkl','rb') as fp:
        rec=pickle.load(fp)
    rec.get_waveforms(phase_select='P',pre_length=600,length=3000, remove_response=False, check_resp_same=True, filter_seismogram=True, f_low=0.05,  f_high=2., rotate=True, rotate_to_rt_in_sac=False, snr=True)
    with open('rec2.pkl', 'wb') as fp:
        pickle.dump(rec, fp, -1) 
if write_sac_cut_file:
    with open('rec2.pkl','rb') as fp:
        rec=pickle.load(fp)
    rec.cut_segment_from_phase(snr=5.0)
   

