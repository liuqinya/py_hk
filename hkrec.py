#!/usr/bin/env python3
"""
 Filename: hkrec.py
 Purpose: receiver function class that can for a given station/net/chan_code
  1) Download event data for events over a set of criteria
  2) Rotate and save the data in SAC format for iter_conv of the hk package
  3) Run the hk package to generate the receiver functions
  4) Visually select all receiver functions for selection
  5) stack all selected receiver functions in the hk package
  Author: Qinya Liu
  TODO: add SNR calculation to guide screening or offer a pre-screening?
"""

import numpy as np
import os,sys,glob,subprocess,copy, pickle

from obspy.clients.fdsn import Client

from obspy.core.stream import Stream
from obspy.core import AttribDict

from obspy import UTCDateTime
from obspy.geodetics.base import locations2degrees
from obspy.geodetics.base import gps2dist_azimuth
from obspy.taup import TauPyModel
from obspy.signal.rotate import rotate2zne

from mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as plt

import sac_utils

# =================
class event_search_par:
    """
    set the event search criteria for receiver function data download
    start/end time, min/max magnitude, min/max epicentral distance and min/max depth
    """
    def __init__(self,stime=UTCDateTime('2002-12-22'),etime=UTCDateTime('2003-12-31'),min_mag=5.5,max_mag=9.0, search_center=(38.54,28.63),
                 min_search_radius=30,max_search_radius=95, min_dep=0.,max_dep=800.):
        mag_type='Mw'
        self.magnitude_range= AttribDict({'mag_type': mag_type,'min_mag':min_mag,'max_mag':max_mag})
        self.time_range= AttribDict({'stime':stime,'etime':etime})
        self.distance_range= AttribDict({'center_lat': search_center[0],'center_lon':search_center[1],
                                         'min_radius_in_deg':min_search_radius, 'max_radius_in_deg':max_search_radius})
        self.depth_range= AttribDict({'min_depth':min_dep,  'max_depth':max_dep})

    def __str__(self):
        # write event search parameters
        output='Events search criteria: \n    with '+self.magnitude_range['mag_type'] +' between ['+str(self.magnitude_range['min_mag'])+', ' \
        +str(self.magnitude_range['max_mag'])+']\n    time between ['+str(self.time_range['stime'])[0:23]+', '+str(self.time_range['etime'])[0:23] \
        +']\n    distance within ['+str(self.distance_range['min_radius_in_deg'])+', '+str(self.distance_range['max_radius_in_deg']) \
        +'] deg from center point ['+ str(self.distance_range['center_lat'])+', '+str(self.distance_range['center_lon'])+'] deg\n    depth between [' \
        +str(self.depth_range['min_depth'])+', '+str(self.depth_range['max_depth'])+'] km\n' 
        return(output)

# ================================================
class waveform_search_par:
    """
        set the waveform search criteria for receiver function data download
        length before/after origin for the initial download
        length before/after P for the processing and 1D model used for phase picking
    """
    def __init__(self,phase_selected='P',length_before_phase=50,length_after_phase=150,length_before_origin=100, length_after_origin=2000, model1d='ak135'):
        self.phase_window= AttribDict({'phase_selected': phase_selected,  'length_before_phase': length_before_phase,  'length_after_phase': length_after_phase})
        self.data_window=AttribDict({'length_before_origin':length_before_origin, 'length_after_origin': length_after_origin})
        self.ref_model1d= model1d
        if self.phase_window['length_before_phase'] < 0:
            sys.exit('set window length_before_phase > 0')
        if phase_selected != 'P' and phase_selected != 'S':
            sys.exit('This is not a phase for receiver functions: '+phase_selected)
        
    def __str__(self):
        output='\nWaveform search criteria: \n' \
        +'    for phases '+self.phase_window['phase_selected'][0]+' between [-'+str(self.phase_window['length_before_phase'])+', '+\
               str(self.phase_window['length_after_phase'])+']'
        return(output)


# ===============================================================
class recfun:
    def __init__(self,sta='KUL',net='XH',chan='BH',loc_code='',data_starttime='2003-06-01',data_endtime='2003-06-30',data_client='IRIS'):
        self.sta=sta
        self.net=net
        self.chan=chan
        self.loc_code=loc_code
        self.event_par=event_search_par(stime=UTCDateTime(data_starttime),etime=UTCDateTime(data_endtime))
        self.waveform_par=waveform_search_par()
        self.client_name=data_client
        print('Request event catalog and waveform data from '+self.client_name+' ...')
        self.client=Client(self.client_name, debug=False)
        self.accepted_event_list=[]
        return
    
# ============================================
    def get_station_info(self,verbose=False):
        """
        todo : add location code
        """
        print('Getting station info ....')
        inv=self.client.get_stations(network=self.net,station=self.sta,
            starttime=self.event_par.time_range.stime, endtime=self.event_par.time_range.etime,
            level='response', channel=self.chan+'?')
        if len(inv) == 0:
            sys.exit('No such station'+self.sta+'.'+self.net)
        else:
            self.inv=inv
            if verbose:
                print('Inventory: \n', inv)
        return
    
# =============================================
    def get_events(self,plot_catalog=False,verbose=False):
         
        print('Search events based on criteria ...')
        self.event_par.distance_range['center_lat']=self.inv[0][0][0].latitude
        self.event_par.distance_range['center_lon']=self.inv[0][0][0].longitude
        if verbose:
            print(self.event_par)
        
        self.event_catalog = self.client.get_events( 
            starttime=self.event_par.time_range['stime'], 
            endtime=self.event_par.time_range['etime'],
            minmagnitude=self.event_par.magnitude_range['min_mag'],
            magnitudetype=self.event_par.magnitude_range['mag_type'],
            latitude=self.event_par.distance_range['center_lat'], 
            longitude=self.event_par.distance_range['center_lon'],
            minradius=self.event_par.distance_range['min_radius_in_deg'], 
            maxradius=self.event_par.distance_range['max_radius_in_deg'],
            mindepth=self.event_par.depth_range['min_depth'], 
            maxdepth=self.event_par.depth_range['max_depth'],
            orderby='time-asc')

        if len(self.event_catalog) == 0:
            sys.exit('    No event within search range. Check your search criteria')
        else:
            print('    Total number of events requested: '+str(len(self.event_catalog)))
            if plot_catalog:
                self.event_catalog.plot()
        return

#===========================================================

    def get_and_process_event_data(self,rotate=True,plot_seismogram=True,write_to_sac=True,verbose=False):

        before_origin_time=self.waveform_par.data_window['length_before_origin']
        after_origin_time=self.waveform_par.data_window['length_after_origin']
        phase=self.waveform_par.phase_window['phase_selected']

        # loop over events to get data
        self.recdata = [Stream() for i in range(len(self.event_catalog))]
        for (iev,ev) in enumerate(self.event_catalog):
            if verbose:
                print('Acquire data for event #', iev + '...')
            ev_time=ev.origins[0].time
            bulk=[]
            # set up bulk search command
            for net in self.inv:
                for sta in net:
                    if  sta.is_active(ev_time): # check channel name later?
                        bulk.append((net.code,sta.code,'*',self.chan+'?',ev_time-before_origin_time,ev_time+after_origin_time))
                #print('with search command:', bulk)

            try:
                stream=self.client.get_waveforms_bulk(bulk, attach_response=True)
            except BaseException:
                print('Skipping event, no data returned')
                continue

            stream.merge(method=1) # maybe add a better way of check the quality of merged data?
            if len(stream) != 3:
                sys.exit('More than 3 traces acquired for for '+str(ev))
    
            if verbose:
                print('Adding trace stats and write to sac ...') # turn write_sac True to see the original waveform
            sac_utils.stream_add_stats(stream,self.inv,ev,iev,write_sac=False,rotate_in_obspy=False)

            print('Plotting waveform for event '+str(ev_time)+', mag='+str(ev.magnitudes[0].mag)+' ...')
            # check if all response files are the same, should be true for modern instruments
            # TODO: add SNR calculation to guide screening or offer a pre-screening?
            sac_utils.process_stream_for_hk(
                stream,event_id=iev,write_to_sac=True,plot_seismogram=True,
                phase='P',time_before_phase=self.waveform_par.phase_window.length_before_phase,
                time_after_phase=self.waveform_par.phase_window.length_after_phase
            )
            self.recdata[iev]=stream
        self.accepted_event_list=sac_utils.accepted_event_list
# ================================================
# allow interactive checking

#         print('Processing data (remove instrument response, apply anti-aliasing filter) ...')
# # water-level of 100 is clearly too big, the default 60 is probably alright.
# # use plot=True to see the effect of water-level. water level 10 is also ok.
# stream.remove_response(water_level=50,output='VEL')
# if filter:
#     stream.detrend("linear")
#     stream.taper(max_percentage=0.05, type="hann")
#     stream.filter('bandpass',freqmin=f_min,freqmax=f_max,corners=2,zerophase=True)
# print('Write data to sac files ...')
        return

    def iter_decon(self):
        return

    def k_stack(self):
        return

    def ccp_stack(self):
        return


# =================================================================
#LQY HERE: seperate get_waveforms from process waveforms, 
    # def get_event_data(self,phase_select='P',pre_length=30,length=1200,model1d='ak135', \
    #     remove_response=False,filter_seismogram=False, f_low=0.01, f_high=10., \
    #     rotate=False, rotate_to_rt_in_sac=False, check_resp_same=False,snr=False, \
    #     pkl=False, pickle_name='rec.pkl'):
    
    #     # set up waveform parameters
    #     wf_par=waveform_search_par(phase_select=phase_select,pre_length=pre_length,length=length,model1d=model1d)
    #     ext='.d' if remove_response else ''
    #     signal_to_noise_ratio=0.
    #     wf_stream=Stream()
        
#stream.sort(keys=['distance','network','station','location','channel'])

# ===============================================


        
        
