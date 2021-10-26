#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Originally written in March-April 2018
rewrite in October 2018 to add P, S picking and rotation
Some notes: 1)UTCDateTime takes only microsec, while the sac header asks nzmsec in millisec (see sac header:
#http://ds.iris.edu/files/sac-manual/manual/file_format.html

Q. Liu
"""

import numpy as np
import os,sys,glob,subprocess,copy, pickle
from obspy.clients.fdsn import Client

from obspy.core.stream import Stream
from obspy.core.util.attribdict import AttribDict

from obspy import UTCDateTime
from obspy.geodetics.base import locations2degrees
from obspy.geodetics.base import gps2dist_azimuth, degrees2kilometers
from obspy.taup import TauPyModel
from obspy.signal.rotate import rotate2zne,rotate_ne_rt

from mpl_toolkits.basemap import Basemap

import matplotlib.pyplot as plt 
from matplotlib.dates import date2num
from matplotlib.widgets import Button
from obspy.taup import TauPyModel, taup_geo

#from sets import Set

radiusOfEarth = 6371.0
accepted_event_list=[]

def stream_add_stats(data_stream,inv,evt,id_evt=0,write_sac=False,rotate_in_obspy=False,taup_model='ak135',plot_seismogram=True,verbose=False):
    """
    this function will add event, station, distance and other sac header
    information for all data stream for one particular event 
    add stats headers including trace.stats.distance/latitude/longitude
    figure out channel code number before select stream ...
    this works for all data stream belongs to one event! 

    Todo: 1) test if there are multiple location codes; right now only work with one stream. update in future to deal
             with multiple channel (total_number_of channels)
          2) check if rotation is fine for BH[12] or BH[EN] when the cmpaz is not 0/90. 
    """
    model = TauPyModel(model=taup_model)

    # loop over all station in the inventory and find the stream for sta.net 
    for net in inv:
        for sta in net:
            str1=data_stream.select(network=net.code,station=sta.code)
            if len(str1) == 0:
                continue
            
            if len(str1) != 3: # does not deal with multi-channel-code case yet
                sys.exit('Problem: missing components', str1)

            # initialize P and S arrival times to be set for str1
            Ptime=0; Stime=0

            for (j,tr) in enumerate(str1):
                # check the consistency of channel code 
                for chan in sta:
                    if tr.stats.channel == chan.code and tr.stats.location == chan.location_code:
                        break
                else:
                    sys.exit('Problem finding channel in inventory for trace ',tr)
                # add trace statistics from inventory and event catalog for sac writing later
                tr.stats.station_coordinates={'latitude':chan.latitude,'longitude':chan.longitude,'elevation':chan.elevation}
                (tr.stats.gcarc,tr.stats.azimuth,tr.stats.back_azimuth)=taup_geo.calc_dist_azi(evt.origins[0].latitude,
                    evt.origins[0].longitude,chan.latitude, chan.longitude,radius_of_planet_in_km=radiusOfEarth,flattening_of_planet=0)
                tr.stats.distance=degrees2kilometers(tr.stats.gcarc)*1000 # in meters
                tr.stats.event_origin=evt.origins[0] # depth in m
                tr.stats.event_mag=evt.magnitudes[0]
                tr.stats.dip=chan.dip
                tr.stats.cmpaz=chan.azimuth
                if verbose:
                    print('tr.stats=',tr.stats)
                if j == 0:
                    # set P and S arrival times and ray angles based on taupmodel(ak135)
                    arrivals = model.get_travel_times(source_depth_in_km=evt.origins[0].depth/1000.,distance_in_degree=tr.stats.gcarc)
                    for k in range(len(arrivals)):  # always take the first P or S arrival
                        if (arrivals[k].name == 'P' or arrivals[k].name == 'p') and Ptime == 0:
                             # ray parameter=rsin(th)/v in s/radians
                            (Ptime, Prayp, Pinc_angle)=(arrivals[k].time, arrivals[k].ray_param, arrivals[k].incident_angle)
                            #  print('P traveltime, ray_p, inc_angle=: ', arrivals[k].time,arrivals[k].ray_param,arrivals[k].incident_angle,evt.origins[0].depth/1000.)
                        elif (arrivals[k].name == 'S' or arrivals[k].name == 's') and Stime == 0:
                            (Stime, Srayp, Sinc_angle)=(arrivals[k].time, arrivals[k].ray_param,arrivals[k].incident_angle)

                    if verbose:
                        print('%s.%s Dist %4.1f   Ptime %4.1f   Stimes %4.1f' % (tr.stats.network, tr.stats.station, tr.stats.distance/1000., Ptime,Stime))
                    (Parr,Sarr)=(evt.origins[0].time+Ptime, evt.origins[0].time+Stime)
                tr.stats.Parr={'arrival_time':Ptime,'time':Parr,'rayp':Prayp,'inc_angle':Pinc_angle}
                tr.stats.Sarr={'arrival_time':Stime,'time':Sarr,'rayp':Srayp,'inc_angle':Sinc_angle}

            if write_sac:
                write_stream_to_sac(str1,ext=str(id_evt))
            # implementation rotation through obspy directly instead of using sac

def process_stream_for_hk(stream, event_id=0, rotate_in_obspy=True, write_to_sac=True,
                          phase='P', time_before_phase=50, time_after_phase=150,
                          plot_seismogram=True, screen_just_P=True,
                          filter=True, filter_fmin=0.05, filter_fmax=2.0,verbose=False):
    """
    This function checks the response files of all traces, detrend, rotate them and cut them and write them out
    for receiver functions.
    """

    if (stream[0].stats.response == stream[1].stats.response and stream[0].stats.response == stream[
        2].stats.response):
        # detrend but without removing instrument response
        stream.detrend('linear')
        # rotate from original Z E N or Z 1 2 to Z R T
        if rotate_in_obspy:
            str2 = rotate_trace_to_zrt(stream)
        else:
            str2=stream
        # cut slice
        if phase == 'P':
            time = str2[0].stats.Parr.time
        elif phase == 'S':
            time = str2[0].stats.Sarr.time
        else:
            sys.exit('Error phase name', phase)
        rf_starttime = time - time_before_phase
        rf_endtime = time + time_after_phase
        if verbose:
            print('slicing stream from ', rf_starttime, ' to ', rf_endtime)

        if filter:
            str2.detrend('linear')
            str2.filter('bandpass', freqmin=filter_fmin, freqmax=filter_fmax, zerophase=True)

        str3 = str2.slice(starttime=rf_starttime, endtime=rf_endtime)
        str3.detrend('linear')
        if plot_seismogram:
            if screen_just_P:
                str_plt=str3
            else:
                str_plt=str2
            # PLOT WITH SOME PRELIMINARY FILTERING FIRST [1, 30 sec] period band; need to be formalized
            plot_stream_with_picks(str_plt, title=str(event_id), add_accept_button=True)
        if write_to_sac:
            write_stream_to_sac(str3, ext=str(event_id) + '.p')
    else:
        # to be added probably in the future when instrument responses are different
        sys.exit('The instrument responses are different for the traces in stream', stream)
    # add data processing here(maybe move rotation here)
    # filter, instrument response cut

    # --------------------------
def rotate_trace_to_zrt(str1):

    """
      This function rotates a stream object (either in ENZ or 12Z) in to a stream object of TRZ
    """ 
  
    for k, tr in enumerate(str1):
        if tr.stats.channel[2] == 'Z':
            tr_z=tr
        elif tr.stats.channel[2] == 'E' or tr.stats.channel[2] == '1':
            tr_e=tr
        elif tr.stats.channel[2] == 'N' or tr.stats.channel[2] == '2':
            tr_n=tr
        else:
            sys.exit('Problem with '+tr.stats.channel)

    if 'tr_e' not in vars():
        sys.exit('No east component')
    else:
        tr_r=tr_e.copy()
    if 'tr_n' not in vars():
        sys.exit('No north component')
    else:
        tr_t=tr_n.copy()
    if 'tr_z' not in vars():
        sys.exit('No vertical component')
    else:
        tr_z2=tr_z.copy()

    # note that even traces labelled as E or N in the channel name may not have the right cmpaz.
    (tr_z2.data, tr_n.data, tr_e.data) = rotate2zne(tr_z.data, tr_z.stats.cmpaz, tr_z.stats.dip,
                                tr_n.data, tr_n.stats.cmpaz, tr_n.stats.dip,
                                tr_e.data, tr_e.stats.cmpaz, tr_e.stats.dip)
    baz=tr_z2.stats.back_azimuth
    #baz=19.29817039741116 vs baz_sac= 1.937192e+01 # the baz in sac is not accurate enough that could cause discrepancy with the rotation results
    # what is in rotate_ne_to_tr (confirmed)
    # ba = radians(ba)
    #r = - e * sin(ba) - n * cos(ba)
    #t = - e * cos(ba) + n * sin(ba)
    (tr_r.data,tr_t.data)=rotate_ne_rt(tr_n.data,tr_e.data,baz)

    # after rotation, inclination is not changed but azimuth is changed
    (tr_r.stats.cmpinc,tr_t.stats.cmpinc)=(90,90)
    (tr_r.stats.cmpaz,tr_t.stats.cmpaz)=((baz+180)%360., (baz+270)%360.) 
    (tr_r.stats.channel,tr_t.stats.channel)=(tr_r.stats.channel[0:2]+'R',tr_t.stats.channel[0:2]+'T')

    str2=Stream(traces=[tr_t,tr_r,tr_z2])
    return str2


def write_stream_to_sac(str1,write_dir='data',ext='',verbose=False):

    if ext != '':
        ext='.'+ext
    if not os.path.isdir(write_dir):
        sys.exit('No such dir to write sac',write_dir)

    for tr in str1:
        sac= AttribDict()
        (sac.kstnm,sac.knetwk,sac.kcmpnm,sac.khole)=(str(tr.stats.station), str(tr.stats.network), str(tr.stats.channel), str(tr.stats.location))
        (sac.stla,sac.stlo,sac.stel) = (tr.stats.station_coordinates.latitude, tr.stats.station_coordinates.longitude, tr.stats.station_coordinates.elevation)

        ev=tr.stats.event_origin
        time=ev.time
        # sac depth is in km
        sac.evla, sac.evlo, sac.evdp, sac.mag=ev.latitude, ev.longitude, ev.depth/1000., tr.stats.event_mag.mag
        sac.evla, sac.evlo, sac.evdp, sac.mag=ev.latitude, ev.longitude, ev.depth/1000., tr.stats.event_mag.mag
        # sac uses millisec while obspy uses microsec.
        sac.nzyear,sac.nzjday,sac.nzhour,sac.nzmin,sac.nzsec,sac.nzmsec=time.year,time.julday,time.hour,time.minute,time.second,time.microsecond/1000
        sac.o=0.
        sac.b=tr.stats.starttime-time # this is very important!!
        sac.kevnm=str(time)
        # dip is from horizontal downward; inc is from vertical downward
        # in SAC component "incidence angle" relative to the vertical 
        sac.cmpaz,sac.cmpinc=tr.stats.cmpaz,tr.stats.dip+90
        sac.gcarc,sac.dist,sac.az,sac.baz= tr.stats.gcarc,tr.stats.distance/1000,tr.stats.azimuth,tr.stats.back_azimuth
        # traveltimes
        sac.a=tr.stats.Parr.arrival_time; sac.ka='P' # cannot add S time because user1 is assigned to ray parameter
        # the ray parameter required by hk code is in sin(th)/v
        (sac.user0, sac.user1) = (tr.stats.Parr.rayp/radiusOfEarth, tr.stats.Sarr.rayp/radiusOfEarth)
        # add sac header to tr.stats
        tr.stats.sac=sac
        # set sac file name
        tr_name=write_dir+'/'+tr.stats.station+'.'+tr.stats.network+'.'+tr.stats.location+'.'+tr.stats.channel+ext+'.sac'
        tr.write(tr_name,format='SAC')
        if verbose:
            print('Writing sac file ...'+tr_name)
                  

def plot_stream_with_picks(str1,title='',add_accept_button=True):
    """
        plot the seismograms from a stream object
        add vertical lines indicating origin time, P and S picks to figure
        note the axes number go from bottom to top in (Z,N,E)
        while the trace plotting number goes from bottom to top (E,N,Z)
    """
    fig = plt.figure()
    str1.plot(fig=fig,show=False)
  # figure out how to add the button with these plotted stream data.

    # loop over subfigures to add picks
    for (j,tr) in enumerate(str1):
        ax=fig.axes[2-j]
        ev_time=tr.stats.event_origin.time
        ax.axvline(date2num(ev_time.datetime),lw=2)
        Parr, Sarr=tr.stats.Parr.time, tr.stats.Sarr.time
        ax.axvline(date2num(Parr.datetime),lw=1,ls='dashed',c='g')
        ax.axvline(date2num(Sarr.datetime),lw=1,ls='dashed',c='m')
       
        if j==0:
        # add picking info on the bottome figure
            text_y=min(tr.data)+(max(tr.data)-min(tr.data))*0.85
            ax.text(date2num(Parr.datetime),text_y,'P')
            ax.text(date2num(Sarr.datetime),text_y,'S')
            ax.text(date2num(ev_time.datetime),text_y,'O')
        elif j==2:
        # put the event information on the top graph
            ax.text(0.48, 0.85, 'Event #'+title+': '+str(ev_time.datetime),transform=ax.transAxes)
            ax.text(0.55, 0.7, 'gcarc='+str(int(tr.stats.gcarc))+'; mag='+str(tr.stats.event_mag.mag),transform=ax.transAxes)
    if add_accept_button:  # is there a better way of writing this?
        # Add an axes to the current figure and make it the current axes.
        # A new axes is added with dimensions rect in normalized (0, 1) units
        ax_accept = plt.axes([0.73, 0.05, 0.1, 0.075])
        bn_accept = Button(ax_accept, 'Accept', color='tomato', hovercolor='green')
        bn_accept.on_clicked(lambda x: accept_event(x, event_id=title))
        ax_close = plt.axes([0.85, 0.05, 0.1, 0.075])
        bn_close = Button(ax_close, 'Close', color='cornflowerblue', hovercolor='green')
        bn_close.on_clicked(close_fig)
    plt.show()

def accept_event(event, event_id='999',verbose=True):
    global accepted_event_list
    if verbose:
        print('Accept event #', event_id)
    try:
        iev=int(event_id)
    except ValueError:
        print('Not an integer', event_id)
    accepted_event_list.append(iev)
    plt.close()

def close_fig(event):
    plt.close()

#=== below are temprary codes that will be deleted in the future ==============
                #     ## display snr (needs fixing and should be after some processing, remean and retrend?)
                    # ts=tr.stats.starttime; t0=evt.origins[0].time
                    # surf_length=tr.stats.distance/1000.*3./10 # properly define this.
                    # t1=Parr; t2=Sarr; t3=Sarr+surf_length
                    # dt=tr.stats.delta
                    # n0=int((t0-ts)/dt); n1=int((t1-ts)/dt)
                    # n2=min(int((t2-ts)/dt),tr.stats.npts-100); n3=min(int((t3-ts)/dt),tr.stats.npts)
                    # print(ts,n0,n1,n2,n3,tr.stats.npts)
                    # pnl_snr=abs(tr.data[n1:n2]).mean()/abs(tr.data[n0:n1]).mean()
                    # print('Pnl SNR    = %4.1f' %  (pnl_snr))
      
                    # surf_snr=abs(tr.data[n2:n3]).mean()/abs(tr.data[n0:n1]).mean()
                    # print('Surf SNR %2d = %4.1f' % (j, surf_snr))
                    # print('max data %f' %(max(tr.data)))
       # if k not in trace_fini_set:  # unprocessed horizontal component
       #                          trace_k_list=[]
       #                          for l in trace_full_set.difference(trace_fini_set) :
       #                              if l != k:  # find the corresponding trace
       #                                  if trace.stats.channel[0:2] == seismo[l].stats.channel[0:2] and trace.stats.location==seismo[l].stats.location:
       #                                      trace_k_list.append(l); trace_fini_set.add(l)
       #                          if len(trace_k_list) !=2: 
       #                              sys.exit('Error number of components for trace: '+trace.stats.station+', '+trace.stats.channel)
       #                          trace2=seismo[trace_k_list[0]]; trace3=seismo[trace_k_list[1]] 
       #                          if len(trace.data) != len(trace2.data) or len(trace.data) != len(trace3.data):
       #                              print('**** Trace length different, skip trace ****',  len(trace.data),  len(trace2.data),  len(trace3.data)); skip_event=True; continue


        # tr.stats.gcarc=locations2degrees(chan.latitude, chan.longitude, evt.origins[0].latitude, evt.origins[0].longitude) # this is in deg; #kilometers2degrees(tr.stats.distance)
                  # 51.4455 176.7083 38.540142 28.633942 85.78551347291221 19.29817039741116 (Checked to be correctin  both gps2dist and taup_geo.calc_az_dist)
                  #51.4455     176.7083     38.54014   28.63394 1.937192e+01 (from SAC is a bit off, not sure why)
                 # tr.stats.back_azimuth=1.937192e+01
               # print('reset back azimuth to 1.937192e+01')
#(tr.stats.distance,tr.stats.azimuth,tr.stats.back_azimuth)=gps2dist_azimuth(evt.origins[0].latitude, evt.origins[0].longitude,chan.latitude, chan.longitude,f=0)
            #   str1 traces are ordered as [E,N,Z], need to check if always true
#
#print('back azimuth ',baz)dd
    #baz=1.937192e+01
    #print('reset baz', baz)

