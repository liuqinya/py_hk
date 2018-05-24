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
from sets import Set

def sac_header_from_event_station(evt,  net,  sta,  phases=('P'), model1D='ak135'):
    
    # note   sac.b=start time - evt.origins[0].time and sac.e needs to be set when the trace has been acquired
    sac=AttribDict()
    # station location: obspy missing station burial information!!
    sac.stla=sta.latitude; sac.stlo=sta.longitude; sac.stel=sta.elevation; sac.stdp=0.
    sac.kstnm=str(sta.code) 
    sac.knetwk=str(net.code)
    # component info set in sac_header_from_channel: sac.kcmpnm=channel.code ; cmpaz, cmpinc,
    # distance/ azimuth
    gcarc = locations2degrees(evt.origins[0].latitude, evt.origins[0].longitude, sta.latitude, sta.longitude)
    azi_baz = gps2dist_azimuth(evt.origins[0].latitude, evt.origins[0].longitude, sta.latitude, sta.longitude)
    sac.dist, sac.az,  sac.baz,  sac.gcarc=  azi_baz[0]/1000.0, azi_baz[1], azi_baz[2], gcarc 
    # event location
    sac.evla=evt.origins[0].latitude; sac.evlo=evt.origins[0].longitude; sac.evdp=evt.origins[0].depth/1000.
    sac.mag=evt.magnitudes[0].mag; time=evt.origins[0].time
    sac.kevnm=str(time)
    sac.nzyear,  sac.nzjday,  sac.nzhour,  sac.nzmin,  sac.nzsec,  sac.nzmsec=time.year, time.julday, time.hour, time.minute, time.second,  time.microsecond/1000
    sac.o=0
    model = TauPyModel(model=model1D)
    if len(phases) > 0:
        arrivals=model.get_travel_times(distance_in_degree=gcarc,source_depth_in_km=evt.origins[0].depth/1000.,phase_list=phases)
        for i, arr in enumerate(arrivals):
            sac['t'+str(i)]=arr.time
            sac['kt'+str(i)]=arr.name
            sac['user0']=arr.ray_param/6371  # Earth's radius
            if arr.name == 'P':
                sac['a']=arr.time
    else:
        arrivals=None
    return sac, arrivals
    
def sac_header_from_channel(chan,  sac):
    sac.kcmpnm=str(chan.code) # convert unicode
    sac.cmpaz=chan.azimuth
    sac.cmpinc=chan.dip+90 # dip is from horizontal downward; inc is from vertical downward
    sac.khole=str(chan.location_code)
    
class event_search_par:
    def __init__(self,stime=UTCDateTime('2016-06-01'),etime=UTCDateTime('2017-08-01'), mag_type='Mw',min_mag=6.0,max_mag=10.0, search_center=(55,-115), \
        min_search_radius=30,max_search_radius=95, min_dep=0.,max_dep=50.):
        self.magnitude_range= AttribDict({'mag_type': mag_type,'min_mag':min_mag,'max_mag':max_mag})
        self.time_range= AttribDict({'stime':stime,'etime':etime})
        self.distance_range= AttribDict({'center_lat': search_center[0],'center_lon':search_center[1], 'min_radius_in_deg':min_search_radius, 'max_radius_in_deg':max_search_radius})  # Alberta
        self.depth_range= AttribDict({'min_depth':min_dep,  'max_depth':max_dep})
        
    def write(self):
        # write event search parameters
        output='Event search criteria: \n    with '+self.magnitude_range['mag_type'] +' between ['+str(self.magnitude_range['min_mag'])+', ' \
        +str(self.magnitude_range['max_mag'])+']\n    time between ['+str(self.time_range['stime'])[0:23]+', '+str(self.time_range['etime'])[0:23] \
        +']\n    distance within ['+str(self.distance_range['min_radius_in_deg'])+', '+str(self.distance_range['max_radius_in_deg']) \
        +'] deg from center point ['+ str(self.distance_range['center_lat'])+', '+str(self.distance_range['center_lon'])+'] deg\n    depth between [' \
        +str(self.depth_range['min_depth'])+', '+str(self.depth_range['max_depth'])+'] km\n' 
        print(output)

class station_search_par:
    def __init__(self, network='*', station='*', channel='BH*,HH*',min_lat=45, max_lat=60, min_lon=-125,max_lon=-110,above_line=(False,-115,49,-120,53,True)):
        self.box_range= AttribDict({'min_lat':min_lat, 'max_lat':max_lat, 'min_lon': min_lon, 'max_lon': max_lon})
        self.channel=channel
        self.network=network
        self.station=station
        
        self.include_line=above_line[0]
        self.line_par= AttribDict({'lon1':above_line[1], 'lat1':above_line[2], 'lon2': above_line[3], 'lat2': above_line[4], 'above': above_line[5], 'aa': 0., 'bb': 0.})
   
    def write(self):
        # write station search parameters
        output='\nStation search criteria:\n    within box range (lat;lon): ('+str(self.box_range['min_lat'])+', '+str(self.box_range['max_lat'])+'; ' \
        +str(self.box_range['min_lon'])+', '+str(self.box_range['max_lon'])+')\n'+'    for channel '+self.channel+'\n    for network '+self.network \
        +'\n    for station '+self.station
        if (self.include_line):
            above_line_dict=['below','above']
            output+='\n    '+above_line_dict[self.line_par['above']]+' line: lat = '+str(self.line_par['aa'])+' * lon + '+str(self.line_par['bb'])
        print(output)

class waveform_search_par:
    def __init__(self,gcarc_min=0.,gcarc_max=180.,phase_select='P',pre_length=10,length=3600,model1d='ak135'):
        # mulitple phases split with ','
        self.dist= AttribDict({'gcarc_min':gcarc_min,  'gcarc_max':gcarc_max})
        self.windows= AttribDict({'phase_select': phase_select.split(','),  'length_before': pre_length,  'length_after': length,  '1d_ref_model': model1d})
        if self.windows['length_before'] < 0:
            sys.exit('set phase window pre_length > 0')
        
    def write(self):
        output='\nWaveform search criteria: \n'+'    within gcarc range: ('+str(self.dist['gcarc_min'])+', '+str(self.dist['gcarc_max'])+')\n' \
        +'    for phases '+self.windows['phase_select'][0]+' between [-'+str(self.windows['length_before'])+', '+str(self.windows['length_after'])+']'
        print(output)

class recdata:
    def __init__(self,client_name='IRIS',data_dir='.'):
        self.client_name=client_name
        print('Initialize IRIS client... (check for event, station, dataselect services)')
        self.client = Client(self.client_name, debug=False)
        print(self.client)
        self.data_dir='.'
        self.event_par=event_search_par()
        self.station_par=station_search_par()
        self.waveform_par=waveform_search_par()
        self.sacfile={}# a dictionary that records a table of net/station/event/channel/sacfile

    def get_events(self,stime=UTCDateTime('2016-06-01'),etime=UTCDateTime('2017-08-01'), mag_type='Mw',min_mag=6.0,max_mag=10.0, \
            search_center=(55,-115),min_search_radius=30,max_search_radius=95, min_dep=0.,max_dep=50.,event_list_txt=None): # default values for Alberta
        
        event=event_search_par(stime=stime,etime=etime, mag_type=mag_type,min_mag=min_mag,max_mag=max_mag, \
            search_center=search_center,min_search_radius=min_search_radius,max_search_radius=max_search_radius, min_dep=min_dep,max_dep=max_dep)
        self.event_par=event
        self.event_par.write()

        print('Request event catalog from '+self.client_name+' ...')
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
            print('    No event within range')
            sys.exit()
        else:
            print('    Total number of events requested: '+str(len(self.event_catalog)))

        # output event list in text file
        if (event_list_txt):
            print('    Writing events to file '+event_list_txt+' ...')
            with open(event_list_txt,'w') as fp:
                for one_event in self.event_catalog:
                    evt=one_event.origins[0]
                    mag=one_event.magnitudes[0]
                    fp.write("xxxxx,%s,%.3f,%.4f,%.1f,xxxx,yyyyy,%s,%.1f\n" % 
                             (str(evt.time)[0:23].replace('T',' '),  
                              evt.latitude,evt.longitude,evt.depth/1000.,self.event_par.magnitude_range['mag_type'],mag.mag))

    def event_catalog_plot(self):
        self.event_catalog.plot()


    def get_stations(self,network='*', station='*', min_lat=45, max_lat=60, min_lon=-125,max_lon=-110,channel='BH*,HH*',above_line=(False,-115,49,-120,53,True), station_list_txt=None):
        # for more channel information, http://scedc.caltech.edu/station/seed.html
        # HH high broadband, triggered data
        
        sta=station_search_par(network=network, station=station, min_lat=min_lat, max_lat=max_lat, min_lon=min_lon,max_lon=max_lon,channel=channel,above_line=above_line)
        if sta.include_line:
            sta.line_par['aa']=(sta.line_par['lat2']*1.0-sta.line_par['lat1'])/(sta.line_par['lon2']-sta.line_par['lon1'])
            sta.line_par['bb']=sta.line_par['lat1']-sta.line_par['aa']*sta.line_par['lon1']
        self.station_par=sta
        self.station_par.write()

        print('\nSearch stations from '+self.client_name+' ...')
        #print(self.station_par.network,  self.station_par.station,self.station_par.channel, self.station_par.box_range['min_lat'], self.station_par.box_range['max_lat'], self.station_par.box_range['min_lon'],  self.station_par.box_range['max_lon']  )
        sta_inv=self.client.get_stations(network=self.station_par.network,  station=self.station_par.station, 
            minlatitude=self.station_par.box_range['min_lat'],maxlatitude=self.station_par.box_range['max_lat'],
            minlongitude=self.station_par.box_range['min_lon'],maxlongitude=self.station_par.box_range['max_lon'],
            channel=self.station_par.channel,level='channel') 

        # need a deep copy here!
        self.station_inventory=copy.deepcopy(sta_inv)
        num_net=0
        if (station_list_txt):
            print('    Writing network/station list file '+station_list_txt+'...')
            fp=open(station_list_txt,'w')
        
        
        # loop over networks
        idel=0
        for i, net in enumerate(sta_inv): 
            sta_lat=[]; sta_lon=[]; jdel=0
            # loop over all stations in this network
            for j, sta in enumerate(net):
                slat=sta.latitude
                slon=sta.longitude
                include_sta=False
                if (self.station_par.include_line):
                    if (self.station_par.line_par['above'] and slat > self.station_par.line_par['aa']*slon+self.station_par.line_par['bb'] or not self.station_par.line_par['above'] and slat < self.station_par.line_par['aa']*slon+self.station_par.line_par['bb']):
                        include_sta=True
                else:
                    include_sta=True
                  # DEBUG                print(i,j,net.code,net[j].code,include_sta,i-idel,j-jdel,len(self.station_inventory.networks[i-idel]))
                if include_sta:
                    sta_lat.append(slat)
                    sta_lon.append(slon)
                    if (station_list_txt):
                        cc='';fp.write('%s,%s,%.4f,%.4f,%.1f,%s,%s,%s\n' % (net.code,sta.code,slat,slon,sta.elevation,sta.start_date,sta.end_date, '/'.join(list(set([cc+str(trace.code) for trace in sta ] ))))) # comp/comp
                else:
                    del self.station_inventory.networks[i-idel].stations[j-jdel]
                    jdel+=1
            # for this network
            if (len(sta_lat) > 0):
                num_net=num_net+1
                #DEBUG print(net.code,len(net))
            else:
                del self.station_inventory.networks[i-idel]
                idel+=1

        print('    Number of networks/stations selected: '+str(num_net)+' '+str(len(self.station_inventory)))
        if station_list_txt:
            fp.close()


    def station_inventory_plot(self):
        # detailed station/network distribution plot
        # setup stereographic basemap.
        # lat_ts is latitude of true scale.
        # lon_0,lat_0 is central point.
        lon0=(self.station_par.box_range['min_lon']+self.station_par.box_range['max_lon'])/2
        lat0=(self.station_par.box_range['min_lat']+self.station_par.box_range['max_lat'])/2
        lats=lat0
        map = Basemap(projection='stere',lon_0=lon0,lat_0=lat0,lat_ts=lats,\
                      llcrnrlat=self.station_par.box_range['min_lat'],urcrnrlat=self.station_par.box_range['max_lat'],\
                      llcrnrlon=self.station_par.box_range['min_lon']-1,urcrnrlon=self.station_par.box_range['max_lon']+2,\
                      rsphere=6371200.,resolution='h',area_thresh=10000)
        # draw coastlines, country boundaries, fill continents.
        map.drawcoastlines()
        map.drawstates()
        map.drawcountries()
        # draw the edge of the map projection region (the projection limb)
        #map.drawmapboundary(fill_color='aqua')
        # draw lat/lon grid lines every 30 degrees.
        map.drawparallels(np.arange(self.station_par.box_range['min_lat'],self.station_par.box_range['max_lat'],5),labels=[True,False,False,True])
        map.drawmeridians(np.arange(self.station_par.box_range['min_lon'],self.station_par.box_range['max_lon'],5),labels=[True,False,False,True])
        map.etopo()

        # set up maximum number of colors for network 3*3*2
        r=np.linspace(0,1,num=3,endpoint=True);g=r;b=np.linspace(0,1,num=2,endpoint=True)
        rv,gv,bv=np.meshgrid(r,g,b)
        nv=rv.size
        rv=rv.reshape(nv,-1);gv=gv.reshape(nv,-1); bv=bv.reshape(nv,-1)   

        # plot by network
        for i, net in enumerate(self.station_inventory): 
            sta_lat=[]; sta_lon=[] 
            sdate=UTCDateTime('2100-01-01'); edate=UTCDateTime('1970-01-01')
            for sta in net:
                sta_lat.append(sta.latitude)
                sta_lon.append(sta.longitude)
                sdate=min(sdate,sta.start_date); edate=max(edate,min(sta.end_date,UTCDateTime.now()))
            x,y=map(sta_lon,sta_lat)
            if (i >= nv):
                sys.exit('Not enough colors for plotting')
            map.scatter(x,y,30,marker='^',color=(rv[i],gv[i],bv[i]),label=str(net.code)+' '+str(sdate)[2:7]+'/'+str(edate)[2:7])
            
        plt.legend(bbox_to_anchor=(1.5, 1))
        plt.show()
        
    def get_waveforms(self,phase_select='P',pre_length=10,length=3600,model1d='ak135', \
        remove_response=False,filter_seismogram=False, f_low=0.001, f_high=10., rotate=False, rotate_to_rt_in_sac=False, check_resp_same=False, snr=False, pkl=False, pickle_name='rec.pkl'):
    # now rotation has been fully tested against sac rotation
        # set up waveform parameters
        wf_par=waveform_search_par(gcarc_min=self.event_par.distance_range['min_radius_in_deg'],gcarc_max=self.event_par.distance_range['max_radius_in_deg'], \
            phase_select=phase_select,pre_length=pre_length,length=length,model1d=model1d)
        self.waveform_par=wf_par
        self.waveform_par.write()
        ext='.d' if remove_response else ''
        signal_to_noise_ratio=0.
        wf_stream=Stream()
        # loop over networks, stations
        for net in self.station_inventory:
            for sta in net:
                sdate=sta.start_date; edate=sta.end_date
                sta_folder=self.data_dir+'/'+str(net.code)+'.'+str(sta.code)
                if not os.path.exists(sta_folder):
                    os.makedirs(sta_folder)
                # loop over events
                for evt in self.event_catalog:
                    skip_event=False; resp_dict=AttribDict()
                    evtorg=evt.origins[0]  # take only the first origin and magnitude
                    evt_starttime=UTCDateTime(evtorg.time); #evt_dep=evtorg.depth/1000. # convert to km
                    evnm=str(evt_starttime)
                    # download seismograms for events within station active window
                    if (evt_starttime < sdate or evt_starttime > edate):
                        print('event not within station active time window ['+str(sdate)+',  '+str(edate)+'] : '+evnm+' for '+sta.code)
                        continue
                    # download seismogram within the given dist range
                    sac_header, arrivals=sac_header_from_event_station(evt, net, sta,  self.waveform_par.windows['phase_select'], model1D=self.waveform_par.windows['1d_ref_model'])
                    if sac_header.gcarc < self.waveform_par.dist['gcarc_min'] or sac_header.gcarc > self.waveform_par.dist['gcarc_max']:
                        print('Station not within gcarc range ['+ str(self.waveform_par.dist['gcarc_min'] ) + ', ' + str(self.waveform_par.dist['gcarc_max'])+'] :'+str(sac_header.gcarc))
                        continue
                    t1=evt_starttime+arrivals[0].time-self.waveform_par.windows['length_before']
                    t2=evt_starttime+arrivals[0].time+self.waveform_par.windows['length_after']
                    try:    
                        seismo=self.client.get_waveforms(network=net.code, station=sta.code, channel=self.station_par.channel, location='*',  \
                            starttime=t1, endtime=t2, attach_response=True)
                    except:
                        print('No data available for '+net.code+' '+sta.code+' '+self.station_par.channel+' '+evnm+ ' within ['+str(t1)+', '+ str(t2)+']' \
                            +' for phase '+arrivals[0].name)
                        continue
                    print('****download data for '+net.code+' '+sta.code+' '+self.station_par.channel+' '+evnm+ ' within ['+str(t1)+' ,'+ str(t2)+'] for phase '+arrivals[0].name+': '+str(len(seismo)))
                    # merge sac files
                    try:
                        seismo.merge(fill_value='interpolate')
                    except:
                        print('Error merging trace for '+net.code+' '+sta.code+' '+self.station_par.channel+' '+evnm)
                        continue
                    if len(seismo) % 3 != 0:
                        sys.exit('Error number of traces not multiples of 3')
                    # add proper sac header ready for output
                    for trace in seismo:
                        sac=sac_header.copy()
                        # obtain channel info from station inventory 
                        channel=next((chan for chan in sta if (chan.code == trace.stats.channel and chan.location_code==trace.stats.location and chan.start_date <= trace.stats.starttime and chan.end_date >= trace.stats.endtime)), None)
                        if not channel:
                            sys.exit('Error finding proper channel info: '+net.code+' '+stat.code+' '+trace.channel)
                        sac_header_from_channel(channel, sac)
                        sac.b, sac.e=trace.stats.starttime-evtorg.time,  trace.stats.endtime-evtorg.time
                        trace.stats.sac = sac
                        trace.stats.back_azimuth=sac.baz
                        loc='.'+trace.stats.location
                        if not loc in resp_dict:
                            resp_dict[loc]=AttribDict()
                        resp_dict[loc][trace.stats.channel]=trace.stats.response
                    if check_resp_same:
                        for loc in resp_dict:
                            for m,chan in enumerate(resp_dict[loc]):
                                if m==0:
                                    r1=resp_dict[loc][chan]
                                elif m==1:
                                    r2=resp_dict[loc][chan]
                                else:
                                    r3=resp_dict[loc][chan]
                            if not (r1 == r2 and r2 == r3):
                                print('***Different instrument response on '+chan+loc); skip_event=True
                    
                    f_high=min(0.5*0.5/trace.stats.delta, f_high)
                    seismo.detrend("linear")
                    seismo.taper(max_percentage=0.05, type="hann")
                    if remove_response:
                        print('    remove instrument response ...')
                    # remove instrument response
                        pre_filter = [f_low*0.9, f_low, f_high,  f_high*1.1]
                        #seismo.detrend("linear")
                        #seismo.taper(max_percentage=0.05, type="hann")
                        seismo.remove_response(output='DISP', pre_filt=pre_filter, water_level=60)  #10^{-60/20}
                    if filter_seismogram:
                        print('    low-pass filter seismogram with f_high='+str(f_high)+' Hz')
                         #   seismo.filter('bandpass', freq1=f_low,  freq2=f_high, corners=4, zerophase=True)
                        seismo.filter('lowpass', freq=f_high, corners=2, zerophase=True) #only lowpass filter works fine with my version
                    file_base="%s.%s.%04d.%03d.%02d%02d%02d" % (net.code,  sta.code,  sac_header.nzyear,  sac_header.nzjday,  sac_header.nzhour,  sac_header.nzmin,  sac_header.nzsec)
                    file_base_list=file_base.split('.')
                    # rename sac files: following the exact naming convention from Wilber 3 requested sac files
                    sac_input={}; sac_output={}; trace_full_set=Set(range(len(seismo))); trace_fini_set=Set([])
                    if rotate and not rotate_to_rt_in_sac:  # rotate in obspy
                        for k, trace in enumerate(seismo):
                            if k not in trace_fini_set:  # unprocessed horizontal component
                                trace_k_list=[]
                                for l in trace_full_set.difference(trace_fini_set) :
                                    if l != k:  # find the corresponding trace
                                        if trace.stats.channel[0:2] == seismo[l].stats.channel[0:2] and trace.stats.location==seismo[l].stats.location:
                                            trace_k_list.append(l); trace_fini_set.add(l)
                                if len(trace_k_list) !=2: 
                                    sys.exit('Error number of components for trace: '+trace.stats.station+', '+trace.stats.channel)
                                trace2=seismo[trace_k_list[0]]; trace3=seismo[trace_k_list[1]] 
                                if len(trace.data) != len(trace2.data) or len(trace.data) != len(trace3.data):
                                    print('**** Trace length different, skip trace ****',  len(trace.data),  len(trace2.data),  len(trace3.data)); skip_event=True; continue
                                trace.data, trace2.data, trace3.data = rotate2zne(trace.data, trace.stats.sac.cmpaz, trace.stats.sac.cmpinc-90., \
                                trace2.data, trace2.stats.sac.cmpaz, trace2.stats.sac.cmpinc-90., \
                                trace3.data, trace3.stats.sac.cmpaz, trace3.stats.sac.cmpinc-90.)  # check if this is ok!
                                # change component names ahead of rotation
                                trace.stats.channel, trace2.stats.channel, trace3.stats.channel= trace.stats.channel[0:2]+'Z', trace2.stats.channel[0:2] +'N', trace3.stats.channel[0:2]+'E'
                                trace.stats.sac.cmpinc,  trace2.stats.sac.cmpinc, trace3.stats.sac.cmpinc = 0.,  90.,  90.  
                                # first correct baz, cmpinc, cmpaz in the sac header because they won't be corrected in rotate ...
                                trace.stats.sac.cmpaz,   trace2.stats.sac.cmpaz,   trace3.stats.sac.cmpaz = 0.,  (trace.stats.back_azimuth+180)%360,  (trace.stats.back_azimuth-90)%360  
                                trace_fini_set.add(k)
                        try:
                            seismo.rotate('NE->RT')  #this seismo should have proper R T components defined
                        except:
                            print('***** Error rotation, skip trace ******'); continue
                    if skip_event:  # no matching horizontal traces, skip to next event
                        continue
                    for trace in seismo:
                        sackey=net.code+'.'+sta.code+'.'+trace.stats.location
                        if not (sackey in self.sacfile):
                            self.sacfile[sackey]=AttribDict()
                        if not evnm in self.sacfile[sackey]:
                            self.sacfile[sackey][evnm]=AttribDict()
                            self.sacfile[sackey][evnm].rot_count=0
                        if not trace.stats.channel in self.sacfile[sackey][evnm]:
                          self.sacfile[sackey][evnm][trace.stats.channel]=AttribDict()
                        if snr:
                            n1=max(int((self.waveform_par.windows.length_before-5.0)*trace.stats.sampling_rate),  1.0)
                            n2=int(self.waveform_par.windows.length_before * trace.stats.sampling_rate) 
                            npts=trace.stats.npts
                            noise=sum(abs(trace.data[0:n1]))/n1; signal=sum(abs(trace.data[n2:npts]))/(npts-n2)  # one can change definition of signal [n2,npts]
                            signal_to_noise_ratio=signal/noise
                            trace.stats.sac.user1=signal_to_noise_ratio
                            self.sacfile[sackey][evnm][trace.stats.channel].snr=signal_to_noise_ratio
                    # write seismograms in sac and rename based on component
                    seismo.write(sta_folder+'/'+file_base+'.SAC', format='sac')
                    for k,  trace in enumerate(seismo):
                        sackey=net.code+'.'+sta.code+'.'+trace.stats.location
                        chan=trace.stats.channel
                        sacfile_base=file_base+'%02d'% (k+1)+'.SAC'
                        trace_base='.'.join(file_base_list[0:2])+'.'+trace.stats.location+'.'+trace.stats.channel+'.'+trace.stats.mseed.dataquality+'.'+'.'.join(file_base_list[2:7])+'.SAC'+ext
                        print('rename '+sta_folder+'/'+sacfile_base+' to '+ sta_folder+'/'+trace_base)
                        os.rename(sta_folder+'/'+sacfile_base,  sta_folder+'/'+trace_base)
                        
                        #wf_stream+=trace
                        if rotate_to_rt_in_sac:
                            loc=trace.stats.location
                            if not loc in sac_input:
                                sac_input[loc]=''; sac_output[loc]=''
                            if  trace.stats.channel[2] != 'Z':    # 1/2 or E/N components
                                sac_input[loc]=sac_input[loc]+' '+sta_folder+'/'+trace_base  # input one horizontal component
                                bname_list=trace_base.split('.'); a=list(bname_list[3]); 
                                if self.sacfile[sackey][evnm].rot_count == 0:
                                    a[2]=u'R'
                                elif self.sacfile[sackey][evnm].rot_count ==1:
                                    a[2]=u'T'
                                else:
                                    sys.exit('Error number of traces in '+trace_base)
                                self.sacfile[sackey][evnm].rot_count+=1
                                bname_list[3]=''.join(a)
                                trace_base='.'.join(bname_list)  # output R and T
                                sac_output[loc]=sac_output[loc]+' '+sta_folder+'/'+trace_base
                            with open('sac.cmd', 'w') as fp:
                                for loc in sac_input:
                                    fp.write('r %s\nrotate\nw%s\n' % (sac_input[loc], sac_output[loc]))
                                fp.write('quit\n')
                            #print('r %s\nrotate\nw%s\n' % (sac_input, sac_output))
                            if os.system('sac < sac.cmd >/dev/null') != 0:
                                sys.exit('Problem rotate seismograms in sac')
                                
                        self.sacfile[sackey][evnm][chan].file=sta_folder+'/'+trace_base
                        self.sacfile[sackey][evnm][chan].sac=trace.stats.sac   
        if pkl:
            with open(pickle_name, 'wb') as fp:
                pickle.dump(self, fp,  -1)
                         #only use the first location available:
  #                      if not self.sacfile[net.code][sta.code][str(evt_starttime)][trace.stats.channel].file:
   #                         self.sacfile[net.code][sta.code][str(evt_starttime)][trace.stats.channel].file=sta_folder+'/'+trace_base   
        #self.stream_data=wf_stream
    def cut_segment_from_phase(self, cut_file='rec_sac_file', snr=0.,before_length=30,after_length=90):
        for i,  phase in enumerate(self.waveform_par.windows.phase_select):
            print('Write sac cut input file '+cut_file+'.'+phase)
            print('Write hk package command lines\n')
            with open(cut_file+'.'+phase, 'w') as fp:
                for sackey in self.sacfile:
                    with open(sackey+'.'+phase+'.hk', 'w') as fpp:
                        for evnm in self.sacfile[sackey]:
                            for chan in self.sacfile[sackey][evnm]:
                                if chan[2] == 'Z': # need to add t file in future?
                                    t1=self.sacfile[sackey][evnm][chan].sac['t'+str(i)]-before_length
                                    t2=self.sacfile[sackey][evnm][chan].sac['t'+str(i)]+after_length
                                    zfile=self.sacfile[sackey][evnm][chan].file
                                elif chan[2]=='R':
                                    rfile=self.sacfile[sackey][evnm][chan].file
                            if self.sacfile[sackey][evnm][chan].snr > snr:
                                fp.write('cut  %.1f   %.1f\nr %s   %s\ncut off\n w %s %s\n'% (t1,  t2, zfile,  rfile,    zfile+'.'+phase,   rfile+'.'+phase))
                                fpp.write('hk/iter_decon -C-2/-10/80 -F1/3/-5 -N100 -T0.1 %s %s\n'%(zfile+'.'+phase, rfile+'.'+phase))
                                #print('write cut  %.1f   %.1f\nr %s   %s\ncut off\n w %s %s\n'% (t1,  t2, zfile,  rfile,    zfile+'.'+phase,   rfile+'.'+phase))
                fp.write('quit\n')
            
  
#        if trace.ts.sac.user1 > self.sacfile[net.code][sta.code][str(evt_starttime)][trace.stats.channel].snr:
#                              self.sacfile[net.code][sta.code][str(evt_starttime)][trace.stats.channel].snr=trace.stats.sac.user1
#                              self.sacfile[net.code][sta.code][str(evt_starttime)][trace.stats.channel].file=sta_folder+'/'+trace_base   

# ###################old code segments #######################################################
#               # download seismograms within depth and magnitude ranges
#                    if evt_dep > self.event_par.depth_range['max_depth'] or evt_dep < self.event_par.depth_range['min_depth']:
#                        print('event not within dep range ['+str(self.event_par.depth_range['min_depth'])+', '+str(self.event_par.depth_range['max_depth'])+'] : ' \
#                            +str(evt_dep))
#                        continue
#                    if evt.magnitudes[0].mag > self.event_par.magnitude_range['max_mag'] or evt.magnitudes[0].mag < self.event_par.magnitude_range['min_mag']:
#                        print('event not within magnitude range ['+str(self.event_par.magnitude_range['min_mag'])+', '+str(self.event_par.magnitude_range['max_mag'])+'] : ' \
#                            +str(evt.magnitudes[0].mag))
#                        continue

# request waveform data
#bulk.append((net.code, sta.code, '*', self.station_par.channel, t1, t2))
#wf_stream=self.client.get_waveforms_bulk(bulk, attach_response=True)

# bname_list[3]=reduce(lambda a, kv: a.replace(*kv), rep_pattern.iteritems(), bname_list[3]) 
#  rep_pattern={'N':'R', 'E':'T', '1':'R', '2':'T'}
# AUG 21 todo: snr goes into waveform_par, phase_select -> waveform_par, extra par for length,  anything else needed by the 
# AUG 22: check why rec.sacfile has only one entry!
# AUG 29: to do: clear how input phases names (as well as in sacfile) are treated because USER0 is expected in hk package
