This is the original HK1.3 package that are listed here for reference. The complete package can be downloaded from Prof. Lupei Zhu's website:

http://www.eas.slu.edu/People/LZhu/home.html  

1. Compile:
Follow instructions in README file and set GMT include/library directories first; Fix the bug related to BOOLEAN (into GMT_LONG), and the warning on fprintf by adding %d into the 'format'.
Add -DBYTE_SWAP to allow byte swapping after reading the sac binary files

2. Test:
 1) ./iter_decon -F1/3/-5 -N100 -C-2/-10/80 -T0.1 example/KUL.z example/KUL.[r,t]
to obtain receiver functions KUL.ri and KUL.ti.
-Fn/f1/f2/order: (n=1/2/3: apply filter to output/src/data; f1: low-pass filter; f2: (<0: pre-signal length; >0: corner freqs of the butterworth filter)
-Cmark/t1/t2: (windowing data [tmarkt+t1, tmark+t2], mark: -5(b), -3(o), -2(a), 0-9 (t0-t9)
-Nniter: number of iterations
-Ttaper: taper the trace with cosine window, taper=0-0.5
saclst a ka f example/KUL.[rtz]
To prepare these *.[rtz] files: download data at [-30, 80] around the P arrival time. Pick the first P arrival and put it into the header (a). Rotate into r,t,z components. What also needed are headers: dista, az, baz, gcarc, user0 (?)

2) ./k_stack -R20/60/1.5/2.0 -I0.5/0.01 -Gexample/hk.grd example/pp.*.ri
to generate example/hk.grd and example/hk.grd.var (variance) 
-Rz1/z2/k1/k2: crustal thickness and Vp/Vs grid
-Idz/dk: search grid spacing
-G output file name
These *.ri files have a, t8, user0, user1(?), user0=0.0418, 0.04929, 0.05789, 007825, 0.06738

3) ./grdmin -D example/hk.grd | ./hk_plt.pl > junk.ps
you can open the junk.ps file with gv junk.ps
hk_plt.pl has been updated to add 'gmtset MEASURE_UNIT inch' and change -JX4/2 to -JX8/4.

Data preparation: waveform data downloaded and organized by events and rotated into the form of event_ID/sta.evid.[r,t,z]. The ray parameter p
needs to be stored in SAC header USER0. IThe first P arrival be picked and stored in SAC header A; 
The time window can be normally chosen to be [-10,80] around t_P with a sampling rate of 10 sample/sec.
