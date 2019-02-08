"""
Simulation pipeline to generate polarized power spectra from PAPER visibilities
Author: Ridhima Nunhokee
"""

Dependencies:
   1) Healpix
   2) Aipy

1) Generate polarized foregrounds using a source catalogue or healpix map:
      python genPolForegrounds.py <catalogue> --jd=<julian_dates> --start=<start_frequency> --stop=<stop_frequency> --chan=<nchan> --xpol=<beam_model_xx> --ypol=<beam_model_yy>
      
     Arguments:
        - source catalogue in text format or healpix map
     Options:
        -julian_dates   : Array containing julian dates
        -start_frequency: Starting frequency in Hz
        -stop_frequency : Stopping frequency in Hz
        -chan           : Number of frequency channels
        -beam_model_xx  : Beam model for xx polarization
        -beam_model_yy  : Beam model for yy polarization

2) Generate per baseline visibility. It genereates visibility for only one baseline at a time
      python genVisibility.py --inname=<input_name> --start=<start_frequency> --stop=<stop_frequency> --chan=<nchan> --filename=<output> --save --stokes --jd=<julian_date> -C <calfile> --radec=<radec>  -a <bl>

      Options:
         -julian_dates   : Array containing julian dates
         -start_frequency: Starting frequency in Hz
         -stop_frequency : Stopping frequency in Hz
         -chan           : Number of frequency channels
         -input_name     : Input name of npz file containing the foregrounds (output of genPolForegrounds.py)
         -output         : Name of output file
         -bl             : Baseline e.g 0_1
         -calfile        : File containg array positions e.g psa819_v008 (exclude .py extension)

3) Combine uvfiles for more that one baselines into a single file
       python combineUV.py <uvfiles> -o <output_file>

       Arguments:
          -per baseline uv files (*.uv)
       Options:
          -output_file : Name of output file
    You can check the output visibilities by plotting the visibilities using aipy's script plot_uv.py.

4) Extracts individual polarization from the combine uvfile
       python pull_pol.py <uvfile> -p <pol>
      
       Arguments:
          -uvfiles (output of combineUV.py)
       Options
          -pol : Polarization (can be I,Q,U,V or xx,xy,yx,yy) depending on what you have simulated

5) Maps (I,Q,U,V) visibility to (xx,yy,xy,yx). For some reason aipy does not want to work with I,Q,U,V so we change the polarization label if the visibilties to (xx,yy,xy,yx). Note that the visibilities are still intact, that is they are for Stokes (I,Q,U,V).
     python mapStokes2xy.py <uvfiles> -p <pol>
     
     Arguments:
       -uvfiles: uvfiles (output of pull_pol.py)

     Options:
       -pol : Polrization (xx,xy,yx,yy)

6) Constructs polarized power spectra, output an npz file
      python genPowerSpectra.py <uvfile> -C <calfile> --jd=<julian_dates> -p <pol> --bw=<frequency bandwidth> --f0=<center frequency>

       Arguments:
         -uvfiles: uvfiles (output of pull_pol.py)

     Options:
         -pol                 : Polarization (xx,xy,yx,yy)
         -julian_dates        : Npz file with the julian dates
         -calfile             : Calfile with the antenna positions
         -frequency bandwidth : Frequncy Bandwidth of the power spectra in Hz
         -center frquency     : Frequency in Hz at which the power spectra is centered

7) Plots waterfall plot of generated power spectra 
      python plotPowerSpectra.py <uvfiles> -p <pol> -C <calfile> --f0=<center frequency> --buff=<buffer>
      
      Arguments:
         -uvfiles: uvfiles

      Options:
         -pol            : Polarization (xx,xy,yx,yy)         
         -calfile        : Calfile with the antenna positions
         -center frquency: Frequency in Hz at which the power spectra is centered
         -buffer         : Buffer in ns
           
  
