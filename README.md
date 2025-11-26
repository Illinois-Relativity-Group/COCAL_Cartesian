#created by Rohan Narasimhan 06/11/24

first softlink/copy all data files to COCAL_output (rnsflu_3D.las, rnsgrids_3D.las etc)
copy bhbhdisk.ct (or whatever colourtable you are using) to your ~/.visit directory (bhbhdisk is in visitconfigs/ctandxmlsbhdisks)

ALWAYS RUN PARAMS WHENEVER YOU OPEN A NEW TERMINAL

params file:

-all the initial settings are in params choose which data you are plotting
-set the size and resolution of the cartesian grid.
-run by typing . params in terminal


run_reader:

-this is the script which runs the interpolator and hdf5 conversion script, 99% of the time you will not need to change anything
- run by typing . run_reader.sh in terminal
- can also run as a process (in the background) using . run_reader.sh >&log.txt&
- this outputs data.txt (variable values), plot.3d (2d plotting), and plot_data.h5 (visit plotting)


visit_plot.py:

-this is the visualization script, YOU DONT RUN THIS SCRIPT
- if the image doesnt look great, the main things to change are Pseudo.min and Psuedo.max, and the iso setting file in visitconfigs/ctandxmlsbhdisks


run_visit_script:
-this is what you run to use visit_plot.py
-make sure the version of visit matches the one you have
-run by typing . run_visit_script.sh in terminal

------------------------------------------additional testing--------------------------------------

plot_2d.py
- plots an xy, xz, and yz slice at a slice value you can set
