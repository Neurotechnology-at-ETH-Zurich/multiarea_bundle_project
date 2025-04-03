# JRCLUST-4.1.1

Added some features/bug-fixes into the JRCLUST-4.1.0 spike-sorting algorithm, which is not maintained anymore. The last official version 4.1.0 can be found here: https://github.com/JaneliaSciComp/JRCLUST

Added features/bug-fixes:
-Trace View can now show band-pass filtered data
-The color of the spikes belonging to a cluster shown in the Trace View matches the color of the same cluster in the Waveform View

I would suggest the following workflow to prevent the noise and artifacts from skewing the clustering process: 
-Run "detect-sort" on your data 
-Open manual curation window
-Deleting the noise/artifact clusters
-Run the prep_for_recluster.m 
-Run "sort" or "recluster" (the latter does not recompute rho and delta, so may not be very accurate since those parameters were calculated with 
the noise spikes included)

JRCLUST was originally developed by [James Jun](https://www.simonsfoundation.org/team/james-jun/) and is currently maintained by [Vidrio Technologies](https://vidriotechnologies.com).

