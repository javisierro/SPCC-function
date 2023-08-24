## R code to compare two sound WAV files using spectrogram cross-correlation.

One sound file, note2 is compared to a reference sound, note2.
This comparison implies that note2 is overlaid onto note1 and the temporal offset can be delimited in maxtoff. The larger the temporal offset, the longer the computing time.
The temporal resolution in each step taken by the SPCC is set as t.resolution object. Increasing temporal resolution increases computing time although it will be more accurate.

The object 'crosscorr' is the SPCC curve with the correlation indexes calculated at each time offset. The classic SPCC usually takes the maximum point in this curve as the SPCC index, although alternative approaches are possible