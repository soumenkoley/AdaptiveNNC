This folder has all the programs necessary for estimating the Wiener filter for a day of data.

buildCorrMat.m is used to estmate the average Rxx and Pxy matrices corresponding to a day of data

getDailyWFilt.m uses these correlation matrices to estimate the Wiener filter corresponding to that
day of data.

general input-output parameters are set using - setParamsIO.m

some parameters include:
gpsStartTime, total signal lengh to be analyzed, filter cutoff frequencies, filter orders,
sampling frequency of target and reference channels etc.

parameters for the Wiener implmentation - setParamsWienerTest.m and setWienStruct.m

RefChaNames.txt - list of the NEB geophones used in the analysis

