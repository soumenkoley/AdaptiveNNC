This folder has all the programs necessary for implementing the numerically stabilized
fast transversal filter implementation of the Recursive Least square algorithm.

Data are downloaded in chunks of 1000 s from the Virgo data server, given it is establishd
that the dynamic Wiener filter implemented every 1000 s serves as a benchmak for dynamic
noise-cancellation.

main code to be run - runFRLSFull.m

parameters for FRLS code is set using - setFRLSStruct.m

general input-output parameters are set using - setParamsIO.m

some parameters include:
gpsStartTime, total signal lengh to be analyzed, filter cutoff frequencies, filter orders,
sampling frequency of target and reference channels etc.

parameters for the Wiener implmentation - setParamsWienerTest.m and setWienStruct.m

function for FTF RLS implementation - doSFTF.m (main program to be modified to change the algorithm)

The codes have been written to test a real time implementation, hence causality is maintained
Several functions like firFilt, firDecimate have been written keeping causality is mind

RefChaNames.txt - list of the NEB geophones used in the analysis

