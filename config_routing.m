% This file lets the end user configure DICOM routing for their site-specific implementation of the Pro-IGNITE trial
% Routing info is defined as ordinary MATLAB variables that will be loaded at runtime

% DICOM destination info is defined in the MATLAB structure: destinations
% You can define multiple destinations like so:
% destinations(1) % contains networking info for first destination
% destinations(2) % contains networking info for second destination
% destinations(3) % contains networking info for third destination
% 
% Each DICOM destination object is defined by 4 parameters:
% ae_title: AE title of the thing that will send DICOMs
% called_ae_title: AE title of the thing that will be triggered to recieve/handle arriving DICOMs
% ip_address: IP address of the machine that will receive DICOMs
% port: Port that will receive DICOMs


% DICOM routing parameters
destinations(1).ae_title = 'PROIGNITE';
destinations(1).called_ae_title = 'ORTHANC';
destinations(1).ip_address = '169.228.56.133';
destinations(1).port = '4242';



