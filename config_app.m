% Configuration file for Pro-IGNITE app
% The path path to this file (or similar) should be the 2nd argument to the entrypoint function pro_ignite.m

% Path to directory that will serve as a temporary staging area for incoming data
params.tmp_directory = '/home/tmp_pro_ignite';

% Time in hours to store data in staging area before deletion
params.lifetime_tmp_directory = 24;

% Container format 
% Typically 'docker' for deployed installations outside of UCSD
% Typically 'singularity' for local CMIG processing at UCSD
params.container = 'docker';

% Path to configuration file for RSI processing
params.rsi_config_file = '/home/ccconlin/work/pro-ignite/config_rsi.m'; 

% Path to configuration file for DICOM networking
params.routing_config_file = '/home/ccconlin/work/pro-ignite/config_routing.m';
