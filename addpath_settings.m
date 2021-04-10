restoredefaultpath
parentpath = pwd;
%addpath([parentpath '/data/',system_type, '/']);
%addpath([parentpath '/figure/',system_type, '/']);
addpath([parentpath '/model/',system_type, '/']);
addpath([parentpath '/simulator/']);
addpath([parentpath '/ATLAS/']);
%addpath([parentpath '/misc/']);

datapath      = ['~/ATLASdata/',system_type, '/'];
%figurepath    = [parentpath '/figure/',system_type, '/fig/'];
if ~isfolder(datapath)
     mkdir(datapath)
end