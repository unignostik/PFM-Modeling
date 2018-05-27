%
% Author: Chase Tillar
% Build Date: 05/26/2018
% Description: Model the PFM elution curve assuming retardation of the
% tracer on the PFM sorbent is goverend by linear partitioning processes.
% The elution curve represents the dimensionless mass remaining in the PFM
% as a function of time
%

clc(); % clear console

% PFM dimensions
rPFM = 1; % radius of PFM (cm)
zPFM = 1; % height of PFM (cm)
sMax = 10; % number of streamtubes in PFM
dy = rPFM/sMax; % width of stream tube (cm)

% sorbent matrix properties 
v = 2.00; % seepage velocity (cm/hr)
p = 0.63; % porosity of sorptive matrix
b = 1.00; % thickness of sorptive matrix (cm)
pb = 0.52*1000; % bulk density of sorptive matrix (mg/mL)

% time constraints
tInit = 0; % initial time of deployment (hr)
tFinal = 10; % total time of deployment (hr)
dt = .5; % time increment (hr)

% tracer properties
cInit = 1000; % inital concentration (mg/L)
kd = 1.20; % tracer-sorbent partitioning coefficient (mL/mg)
rd = 1 + (pb*kd/p); % retardation factor

mInitPFM = 0;

for t=0:dt:tFinal % iterate through time by increments of dt
    
    mRPFM = 0; % reset mass in PFM at time t
    xf = v*t/rd; % shockfront position seperating two regions
    
    for tube=1:sMax % iterate through each streamtube 
        
        y = (tube-1)*dy; % distance of streamtube from center of PFM
        xd = 2*(rPFM.^2 - y.^2).^(1/2); % length of given streamtube
        
        mInit = (cInit/1000)*p*zPFM*dy*xd; % initial mass in streamtube
        
        if xf < xd % position inside stream tube
            mRtube = (cInit/1000)*p*zPFM*dy*(xd - xf); % mass remaining from streamtube
        else
            mRtube = 0;
        end 
        
        mInitPFM = mInitPFM + mInit; % intial mass in PFM
        mRPFM = mRPFM + mRtube; % mass remaining in PFM at time t 
        dMR = (mRPFM/mInitPFM); % dimensionless mass remaining in PFM
    end 
    
    % display data
    disp(t); 
    disp(dMR);
end
