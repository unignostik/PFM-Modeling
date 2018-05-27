%
% Author: Chase Tillar
% Build Date: 05/26/2018
% Description: Model the PFM elution curve assuming retardation of the
% tracer on the PFM sorbent is goverend by Fruendlich partitioning processes.
% The elution curve represents the dimensionless mass remaining in the PFM
% as a function of time
%

clc(); % clear console

% PFM dimensions
rPFM = 1; % radius of PFM (cm)
zPFM = 1; % height of PFM (cm)
sMax = 10; % number of streamtubes in PFM
dy = rPFM/sMax; % width of streamtube (cm)

% sorbent matrix properties 
v = 2.00; % seepage velocity (cm/hr)
p = 0.63; % porosity of sorptive matrix
b = 1.00; % thickness of sorptive matrix (cm)
pb = 0.52*1000; % bulk density of sorptive matrix (mg/mL)

% time constraints
tFinal = 10; % total time of deployment (hr)
dt = .5; % time increment (hr)

cInit = 1000; % inital concentration (mg/L)
mInitPFM = 0; % variable declaration for initial mass in PFM (mg)

for kf=1.20:0.25:1.50 % vary partitioning coefficient (mL/mg)
    for m=0.05:0.05:1.00 % vary Fruendlich exponent 
    rd = (p+pb*kf*cInit.^(m-1))/p; % retardation factor
    
        for t=0:dt:tFinal % iterate through time by increments of dt
            mRPFM = 0; % reset mass in PFM at time t
            vc = p*v/(p+pb*kf*m*cInit.^(m-1)); % fluid velocity at given concentration
            xb = vc*t; % position of fluid at time t
            
            for tube=1:sMax % iterate through each streamtube 
                y = (tube-1)*dy; % distance of streamtube from center of PFM (cm)
                xd = 2*(rPFM.^2 - y.^2).^(1/2); % length of given streamtube (cm)
                mInit = (cInit/1000)*p*zPFM*dy*xd; % initial mass in streamtube (mg)
                
                if xb<xd % position inside stream tube
                     % mass remaining in streamtube (mg)
                else
                     % mass remaining in streamtube (mg)
                end
                
                mInitPFM = mInitPFM + mInit; % initial mass in PFM (mg)
                mRPFM = mRPFM + mRtube; % mass remaining in PFM at time t (mg)
            end 
        dMR = (mRPFM/mInitPFM); % dimensionless mass remaining in PFM
        % display data
        disp(t); 
        disp(dMR);
        end
    end
end
