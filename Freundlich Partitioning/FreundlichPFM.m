%
% Author: Chase Tillar
% Build Date: 05/26/2018
% Last Update: 07/12/2018
% Description: Model the PFM elution curve assuming retardation of the
% tracer on the PFM sorbent is goverend by Fruendlich partitioning processes;
% the elution curve represents the dimensionless mass remaining in the PFM
% as a function of time. Freundlich parameters 'kf' and 'm' are varied and
% output is written to a file of .csv format in the directory the script is ran.
%

clc(); % clear console

% MARK: PFM dimensions
rPFM = 1; % radius of PFM (cm)
zPFM = 1; % height of PFM (cm)
sMax = 10; % number of streamtubes in PFM
dy = rPFM/sMax; % width of stream tube (cm)

% MARK: sorbent matrix properties 
v = 2.00; % seepage velocity (cm/hr)
p = 0.63; % porosity of sorptive matrix
b = 1.00; % thickness of sorptive matrix (cm)
pb = 0.52*1000; % bulk density of sorptive matrix (mg/mL)

% MARK: time constraints
tFinal = 10; % total time of deployment (hr)
dt = .5; % time increment (hr)

% MARK: tracer parameters
cInit = 1000; % inital concentration (mg/L)
mInitPFM = 0; % variable declaration for initial mass in PFM (mg)

% MARK: file handling for script output
now = clock; % current time and date
filePath = pwd+"/"; % variable declaration for file name
for i=1:6 % add date and time to file name
    filePath = filePath+now(i)+ "-";
end
filePath = filePath+".csv"; % add file extension
file = fopen(filePath, 'w+'); % open file for writing
fprintf(file, 'kf, m, time, dMR\n'); % write header to file

% MARK: modeling start
for kf=1.20:0.25:1.50 % vary partitioning coefficient (mL/mg)
    for m=0.05:0.05:1.00 % vary Fruendlich exponent 
        % write data
        format = '%2.3f, %2.3f,';
        fprintf(file, format, kf, m);   
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
                     c = (n*v*t-n*x)/(pd*k*m*x).^(1/m-1)
                    mRtube = z*dy*n*integral(c,0,xd)+pb*kf*integral(c.^m,0,xd)+(p*cInit+pb*kf*(cInit).^m)*(xd-xb)
                    mRtube += mRtube
                else
                     % mass remaining in streamtube (mg)
                     mRtube = z*dy*n*integral(c,0,xd)+pb*kf*integral(c.^m,0,xd)
                     mRtube += mRtube
                end
                mInitPFM = mInitPFM + mInit; % initial mass in PFM (mg)
                mRPFM = mRPFM + mRtube; % mass remaining in PFM at time t (mg)
            end 
        dMR = (mRPFM/mInitPFM); % dimensionless mass remaining in PFM
        % write data
        format = '%5.5f, %3.5f\n';
        fprintf(file, format, t, dMR);
        end
        fprintf(file,'\n');
    end
end
fclose(file); % close file
