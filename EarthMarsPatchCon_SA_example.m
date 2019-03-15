close all, clear all, clc

%% Patched conics optimization
% This code calculates the optimal patched-conic orbital transfer from
% Earth to Mars over a 4 year period of departure dates (2020-2024) and a 6
% year arrival window for Mars arrival dates (2021-2027).

% Assumptions:
% - No mid-course correction burn(s)
% - No orbital plane change necessary at departure or arrival
% - Burn maneuvers are impulsive
% - Patched conic two-body orbital dynamics
% - Negligible drag, SRP, n-body, or other perturbational effects

% Planetary bodies constants
mu_Earth = 398600; %[km^3/s^2]
soi_Earth = 925000; %[km]
mu_Mars = 42828; %[km^3/s^2]
soi_Mars = 577000; %[km]
mu_Sun = 132712000000; %[km^3/s^2]
r_Earth = 6378; %[km]
r_Mars = 3396; %[km]

% Circular Parking Orbit Characteristics
altPark_Earth = 400; %[km] altitude of circular parking orbit @ departure
altPark_Mars = 500; %[km] altitude of circular parking orbit @ arrival
rPark_Earth = r_Earth + altPark_Earth; %[km] orbital radius @ departure 
rPark_Mars = r_Mars + altPark_Mars; %[km] orbital radius @ arrival 

% Departure and Arrival times
depTimeJD = [2458849.5,2460310.5]; % Departure time Julian Dates
% 01-January-2020 00:00:00 UTC to 01-January-2024 00:00:00 UTC
arrTimeJD = [2459215.5,2461406.5]; % Arrival time Julian Dates
% 01-January-2021 00:00:00 UTC to 01-January-2027 00:00:00 UTC

% Body structure inputs for Patched Conics function
depBody.ID = 3; % Earth planet ID, per AERO557planetcoe_and_sv.m doc
depBody.mu = mu_Earth;
arrBody.ID = 4; % Mars planet ID, per AERO557planetcoe_and_sv.m doc
arrBody.mu = mu_Mars;

% Set Anonymous function with inputs for PatConTransferFunc.m
pctFunc = @(x)PatConTransferFunc(x,rPark_Earth,rPark_Mars,...
    depBody,arrBody,mu_Sun);

% Porkchop Plot generation
numTimePts = 100;
depJDVec = linspace(depTimeJD(1),depTimeJD(2),numTimePts);
arrJDVec = linspace(arrTimeJD(1),arrTimeJD(2),numTimePts);
dv = zeros([length(depJDVec),length(arrJDVec)]);
for i1 = 1:length(depJDVec)
    for j1 = 1:length(arrJDVec)
        dv(i1,j1) = pctFunc([depJDVec(i1);arrJDVec(j1)]);       
    end
end
%Plot Figure
figure('units','normalized','outerposition',[0 0 1 1])
contourf(depJDVec-depJDVec(1),arrJDVec-arrJDVec(1),dv',0:2:20,'k',...
    'showtext','on','linewidth',1);
grid on, grid minor, hold on;
xlabel('Departure Date: Days since 01-Jan-2020 00:00:00 (UT)')        
ylabel('Arrival Date: Days since 01-Jan-2021 00:00:00 (UT)')
title('Total Earth-Mars Mission \DeltaV [km/s]')
axPCPlot = gca;
axPCPlot.FontSize = 16;

% Set Initial Guess
minDV = min(min(dv));
[minRow,minCol] = find(minDV == dv);
x0 = [depJDVec(minRow);arrJDVec(minCol)];
plot(x0(1)-depJDVec(1),x0(2)-arrJDVec(1),'rs','markersize',10,...
    'markerfacecolor','r');
% % % % % text(x0(1)-depJDVec(1)+20,x0(2)-arrJDVec(1)+(1024-671),...
% % % % %     sprintf('Depart: %.1f \nArrive: %.1f\nTotal dV: %.2f',...
% % % % %     x0(1)-depJDVec(1),x0(2)-depJDVec(2),minDV),...
% % % % %     'Fontsize',16,'fontweight','bold','color','r','HorizontalAlignment',...
% % % % %     'center');
% initial guess at brute-force minimum delta-V

% Set simulated annealing solver options
saOpts = optimoptions(@simulannealbnd,...
    'MaxIterations',1000,...
    'FunctionTolerance',1e-8,...
    'Display','final',...
    'TemperatureFcn','temperaturefast',...
    'PlotFcn',{@saplotbestf,@saplottemperature,@saplotf,@saplotstopping});
lowBnd = [depTimeJD(1);arrTimeJD(1)]; %lower bound of solution values
uppBnd = [depTimeJD(2);arrTimeJD(2)]; %upper bound of solution values

% Solve for trajectory
% input('Press Enter to begin Simulated Annealing:')
[xOut,dvOut] = simulannealbnd(pctFunc,...
    x0,lowBnd,uppBnd,saOpts);
plot(axPCPlot,xOut(1)-depJDVec(1),xOut(2)-arrJDVec(1),'bo',...
    'markersize',10, 'linewidth',2);

% Calculate guess patched conic information
[depDatTim_g,arrDatTim_g,R_depart,V_depart,R_arrive,V_arrive,tof_g,~,...
    typeTraj_g,minDV1and2_g] = PatConTransferCalc(x0,...
    rPark_Earth,rPark_Mars,depBody,arrBody,mu_Sun);

% Display guess results
fprintf('\n')
disp('----- Simulated Annealing Results: -----')
disp('*** Original Brute Force Solution ***')
disp(['Departure Julian Date:   ',num2str(x0(1))])
disp(['Departure Date/Time (UT) ',sprintf('%s',depDatTim_g)])
disp(['Arrival Julian Date:     ',num2str(x0(2))])
disp(['Arrival Date/Time (UT)   ',sprintf('%s',arrDatTim_g)])
disp(['Time of Flight (days):   ',num2str(tof_g)])
disp(['Type of Trajectory:      ',num2str(typeTraj_g)])
disp(['Departure Delta-V [km/s]:',num2str(minDV1and2_g(1))])
disp(['Arrival Delta-V [km/s]:  ',num2str(minDV1and2_g(2))])
disp(['Total Delta-V Required:  ',num2str(sum(minDV1and2_g))])
fprintf('\n')

% Calculate optimized patched conic information
[depDatTim_o,arrDatTim_o,R_depart_o,V_depart_o,R_arrive_o,V_arrive_o,...
    tof_o,~,typeTraj_o,minDV1and2_o] = PatConTransferCalc(xOut,...
    rPark_Earth,rPark_Mars,depBody,arrBody,mu_Sun);
percSaved = abs(sum(minDV1and2_g)-sum(minDV1and2_o))/sum(minDV1and2_g);
disp('*** Optimized Simulated Annealing Solution ***')
disp(['Departure Julian Date:   ',num2str(xOut(1))])
disp(['Departure Date/Time (UT) ',sprintf('%s',depDatTim_o)])
disp(['Arrival Julian Date:     ',num2str(xOut(2))])
disp(['Arrival Date/Time (UT)   ',sprintf('%s',arrDatTim_o)])
disp(['Time of Flight (days):   ',num2str(tof_o)])
disp(['Type of Trajectory:      ',num2str(typeTraj_o)])
disp(['Departure Delta-V [km/s]:',num2str(minDV1and2_o(1))])
disp(['Arrival Delta-V [km/s]:  ',num2str(minDV1and2_o(2))])
disp(['Total Delta-V Required:  ',num2str(sum(minDV1and2_o))])
disp(['% of Guess Delta-V Saved:',num2str(percSaved)])
fprintf('\n')