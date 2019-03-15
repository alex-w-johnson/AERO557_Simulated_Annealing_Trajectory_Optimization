function [depDatTim,arrDatTim,...
    R_depart,V_depart,...
    R_arrive,V_arrive,...
    tof,minDV,typeTraj,minDV1and2] = ...
    PatConTransferFunc(x,...
    rParkDep, rParkArr,...
    depBodyStruct, arrBodyStruct, mu_cent)
% PCTRANSFERFUNC Function to use in simulannealbnd to find optimal
% patched conics transfer window to and from two solar system bodies that
% orbit the Sun using a Lambert's problem solver. Ephemeris data for the
% bodies are calculated using the AERO557planetcoe_and_sv() function.
% 
% =========================================================================
% 
%   Inputs:
% 
%     x = solution input 2-element column vector, whose elements are:
%         depDate =(1){num} orbit transfer departure Julian Date
%         arrDate =(2){num} orbit transfer arrival Julian Date
%     
%     rParkDep ={num} circular orbit radius of departure parking orbit [km]
% 
%     rParkArr ={num} circular orbit radius of arrival parking orbit [km]
% 
%     depBodyStruct = struct containing departure body information:
%         .ID ={int} planet_id - planet identifier (per 
%           AERO557planetcoe_and_sv.m documentation):
%             1 = Mercury
%             2 = Venus
%             3 = Earth
%             4 = Mars
%             5 = Jupiter
%             7 = Uranus
%             8 = Neptune
%             9 = Pluto
%         .mu ={num} gravitational parameter of departure body [km^3/s^2]
% 
%     arrBodyStruct = struct containing departure body information:
%         .ID ={int} planet_id - planet identifier (per 
%           AERO557planetcoe_and_sv.m documentation):
%             1 = Mercury
%             2 = Venus
%             3 = Earth
%             4 = Mars
%             5 = Jupiter
%             7 = Uranus
%             8 = Neptune
%             9 = Pluto
%         .mu ={num} gravitational parameter of arrival body [km^3/s^2]
% 
%     mu_cent ={num} gravitational parameter of central body [km^3/s^2]
% 
% =========================================================================
% 
% Outputs:
% 
%     depDatTim ={datetime} (UT) date and time of departure on trajectory
% 
%     arrDatTim ={datetime} (UT) date and time of arrival on trajectory
% 
%     R_depart,V_depart ={num} vectors of heliocentric inertial 
%         position [km] and velocity [km/s] at departure
% 
%     R_arrive,V_arrive ={num} vectors of heliocentric inertial 
%         position [km] and velocity [km/s] at arrival
% 
%     tof ={num} time of flight of patched conics trajectory [days]
% 
%     minDV ={num} patched-conic delta-V [km/s] required to perform
%         departure and arrival burns
% 
%     typeTraj ={int} integer indicating type of interplanetary trajectory
%         used by Lambert solution for trajectory
% 
%     minDV1and2 ={num} 2 element row vector of delta-V [km/s] required
%         for: 
%             -departure burn (1) 
%             -arrival burn   (2)
% 
% =========================================================================
% 
% Notes regarding below code:
%     This function uses the Lambert-targeter for ballistic flights written
%     by Rody Oldenhuis, based on methodology described by Izzo, Lancaster,
%     Blanchard, and Gooding. See lambert.m documentation for relevant
%     references and usage information. This method is faster and more 
%     robust compared to other solvers, making it well-suited for
%     optimization.
% 
%     This function additionally uses a planetary ephemeris code that is 
%     far faster than the Aerospace Toolbox's planetEphemeris() code.
%     Refactoring of the code using planetEphemeris() is possible, but has
%     been found to significantly increase the computation time per
%     iteration, which may slow the annealing process unless computational
%     time can be saved by reducing the time of other calculations.

% check inputs structure
if ~(iscolumn(x) && length(x) == 2)
    error('Input x must be 2-element row vector!');
end
numArgs = 6;
if nargin < numArgs
    error('Not enough input arguments.')
elseif nargin > numArgs
    error('Too many input arguments.')
end

% Begin Calculating Planetary Ephemeris Data
depDate = x(1);
arrDate = x(2);

% Ephemeris outputs are in [km] for distance in the heliocentric inertial 
% frame, and [km/s] for velocity in the heliocentric inertial frame.

% Time formatting for ephemeris data input
depDatTim = datetime(depDate,'ConvertFrom','juliandate');
arrDatTim = datetime(arrDate,'ConvertFrom','juliandate');
[yDep,MDep,dDep] = ymd(depDatTim);
[hDep,mDep,sDep] = hms(depDatTim);
[yArr,MArr,dArr] = ymd(arrDatTim);
[hArr,mArr,sArr] = hms(arrDatTim);

% Ephemeris data calculation
[~,R_depart,V_depart,~] = AERO557planetcoe_and_sv(depBodyStruct.ID,...
    yDep,MDep,dDep,hDep,mDep,sDep);
[~,R_arrive,V_arrive,~] = AERO557planetcoe_and_sv(arrBodyStruct.ID,...
    yArr,MArr,dArr,hArr,mArr,sArr);
R_depart = R_depart';V_depart = V_depart';
R_arrive = R_arrive';V_arrive = V_arrive';

% Calculate Lambert trajectory
tof = abs(x(1)-x(2)); % Transfer Time of Flight in Days
[V1type1,V2type1,~,~] = lambert(R_depart,R_arrive,tof,0,mu_cent);  % Type 1 
[V1type2,V2type2,~,~] = lambert(R_depart,R_arrive,-tof,0,mu_cent); % Type 2

% Calculate V_infinity for both types of trajectory
Vinf1t1 = V1type1 - V_depart; % V_infinity at departure body, type 1 traj.
Vinf2t1 = V2type1 - V_arrive; % V_infinity at arrival body, type 1 traj.
Vinf1t2 = V1type2 - V_depart; % V_infinity at departure body, type 2 traj.
Vinf2t2 = V2type2 - V_arrive; % V_infinity at arrival body, type 2 traj.
v_inf1t1 = norm(Vinf1t1); % (scalar) v_infinity at departure body, type 1
v_inf2t1 = norm(Vinf2t1); % (scalar) v_infinity at arrival body, type 1
v_inf1t2 = norm(Vinf1t2); % (scalar) v_infinity at departure body, type 2
v_inf2t2 = norm(Vinf2t2); % (scalar) v_infinity at arrival body, type 2

% Calculate delta-v for departure and arrival for both types of
% trajectories

dV1t1 = calcBurnDv(depBodyStruct.mu,rParkDep,v_inf1t1); %dep. DV, type 1
dV2t1 = calcBurnDv(arrBodyStruct.mu,rParkArr,v_inf2t1); %arr. DV, type 1
dV1t2 = calcBurnDv(depBodyStruct.mu,rParkDep,v_inf1t2); %dep. DV, type 2
dV2t2 = calcBurnDv(depBodyStruct.mu,rParkDep,v_inf2t2); %arr. DV, type 2

% Calculate total delta-v for type 1 and type 2 trajectories
dVtype1 = [dV1t1,dV2t1];
dVtype2 = [dV1t2,dV2t2];
dVtot1 = sum(dVtype1);
dVtot2 = sum(dVtype2);

dVMat = [dVtype1;dVtype2];

[minDV,typeTraj] = min([dVtot1,dVtot2]);
minDV1and2 = dVMat(typeTraj,:);

% Function for calculating departure or arrival burn delta-V
    function [deltaV] = calcBurnDv(mu,rp,vinf)
        % Calculates 2-body patched conic delta-v for injection or
        % departure burn from/to a hyperbolic orbit to/from a circular
        % parking orbit
        
        vCirc = sqrt(mu/rp); % Velocity of circular parking orbit [km/s]
        % Curtis, Orbital Mechanics for Engineering Students 3rd Ed. 
        % Eq. 8.41
        
        vp = sqrt(vinf^2+2*mu/rp); % Velocity at periapsis of hyperbolic orbit [km/s]
        % Eq. 8.40
        
        deltaV = vp - vCirc; % delta-V required to get from circular
        % parking orbit to hyperbolic orbit or vice-versa
    end

end

