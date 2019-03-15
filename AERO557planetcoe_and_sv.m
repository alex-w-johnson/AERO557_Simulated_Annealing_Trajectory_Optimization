%function used to calc the coes and positions for each planet

% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
function [coe, r, v, jd] = AERO557planetcoe_and_sv ...
(planet_id, year, month, day, hour, minute, second)
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%{
This function calculates the orbital elements and the state
vector of a planet from the date (year, month, day)
and universal time (hour, minute, second).
mu - gravitational parameter of the sun (km^3/s^2)
deg - conversion factor between degrees and radians
pi - 3.1415926...
coe - vector of heliocentric orbital elements
[h e RAAN inc w TA a w_hat L M E],
where
h = angular momentum (km^2/s)
e = eccentricity
RA = right ascension (deg)
incl = inclination (deg)
w = argument of perihelion (deg)
TA = true anomaly (deg)
a = semimajor axis (km)
w_hat = longitude of perihelion ( = RA + w) (deg)
Appendix D Page 85 of 101 10/27/09 9:07 AM
L = mean longitude ( = w_hat + M) (deg)
M = mean anomaly (deg)
E = eccentric anomaly (deg)
planet_id - planet identifier:
1 = Mercury
2 = Venus
3 = Earth
4 = Mars
5 = Jupiter
7 = Uranus
8 = Neptune
9 = Pluto
year - range: 1901 - 2099
month - range: 1 - 12
day - range: 1 - 31
hour - range: 0 - 23
minute - range: 0 - 60
second - range: 0 - 60
j0 - Julian day number of the date at 0 hr UT
ut - universal time in fractions of a day
jd - julian day number of the date and time
J2000_coe - row vector of J2000 orbital elements from Table 9.1
rates - row vector of Julian centennial rates from Table 9.1
t0 - Julian centuries between J2000 and jd
elements - orbital elements at jd
r - heliocentric position vector
v - heliocentric velocity vector
User M-functions required: J0, kepler_E, sv_from_coe
User subfunctions required: planetary_elements, zero_to_360
%}
% --------------------------------------------------------------------
mu = 1.327124e11; %km3/s2
deg = pi/180;
%...Equation 5.48:
jd = JDcalc(year,month,day,hour,minute,second);

%...Equation 8.93a:
t0 = (jd - 2451545)/36525;

%...Equation 8.93b:
%...Obtain the data for the selected planet:
elements = AERO451planetary_elements2(planet_id, t0);
a = elements(1);
e = elements(2);
%...Equation 2.71:
h = sqrt(mu*a*(1 - e^2));
%...Reduce the angular elements to within the range 0 - 360 degrees:
inc = elements(3);
RAAN = zero_to_360(elements(4));
w_hat = zero_to_360(elements(5));
L = zero_to_360(elements(6));
w = zero_to_360(w_hat - RAAN);
M = zero_to_360((L - w_hat));
%...Algorithm 3.1 (for which M must be in radians)
% %use an initial estimate E
ecc=e;
M = M*pi/180;
if (M <pi)
     E = M+ecc/2;
else
     E=M-ecc/2;
end
% 
% %set an error tolerance
 tol = 1e-6;
 nmax = 50;
% 
% %Use Newton's to iterate
 ratio = 1;
 m=0;
 count =1;
 while (abs(ratio)>tol)&&(m<=nmax)
     m=m+1;
     [funcFE,funcFEdot]=FfuncE(E,M,ecc);
     ratio = funcFE/funcFEdot;
     E=E-ratio;
     count = count+1;
 end

%...Equation 3.13 (converting the result to degrees):
TA = zero_to_360...
(2*atan(sqrt((1 + e)/(1 - e))*tan(E/2))/deg);
coe = [h e RAAN inc w TA a w_hat L M E];
%...Algorithm 4.5 (for which all angles must be in radians):
[r,v]=...
    COEStoRV2sun(h,ecc,inc,RAAN, w,TA);

%% COEStoRV with sun
function [rvect,vvect]=...
    COEStoRV2sun(h,ecc,inc,raan, omega,theta)

%compute r and v from COES
mu = 1.327124e11; %km3/s2

theta = theta*pi/180;
omega = omega*pi/180;
raan = raan*pi/180;

%T = (2*pi/sqrt(mu))*a^(1.5);

%energy = -mu/(2*a);

%n = mu^(0.5)/a^(1.5);

%p = a*(1-ecc^2);

%h= (mu*p)^(1/2);

%E = 2*atan(sqrt((1-ecc/1+ecc))*tan(theta/2));

%book step
rvectx = (h^2/mu)*(1/(1+ecc*cos(theta))).*[cos(theta);sin(theta);0];
vvectx = (mu/h).*[-sin(theta);ecc+cos(theta);0];

%matrix conversion back into geocentric
term1 =[cos(omega) sin(omega) 0;...
        -sin(omega) cos(omega) 0; 0 0 1];
term2 =[1 0 0; 0 cosd(inc) sind(inc); 0 -sind(inc) cosd(inc)];
term3 =[cos(raan) sin(raan) 0;...
        -sin(raan) cos(raan) 0; 0 0 1];
    
convmat = term1*term2*term3;
%invconvmat = inv(convmat);

rvect = convmat\rvectx;
vvect = convmat\vvectx;
end


%% planetary elements

function [planet_coes] = AERO451planetary_elements2(planet_id,T)
% Planetary Ephemerides from Meeus (1991:202-204) and J2000.0
% Output:
% planet_coes
% a = semimajor axis (km)
% ecc = eccentricity
% inc = inclination (degrees)
% raan = right ascension of the ascending node (degrees)
% w_hat = longitude of perihelion (degrees)
% L = mean longitude (degrees)

% Inputs:
% planet_id - planet identifier:
% 1 = Mercury
% 2 = Venus
% 3 = Earth
% 4 = Mars
% 5 = Jupiter
% 6 = Saturn
% 7 = Uranus
% 8 = Neptune

if planet_id == 1
    a = 0.387098310; % AU but in km later
    ecc = 0.20563175 + 0.000020406*T - 0.0000000284*T^2 - 0.00000000017*T^3;
    inc = 7.004986 - 0.0059516*T + 0.00000081*T^2 + 0.000000041*T^3; %degs
    raan = 48.330893 - 0.1254229*T-0.00008833*T^2 - 0.000000196*T^3; %degs
    w_hat = 77.456119 +0.1588643*T -0.00001343*T^2+0.000000039*T^3; %degs
    L = 252.250906+149472.6746358*T-0.00000535*T^2+0.000000002*T^3; %degs
elseif planet_id == 2
    a = 0.723329820; % AU
    ecc = 0.00677188 - 0.000047766*T + 0.000000097*T^2 + 0.00000000044*T^3;
    inc = 3.394662 - 0.0008568*T - 0.00003244*T^2 + 0.000000010*T^3; %degs
    raan = 76.679920 - 0.2780080*T-0.00014256*T^2 - 0.000000198*T^3; %degs
    w_hat = 131.563707 +0.0048646*T -0.00138232*T^2-0.000005332*T^3; %degs
    L = 181.979801+58517.8156760*T+0.00000165*T^2-0.000000002*T^3; %degs
elseif planet_id == 3 
    a = 1.000001018; % AU
    ecc = 0.01670862 - 0.000042037*T - 0.0000001236*T^2 + 0.00000000004*T^3;
    inc = 0.0000000 + 0.0130546*T - 0.00000931*T^2 - 0.000000034*T^3; %degs
    raan = 0.0; %degs
    w_hat = 102.937348 + 0.3225557*T + 0.00015026*T^2 + 0.000000478*T^3; %degs
    L = 100.466449 + 35999.372851*T - 0.00000568*T^2 + 0.000000000*T^3; %degs
elseif planet_id == 4
    a = 1.523679342; % AU
    ecc = 0.09340062 + 0.000090483*T - 0.00000000806*T^2 - 0.00000000035*T^3;
    inc = 1.849726 - 0.0081479*T - 0.00002255*T^2 - 0.000000027*T^3; %degs
    raan = 49.558093 - 0.2949846*T-0.00063993*T^2 - 0.000002143*T^3; %degs
    w_hat = 336.060234 +0.4438898*T -0.00017321*T^2+0.000000300*T^3; %degs
    L = 355.433275+19140.2993313*T+0.00000261*T^2-0.000000003*T^3; %degs
elseif planet_id == 5
    a = 5.202603191 + 0.0000001913*T; % AU
    ecc = 0.04849485+0.000163244*T - 0.0000004719*T^2 + 0.00000000197*T^3;
    inc = 1.303270 - 0.0019872*T + 0.00003318*T^2 + 0.000000092*T^3; %degs
    raan = 100.464441 + 0.1766828*T+0.00090387*T^2 - 0.000007032*T^3; %degs
    w_hat = 14.331309 +0.2155525*T +0.00072252*T^2-0.000004590*T^3; %degs
    L = 34.351484+3034.9056746*T-0.00008501*T^2+0.000000004*T^3; %degs
elseif planet_id == 6
    a = 9.5549009596 - 0.0000021389*T; % AU
    ecc = 0.05550862 - 0.000346818*T -0.0000006456*T^2 + 0.00000000338*T^3;
    inc = 2.488878 + 0.0025515*T - 0.00004903*T^2 + 0.000000018*T^3; %degs
    raan = 113.665524 - 0.2566649*T-0.00018345*T^2 + 0.000000357*T^3; %degs
    w_hat = 93.056787 +0.5665496*T +0.00052809*T^2-0.000004882*T^3; %degs
    L = 50.077471+1222.1137943*T+0.00021004*T^2-0.000000019*T^3; %degs
elseif planet_id == 7
    a = 19.218446062-0.0000000372*T+0.00000000098*T^2; % AU
    ecc = 0.04629590 - 0.000027337*T + 0.0000000790*T^2 + 0.00000000025*T^3;
    inc = 0.773196 - 0.0016869*T + 0.00000349*T^2 + 0.00000000016*T^3; %degs
    raan = 74.005947 + 0.0741461*T+0.00040540*T^2 +0.000000104*T^3; %degs
    w_hat = 173.005159 +0.0893206*T -0.00009470*T^2+0.000000413*T^3; %degs
    L = 314.055005+428.4669983*T-0.00000486*T^2-0.000000006*T^3; %degs
elseif planet_id == 8
    a = 30.110386869-0.0000001663*T+0.00000000069*T^2; % AU
    ecc = 0.00898809 + 0.000006408*T -0.0000000008*T^2;
    inc = 1.769952 +0.0002557*T +0.00000023*T^2 -0.0000000000*T^3; %degs
    raan = 131.784057 - 0.0061651*T-0.00000219*T^2 - 0.000000078*T^3; %degs
    w_hat = 48.123691 +0.0291587*T +0.00007051*T^2-0.000000000*T^3; %degs
    L = 304.348665+218.4862002*T+0.00000059*T^2-0.000000002*T^3; %degs
end

planet_coes = [a;ecc;inc;raan;w_hat;L];
%Convert to km:
au = 149597870;
planet_coes(1) = planet_coes(1)*au;
end

%% JD calc

function JD1 = JDcalc(Y,M,D,UThr,UTmin,UTsec)


UT=UThr+UTmin/60+UTsec/60;


JD1 = 367*Y - floor((7*(Y+floor((M+9)/12)))/4) + ...
    floor((275*M)/9) + D + 1721013.5 + UT/24; 
end

%% FfuncE
function [funcF,funcFdot] =FfuncE(E,M,ecc)
    funcF = E -ecc*sin(E)-M;
    funcFdot = 1-ecc*cos(E);
end

%% zero to 360
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~
function y = zero_to_360(x)
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~
%{
This function reduces an angle to lie in the range 0 - 360 degrees.
x - the original angle in degrees
y - the angle reduced to the range 0 - 360 degrees
%}
% ---------------------------
if x >= 360
x = x - fix(x/360)*360;
elseif x < 0
x = x - (fix(x/360) - 1)*360;
end
y = x;
end %zero_to_360

end %main