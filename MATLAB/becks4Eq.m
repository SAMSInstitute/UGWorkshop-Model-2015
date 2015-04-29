function xp = becks4Eq(t,x)

xp = zeros(4,1);

global N0
global D

% growth rates (1/sec)
% values below are in 1/day; divide by 24*3600 to get seconds
munr = 12.0 ;  % original
% munr = 12.157197177976418 ;  % original
% munr = 6 ;
munc = 6.0 ;   % original
% munc = 2.5009944119115524 ;
% mupr = 1.6919353004300042 ;
% mupc = 0.13790519501669446 ;
mupr = 2.2 ;  % original
mupc = 2.2/8.0 ;  % original
% mupc = 2.2 ;
%
% mass (g)
mr = 1.6e-12 ;
mc = 8.2e-12 ;
mp = 9.8e-10 ;
%
% half-saturation constants (g/cc)
knr = 8.0e-6 ;
knc = 8.0e-6 ;
kpr = 1.0e-6 ;
kpc = 1.0e-6 ;
%
% death rates (1/sec)
% values below are in 1/day; divide by 24*3600 to get seconds
dr = 0.5/10; % original
% dr = 0.03535901442365246;
dc = 0.25/100; % original
% dc = 0.002636146332980288
% dr = 0.5 ;
% dc = 0.25/10 ;
dp = 0.08/100 ; % original
% dp = 4.96579147003975e-4;
%
% initial nutrient (g/cc)
% global d;
% N0 = 2.3e-4 ; % original
% N0 = 2.3e-4 ;
% N0 = d;
%
% yield coefficients
ypr = 0.12 ;
ypc = 0.12 ;
ynr = 0.1 ;
ync = 0.1 ;
%
% dilution rate
% values below are in 1/day; divide by 24*3600 to get seconds
% D = 0.14512023872492155;
% D = 0.14660367473913746;
% D = 0.1451;  % 
% D = 0.6 ;  % original
%
%
% Differential equations

xp(1,1) = x(1,1)*(munr*x(4,1)/(x(4,1)+knr) - dr) - mupr/ypr*mp/mr*x(1,1)/(kpr/mr+x(1,1))*x(3,1) - D*x(1,1) ;

xp(2,1) = x(2,1)*(munc*x(4,1)/(x(4,1)+knc) - dc) - mupc/ypc*mp/mc*x(2,1)/(kpc/mc+x(2,1))*x(3,1) - D*x(2,1) ;

xp(3,1) = x(3,1)*(mupr*x(1,1)/(kpr/mr+x(1,1)) + mupc*x(2,1)/(kpc/mc+x(2,1)) - dp) - D*x(3,1) ;

xp(4,1) = D*N0 - x(1,1)*munr*mr/ynr*x(4,1)/(knr+x(4,1)) - x(2,1)*munc*mc/ync*x(4,1)/(knc+x(4,1)) - D*x(4,1);


