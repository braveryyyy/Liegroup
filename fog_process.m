clear all;
ggpsvars

psinstypedef('test_SINS_GPS_tightly_def');

load foginsgps_tc.mat; % vars: imu ephs obss t0

[nn, ts, nts] = nnts(1, diff(imu(1:2,end))); % imuplot(imu);

vp = getgnssvp(ephs, obss, obss(1,1), 1);

att = aligni0(datacut(imu,0,200), vp(4:6));

avp0 = [att; vp]; % avp=inspure(imu, [att;vp]);
%  davp = avperrset([5/0.0175;5/0.0175;60/0.0175], [0.01;0.01;0.01]*0, [1;1;1]*0);  %% large misalignment
 davp = avperrset([1/0.0175;1/0.0175;1/0.0175], [0.01;0.01;0.01]*0, [1;1;1]*0);%% small misalignment
   avp0 = avpadderr(avp0, davp);
phi = [5;5;60]*glv.deg;
imuerr1 = imuerrset(1, 1000, 0.1, [10,10,50]);
wvn = [0.01;0.01;0.01];
ins = insinit(avp0, ts);
%      [att0v, attkv0, xkpkv] = alignvn(datacut(imu,0,200),  avp0(1:3), ins.pos, phi, imuerr1, wvn,1);%% Traditional Error Model
   [att0v, attkv1, xkpkv] = alignvn_lq(datacut(imu,0,200), avp0(1:3), ins.pos, phi, imuerr1, wvn,ins,1);%% Lie Group Error Model
