clear all;
glvs;

ts = 1/100;

load cpt.mat; 

avp0 = getat(avpr,800);

ins = insinit(avp0(1:9), ts);

avperr = avperrset([200;300], 10, 100);

imuerr = imuerrset(100, 10000, 0.1, 105);

Pmin = [avperrset([0.2,1.0],0.01,0.2); gabias(0.01, [10,100]); [0.01;0.01;0.01]; 0.001].^2;

Rmin = vperrset(0.01, 0.03).^2;


%% Combined Error Model
     [avp0, xkpk1, zkrk1, sk1,ins,kf,count] = sinsgps_com_lq(imu(800/ts:1400/ts,:), gps, ins, avperr, imuerr, rep3(1), 0.1, vperrset(1,10), Pmin, Rmin, 'avped');
     [avp1, xkpk1, zkrk1, sk1] = sinsgps(imu(800/ts+count:1400/ts,:), gps, ins, avperr, imuerr, rep3(1), 0.1, vperrset(1,10), Pmin, Rmin, 'avped');
 %%    
avp1=[avp0;avp1];
 avperr=avpcmpplot(avpr, avp1);
save('cpt_tra','avperr');