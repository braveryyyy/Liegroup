function [att0, attk, xkpk,ins] = alignvn_lq(imu, qnb, pos, phi0, imuerr, wvn,ins, isfig)
% SINS initial align uses Kalman filter with vn as measurement.
% Kalman filter states: 
%    [phiE,phiN,phiU, dvE,dvN,dvU, ebx,eby,ebz, dbx,dby,dbz]'.
%
% Prototype: [att0, attk, xkpk] = alignvn(imu, qnb, pos, phi0, imuerr, wvn, isfig)
% Inputs: imu - IMU data
%         qnb - coarse attitude quaternion
%         pos - position
%         phi0 - initial misalignment angles estimation
%         imuerr - IMU error setting
%         wvn - velocity measurement noise (3x1 vector)
%         isfig - figure flag
% Outputs: att0 - attitude align result
%         attk, xkpk - for debug
%
% Example:
%	avp0 = [[0;0;0], zeros(3,1), glv.pos0];
%	imuerr = imuerrset(0.03, 100, 0.001, 10);
%	imu = imustatic(avp0, 1, 300, imuerr);
%	phi = [.5; .5; 5]*glv.deg;
%	wvn = [0.01; 0.01; 0.01];
%	[att0, attk, xkpk] = alignvn(imu, avp0(1:3)', avp0(7:9)', phi, imuerr, wvn);
%



global glv
    if nargin<4,  phi0 = [1.5; 1.5; 3]*glv.deg;  end
    if nargin<5,  imuerr = imuerrset(0.01, 100, 0.001, 1);  end
    if nargin<6,  wvn = 0.01;  end;  if length(wvn)==1, wvn=repmat(wvn,3,1); end
    if nargin<7,  isfig = 1; end
    if length(qnb)==3, qnb=a2qua(qnb); end  %if input qnb is Eular angles.
    [nn, ts, nts] = nnts(2, diff(imu(1:2,end)));
    len = fix(length(imu)/nn)*nn;
    eth = earth(pos); vn = zeros(3,1); Cnn = rv2m(-eth.wnie*nts/2);
    kf = avnkfinit(nts, pos, phi0, imuerr, wvn);
    kf.lqmode='right';   %%  "kf.lqmode='right'" denotes right-invariant error model, and "kf.lqmode='left'" denotes left-invariant error model.
    [attk, xkpk] = prealloc(fix(len/nn), 7, 2*kf.n+1);
    ki = timebar(nn, len, 'Initial align using vn as meas.');
    Cnb = q2mat(qnb);
    for k=1:nn:len-nn+1
        wvm = imu(k:k+nn-1,1:6);  t = imu(k+nn-1,end);
        ins.wvm=wvm;
        ins = insupdate(ins, wvm);
        [phim, dvbm] = cnscl(wvm);
        dvn = Cnn*Cnb*dvbm;
        vn = vn + dvn + eth.gn*nts;
        Cnb=rv2m(-eth.wnin*nts)*Cnb*rv2m(phim);
        ins.Cnb=Cnb;
        if strcmp(kf.lqmode,'right')==1
         [kf.Phikk_1,kf.Gammak]=etm_static_right(ins);
        else
        [kf.Phikk_1,kf.Gammak]=etm_static_left(ins);
        end
        if strcmp(kf.lqmode,'right')==1
             kf.Hk=[zeros(3,3),-eye(3),zeros(3,6)];
        else
             kf.Hk=[zeros(3,3),-ins.Cnb,zeros(3,6)]; 
        end
        kf = kfupdate(kf, vn);
        
        if strcmp(kf.lqmode,'right')==1
               tem0=kf.xk(1:3);
                Q=askew(tem0);
                delat_R=eye(3)+sin(norm(tem0))/norm(tem0)*Q+(1-cos(norm(tem0)))/(norm(tem0))^2*Q^2;
                Cnb=delat_R*Cnb;  
                vn=vn-cross(vn,kf.xk(1:3))+kf.xk(4:6);
        else
        tem0=0.91*kf.xk(1:3);
        Q=askew(tem0);
        delat_R=eye(3)+sin(norm(tem0))/norm(tem0)*Q+(1-cos(norm(tem0)))/(norm(tem0))^2*Q^2;
        tem_Cnb=Cnb;
        Cnb=Cnb*delat_R;
        ins.Cnb=Cnb;
        vn=vn+tem_Cnb*0.91*kf.xk(4:6)
        end
                
         ins.vn=vn;
         kf.xk(1:3) = 0.09*kf.xk(1:3);
         
         kf.xk(4:6) = 0.09*kf.xk(4:6);
        attk(ki,:) = [m2att(Cnb); vn; t]';
        xkpk(ki,:) = [kf.xk; diag(kf.Pxk); t];
        ki = timebar;
%         if k>=15000
%             break;
%         end
    end
    ins.att=m2att(Cnb);
    attk(ki:end,:) = []; xkpk(ki:end,:) = [];
    att0 = attk(end,1:3)';
    resdisp('Initial align attitudes (arcdeg)', att0/glv.deg);
    if isfig, avnplot(attk, xkpk); end
    
function kf = avnkfinit(nts, pos, phi0, imuerr, wvn)
    eth = earth(pos); wnie = eth.wnie;
    kf = []; kf.s = 1; kf.nts = nts;
	kf.Qk = diag([imuerr.web; imuerr.wdb; zeros(6,1)])^2*nts;
    kf.Gammak = 1;
	kf.Rk = diag(wvn)^2/nts;
	kf.Pxk = diag([phi0; [1;1;1]; imuerr.eb; imuerr.db])^2;
	Ft = zeros(12); Ft(1:3,1:3) = askew(-wnie); kf.Phikk_1 = eye(12)+Ft*nts;
	kf.Hk = [zeros(3),eye(3),zeros(3,6)];
%     [kf.m, kf.n] = size(kf.Hk);
%     kf.I = eye(kf.n);
%     kf.xk = zeros(kf.n, 1);
%     kf.adaptive = 0;
%     kf.xconstrain = 0; kf.pconstrain = 0;
%     kf.fading = 1;
    kf = kfinit0(kf, nts);

function avnplot(attk, xkpk)
global glv
    if glv.isfig==0, return; end
    t = attk(:,end);
    myfigure;
	subplot(421); plot(t, attk(:,1:2)/glv.deg); xygo('pr'); title('Xi');
	subplot(423); plot(t, attk(:,3)/glv.deg); xygo('y');
	subplot(425); plot(t, xkpk(:,7:9)/glv.dph); xygo('eb'); 
	subplot(427); plot(t, xkpk(:,10:12)/glv.ug); xygo('db'); 
	subplot(422); plot(t, sqrt(xkpk(:,13:15))/glv.min); xygo('phi'); title('\surdPii')
	subplot(424); plot(t, sqrt(xkpk(:,16:18))); xygo('dV');
	subplot(426); plot(t, sqrt(xkpk(:,19:21))/glv.dph); xygo('eb');
 	subplot(428); plot(t, sqrt(xkpk(:,22:24))/glv.ug); xygo('db');   
