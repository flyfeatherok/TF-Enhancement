%% Non-stationary signal WVD and TF image test data generator
function Test_Signal_Gen()
%% Parameter Configuration
Show=0;%Show the WVD and TF image in generation process? 0-no,1-yes. Notice: 'yes' will make the process slow
SignalMode=1;%Generator Type: 1-FH,2-LFM,3-SFM
Components=1;%Singal component number Components
SNR=2;%SNR
T=0.4;%Singal duration
T_t=0.04;%Frequency modulation rate cycle or frequency hopping cycle
FFTn=256;%WVD size
dT=0.0001;%Sampling time
r=1/(2*FFTn*dT);%resolution
T_n=T/T_t;%Time period
time=ones(1,T_n)*T_t;
FF=[40,70;
    80,30];%The frequency modulation variation range of SFM signals.
FA=[20*r,50*r;
    40*r,80*r];%The amplitude variation range of frequency modulation of SFM signals.
%i.e., the No.1 SFM signal start from phase j*2*pi*FA(1,1)*sin(2*pi*FF(1,1)*t), end by phase j*2*pi*FA(1,2)*sin(2*pi*FF(1,2)*t),
%  and the N0.2 SFM signal start from phase j*2*pi*FA(2,1)*sin(2*pi*FF(2,1)*t), end by phase j*2*pi*FA(2,2)*sin(2*pi*FF(2,2)*t).
%Notice: FA should not larger than FFTn/2
path='E:\TF_Test\';
mkdir(path);
SignalType=[{'FH signal'},{'LFM siganl'},{'SFM siganl'}];
%% Display
disp('Signal Generator Parameter Configuration:')
disp(['--------Singal Type: ', SignalType{SignalMode}])
disp(['--------Singal components: ', num2str(Components)])
if SignalMode==1
    disp(['--------Frequency hopping cycle: ', num2str(T_t), ' second'])
end
if SignalMode==2
    disp(['--------Frequency modulation rate cycle: ', num2str(T_t), ' second'])
end
if SignalMode==3
    disp('--------Frequency modulation variation range: ')
    for i=1:Components
        disp(['----------------Component ',num2str(Components),': ',num2str(FF(i,1)),'-',num2str(FF(i,2)),' Hz'])
    end
    disp('--------Amplitude variation range of frequency modulation: ')
    for i=1:Components
        disp(['----------------Component ',num2str(Components),': ',num2str(FA(i,1)),'-',num2str(FA(i,2)),' Hz'])
    end
end
disp(['--------Singal duration: ', num2str(T), ' second'])
disp(['--------Singal SNR: ', num2str(SNR), 'dB'])
disp(['--------Sampling rate: ', num2str(1/dT), 'Hz'])
disp(['--------WVD size: ', num2str(FFTn)])
disp(['--------Singal save path: ', path])
disp('........Generating........')

%% FH signal
if SignalMode==1
    tic
    TFSavePath=[path 'FH_Test.dat'];
    fid = fopen(TFSavePath,'wb');
    fd=rand(Components,T_n)*(FFTn-1)*r;
    [TF,IMG] = FP_TF_IMG(fd,SNR,FFTn,dT,r,time,Components,Show);
    fwrite(fid,TF,'float');
    fclose(fid);
    imwrite(IMG,[path 'FH_Test.bmp'],'bmp');
    figure(2);imshow(IMG);title('TF Image')
    toc
end

%% LFM signal
if SignalMode==2
    tic
    TFSavePath=[path 'LFM_Test.dat'];
    fid = fopen(TFSavePath,'wb');
    jerk=zeros(Components,T_n);
    fd=rand(Components,T_n)*(FFTn-1)*r;
    a=[zeros(Components,1),(fd(:,2:end)-fd(:,1:end-1))/T_t];
    [TF,IMG] = TF_WVD_LFM(jerk,a,fd,SNR,FFTn,dT,r,time,Components,Show);
    fwrite(fid,TF,'float');
    fclose(fid);
    imwrite(IMG,[path 'LFM_Test.bmp'],'bmp');
    figure(2);imshow(IMG); title('TF Image')
    toc
end

%% SFM signal
if SignalMode==3
    tic
    TFSavePath=[path 'SFM_Test.dat'];
    fid = fopen(TFSavePath,'wb');
    Dead_Zone=max(FA')/r;
    for i=1:Components
        fd(i,:)=Dead_Zone(i)*r+rand(1,T_n).*(FFTn-1-2*Dead_Zone(i))*r;
    end
    a=[zeros(Components,1),(fd(:,2:end)-fd(:,1:end-1))/T_t];
    [TF,IMG] = SFM_TF_IMG(a,fd,FF,FA,SNR,FFTn,dT,r,time,Components,Show);
    fwrite(fid,TF,'float');
    fclose(fid);
    imwrite(IMG,[path 'SFM_Test.bmp'],'bmp');
    figure(2);imshow(IMG);title('TF Image')
    toc
end

disp('--------DONE!---------')
end

%% FH TF
function [TF,IMG] = FP_TF_IMG(fd,SNR,FFTn,dT,r,time,signalnum,Show)
fs=1/dT;fc=0;
ph0=0;%pi/4;
T=sum(time);
t=1/fs:1/fs:T;
noise=1/(10^(SNR/20));
[sc,ss,IF(1,:)]=Hopping_carrier_signal(ph0,fs,fc,fd(1,:),time);
for i=2:signalnum
    [scx,ssx,IF(i,:)]=Hopping_carrier_signal(ph0,fs,fc,fd(i,:),time);
    sc=sc+scx;
    ss=ss+ssx;
end
sc=sc+randn(1,length(t))*noise;
ss=ss+randn(1,length(t))*noise;
s=complex(sc,ss);
%% TF image
T=sum(time);
N=round(IF/r)+1;
IMG=uint8(zeros(FFTn,int32(T/dT)));
for i=1:signalnum
    for w=1:T/dT
        IMG(N(i,w),w)=IMG(N(i,w),w)+255;
    end
end
%% WVD
TF=zeros(FFTn,1);
for i=1:int32(T/dT)-FFTn
    sx=s(i:i+FFTn-1);
    tfr=WVD(sx);
    Mtfr=max(max(abs(tfr)));
    TFR=tfr/Mtfr;
    TF=[TF TFR'];
    if Show==1
        figure(1);subplot(121);imshow(255-uint8(255*TFR));title('WVD');subplot(122);imshow(255-IMG(:,i:i+FFTn-1));title('TF Image');pause(0.001);
    end
end
TF=TF(:,2:end);
end

%% LFM TF
function [TF,IMG] = TF_WVD_LFM(jerk,a,fd,SNR,FFTn,dT,r,time,signalnum,Show)
fs=1/dT;fc=0;
h=1;ph0=0;%pi/4;
t=1/fs:1/fs:sum(time);
noise=1/(10^(SNR/20));
[sc,ss,IF(1,:),~]=carrier_signal(a(1,:),jerk(1,:),ph0,fs,fc,fd(1,:),time,h);
for i=2:signalnum
    [scx,ssx,IF(i,:),~]=carrier_signal(a(i,:),jerk(i,:),ph0,fs,fc,fd(i,:),time,h);
    sc=sc+scx;
    ss=ss+ssx;
end
sc=sc+randn(1,length(t))*noise;
ss=ss+randn(1,length(t))*noise;
s=complex(sc,ss);
%% TF image
T=sum(time);
N=round(IF/r)+1;
IMG=uint8(zeros(FFTn,int32(T/dT)));
for i=1:signalnum
    for w=1:T/dT
        IMG(N(i,w),w)=IMG(N(i,w),w)+255;
    end
end
%% WVD
TF=zeros(FFTn,1);
for i=1:int32(T/dT)-FFTn
    sx=s(i:i+FFTn-1);
    tfr=WVD(sx);
    Mtfr=max(max(abs(tfr)));
    TFR=tfr/Mtfr;
    TF=[TF TFR'];
    if Show==1
        figure(1);subplot(121);imshow(255-uint8(255*TFR));title('WVD');subplot(122);imshow(255-IMG(:,i:i+FFTn-1));title('TF Image');pause(0.001);
    end
end
TF=TF(:,2:end);
end

%% SFM TF
function [TF,IMG] = SFM_TF_IMG(a,fd,FF,FA,SNR,FFTn,dT,r,time,signalnum,Show)
fs=1/dT;fc=0;
ph0=0;%pi/4;
t=1/fs:1/fs:sum(time);
noise=1/(10^(SNR/20));
[s,IF(1,:)]=SFM_carrier_signal(ph0,a(1,:),fs,fc,fd(1,:),FF(1,:),FA(1,:),time);
for i=2:signalnum
    [sx,IF(i,:)]=SFM_carrier_signal(ph0,a(i,:),fs,fc,fd(i,:),FF(i,:),FA(i,:),time);
    s=s+sx;
end
s=s+complex(randn(1,length(t))*noise,randn(1,length(t))*noise);
%% TF image
T=sum(time);
N=round(IF/r)+1;
IMG=uint8(zeros(FFTn,int32(T/dT)));
for i=1:signalnum
    for w=1:T/dT
        IMG(N(i,w),w)=IMG(N(i,w),w)+255;
    end
end
%% WVD
TF=zeros(FFTn,1);
for i=1:int32(T/dT)-FFTn
    sx=s(i:i+FFTn-1);
    tfr=WVD(sx);
    Mtfr=max(max(abs(tfr)));
    TFR=tfr/Mtfr;
    TF=[TF TFR'];
    if Show==1
        figure(1);subplot(121);imshow(255-uint8(255*TFR));title('WVD');subplot(122);imshow(255-IMG(:,i:i+FFTn-1));title('TF Image');pause(0.001);
    end
end
TF=TF(:,2:end);
end

%% Hopping carrier signal
function [sc,ss,dop]=Hopping_carrier_signal(ph0,fs,fc,fd,time)
T=length(time);
t=0:1/fs:time(1);
t=t(2:end);
Ph=2*pi*fd(1)*t+ph0;
sc=cos(2*pi*fc*t+Ph);
ss=sin(2*pi*fc*t+Ph);
dop=ones(1,length(t))*fd(1);
for k=2:T
    t=0:1/fs:time(k);
    t=t(2:end);
    phx=Ph(end);
    dopx=ones(1,length(t)).*fd(k);
    scx=cos(2*pi*(fc+fd(k))*t+phx);
    ssx=sin(2*pi*(fc+fd(k))*t+phx);
    sc=[sc scx];
    ss=[ss ssx];
    dop=[dop dopx];
end
end

%% Sinusoidal Frequency Modulation carrier signal 
function [s,dop]=SFM_carrier_signal(ph0,a,fs,fc,fd,FF,FA,time)
%% Pure SFM
T=sum(time);
t=0:1/fs:T;
t=t(2:end);
FFa=(FF(end)-FF(1))/T;
FAs=FA(1)+(FA(end)-FA(1))/length(t):(FA(end)-FA(1))/length(t):FA(end);
dopS=FAs.*cos(2*pi*FF(1)*t+pi*FFa*t.^2);
PhS=FAs./(FF(1)+FFa*t).*sin(2*pi*FF(1)*t+pi*FFa*t.^2)+ph0;
scS=cos(PhS);
ssS=sin(PhS);
j1=complex(scS,ssS);
%% LFM
T=length(time);
t=0:1/fs:time(1);
t=t(2:end);
Dop=fd(1)+a(1)*t;
Ph=2*pi*fd(1)*t+2*pi*( a(1)*t.^2/2 )+ph0;
sc=cos(2*pi*fc*t+Ph);
ss=sin(2*pi*fc*t+Ph);
for k=2:T
    t=0:1/fs:time(k);
    t=t(2:end);
    Dopx=a(k)*t+Dop(end);
    phx=2*pi*( a(k)*t.^2/2 )+Ph(end);
    scx=cos(2*pi*(fc+Dop(end))*t+phx);
    ssx=sin(2*pi*(fc+Dop(end))*t+phx);
    sc=[sc scx];
    ss=[ss ssx];
    Dop=[Dop Dopx];
end
j2=complex(sc,ss);
%% 
dop=dopS+Dop;
s=j1.*j2;
end

%% carrier signal 
function [sc,ss,Dop,Ph,ar]=carrier_signal(a,jerk,ph0,fs,fc,fd0,time,h)
T=length(time);
t=0:1/fs:time(1);
t=t(2:end);
ar=a(1)+jerk(1)*t;
Dop=fd0(1)+h*a(1)*t+0.5*h*jerk(1)*t.^2;
Ph=2*pi*fd0(1)*t+2*pi*( h*a(1)*t.^2/2 )+2*pi*(h*jerk(1)*t.^3/6)+ph0;
sc=cos(2*pi*fc*t+Ph);
ss=sin(2*pi*fc*t+Ph);
for k=2:T
    t=0:1/fs:time(k);
    t=t(2:end);
    arx=a(k)+jerk(k)*t;
    Dopx=h*a(k)*t+0.5*h*jerk(k)*t.^2+Dop(end);
    phx=2*pi*( h*a(k)*t.^2/2 )+2*pi*(h*jerk(k)*t.^3/6)+Ph(end);
    scx=cos(2*pi*(fc+Dop(end))*t+phx);
    ssx=sin(2*pi*(fc+Dop(end))*t+phx);
    sc=[sc scx];
    ss=[ss ssx];
    ar=[ar arx];
    Ph=[Ph phx];
    Dop=[Dop Dopx];
end
end

%% WVD
function tfr=WVD(S)
L=length(S);
for i=1:L
    taumax=min([i-1,L-i,round(L/2)-1]);
    tau=-taumax:taumax;
    indices= rem(L+tau,L)+1;
    tfr(indices,i) = S(i+tau) .* conj(S(i-tau));
end; 
tfr= fft(tfr); 
tfr= real(tfr);
end