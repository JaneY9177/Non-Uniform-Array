%% random spacing linear array side lobe Optimization

%clear;

%% Configuration
wavelength=1;
k=2*pi/wavelength;   % wave number
angleStep=0.05;
theta=0:angleStep:360;
elementNumber=16;
spacingMIN=0.5; % minimum space between two elements
spacingMAX=2;   % maximum space between two elements
mainbeam=5;     % width of the main lobe
thetaM=90;      % location of the main lobe

%% random spacing
%     d=spacingMIN + (spacingMAX-spacingMIN).*rand(1,elementNumber-1);
%     d=[0,d];

d=spacingMIN + (spacingMAX-spacingMIN).*rand(1,elementNumber/2-1);
d=[0,d,spacingMIN*2 + (spacingMAX*2-spacingMIN*2).*rand(1),d];

for nn=2:length(d)
    d(nn)=d(nn-1)+d(nn);
end

%% check results
% load('resultd.mat');
% load('resultw.mat')
% [r,c]=size(resultd);
% nn=21;
% d=resultd(nn,:);
%w=resultw(nn,:)';

figure(1);
plot(d,zeros(1,length(d)),'x');
axis([0,d(length(d)),-1,1]);

%% Array factor
A=zeros(length(theta),elementNumber);
for nn=1:length(d)
    A(:,nn)=exp(1i*k*d(nn)*cosd(theta));
end

%% Main lobe
A_M=zeros(1,elementNumber);
for nn=1:length(d)
    A_M(:,nn)=exp(1i*k*d(nn)*cosd(thetaM));
end

%% Side lobe
theta_SL=[0:angleStep:thetaM-mainbeam/2,thetaM+mainbeam/2:angleStep:180];
%theta_ML=90+mainbeam/2:angleStep:180;

A_SL=zeros(length(theta_SL),elementNumber);
for nn=1:length(d)
    A_SL(:,nn)=exp(1i*k*d(nn)*cosd(theta_SL));
end

%% Optimization
cvx_begin quiet
variable w(elementNumber) complex
minimize( max(abs(A_SL*w)) )
subject to
A_M*w==1;
cvx_end

%% Plot result
%w=ones(1,elementNumber)';
figure(2);
plot(theta,20*log10(abs(A*w))-max(20*log10(abs(A*w))));
axis([0,180,-30,0]);
