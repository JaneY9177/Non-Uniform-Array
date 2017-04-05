%%
%  Random Spacing Linear Array Side Lobe Optimization
%
%  Version 1
%  Zhengyu Peng
%%

clear;

%% Configuration
wavelength=1;
k=2*pi/wavelength;   % wave number
angleStep=0.05;
theta=0:angleStep:360;
elementNumber=16;
spacingMIN=0.5; % minimum space between two elements
spacingMAX=3.5;   % maximum space between two elements
mainbeam=5;     % width of the main lobe
thetaM=90;      % location of the main lobe

%% random spacing
% d0, d1, d2, d3, d4, d5, d6, d7, d8, d9, d10, d11, d12, d13, d14

% d=spacingMIN + (spacingMAX-spacingMIN).*rand(1,elementNumber-1);
% d=[0,d];

%% random symmetrical spacing
% d0, d1, d2, d3, d4, d5, d6, d7, d6, d5, d4, d3, d2, d1, d0

% d=[spacingMIN/2 + (spacingMAX/2-spacingMIN/2).*rand(1), spacingMIN + (spacingMAX-spacingMIN).*rand(1,elementNumber/2-1)];
% d1=fliplr(d);
% d=[0,d1(1:length(d1)-1),d1(length(d1))+d(1),d(2:length(d))];

%% random spacing + shift
% d0, d1, d2, d3, d4, d5, d6, d7, d0, d1, d2, d3, d4, d5, d6

d=spacingMIN + (spacingMAX-spacingMIN).*rand(1,elementNumber/2-1);
d=[0,d,spacingMIN*2 + (spacingMAX*2-spacingMIN*2).*rand(1),d];

%%
for nn=2:length(d)
    d(nn)=d(nn-1)+d(nn);
end

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
cvx_begin %quiet
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
