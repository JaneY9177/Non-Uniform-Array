%clear;

%% Configuration
wavelength=1;
k=2*pi/wavelength;   % wave number
angleStep=0.05;
theta=0:angleStep:360;
elementNumber=16;
spacingMIN=0.5;
spacingMAX=5;
mainbeam=7;
thetaM=90;

%% random search
% d=[spacingMIN/2 + (spacingMAX/2-spacingMIN/2).*rand(1), spacingMIN + (spacingMAX-spacingMIN).*rand(1,elementNumber/2-1)];
% %d=[spacingMIN/2,spacingMIN*ones(1,elementNumber/2-1)];
% 
% for nn=2:length(d)
%     d(nn)=d(nn-1)+d(nn);
% end

%% check results
load('resultd.mat');
load('resultw.mat')
[r,c]=size(resultd);
nn=21;
d=resultd(nn,:);
w=resultw(nn,:)';

%% Array factor
A=zeros(length(theta),elementNumber);

for nn=1:length(d)
    A(:,nn)=exp(-1i*k*d(length(d)-nn+1)*cosd(theta));
end

for nn=1:length(d)
    A(:,nn+length(d))=exp(1i*k*d(nn)*cosd(theta));
end

%% Main lobe
A_M=zeros(1,elementNumber);
for nn=1:length(d)
    A_M(:,nn)=exp(-1i*k*d(length(d)-nn+1)*cosd(thetaM));
end

for nn=1:length(d)
    A_M(:,nn+length(d))=exp(1i*k*d(nn)*cosd(thetaM));
end

%% Side lobe
theta_SL=[0:angleStep:thetaM-mainbeam/2,thetaM+mainbeam/2:angleStep:180];
%theta_ML=90+mainbeam/2:angleStep:180;

A_SL=zeros(length(theta_SL),elementNumber);
for nn=1:length(d)
    A_SL(:,nn)=exp(-1i*k*d(length(d)-nn+1)*cosd(theta_SL));
end

for nn=1:length(d)
    A_SL(:,nn+length(d))=exp(1i*k*d(nn)*cosd(theta_SL));
end

%% Optimization
cvx_begin
variable w(elementNumber) complex
minimize( max(abs(A_SL*w)) )
subject to
A_M*w==1;
cvx_end

%% Plot result
plot(theta,20*log10(abs(A*w))-max(20*log10(abs(A*w))));
axis([0,180,-30,0]);
