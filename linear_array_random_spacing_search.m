%clear;

%% Configuration
wavelength=1;
k=2*pi/wavelength;   % wave number
angleStep=0.05;
theta=0:angleStep:360;
elementNumber=16;
spacingMIN=0.5;
spacingMAX=2;
mainbeam=5;
thetaM=90;

count=1;
%% random search
for m=1:10000
    d=spacingMIN + (spacingMAX-spacingMIN).*rand(1,elementNumber-1);
    d=[0,d];
    % d=[spacingMIN/2,spacingMIN*ones(1,elementNumber/2-1)];
    
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
    cvx_begin
    variable w(elementNumber) complex
    minimize( max(abs(A_SL*w)) )
    subject to
    A_M*w==1;
    cvx_end
    
    if max(abs(A_SL*w))<10^(-10/20)
        resultd(count,:)=d;
        resultw(count,:)=w';
        resultSLL(count)=max(abs(A_SL*w));
        save('resultd.mat', 'resultd');
        save('resultw.mat', 'resultw');
        save('SLL.mat', 'resultSLL');
        count=count+1;
        %plot(theta,20*log10(abs(A*x(2:9))));
        %axis([0,180,-70,0]);
        %hold on;
    end
end

%% Plot result
%w=ones(1,elementNumber)';
% plot(theta,20*log10(abs(A*w))-max(20*log10(abs(A*w))));
% axis([0,180,-30,0]);
