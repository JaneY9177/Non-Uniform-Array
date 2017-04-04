%%
%  Verification for Random Spacing Linear Array Side Lobe Optimization Steering to Other Angles
%
%  Version 1
%  Zhengyu Peng
%%

%% Configuration
wavelength=1;
k=2*pi/wavelength;   % wave number
angleStep=0.05;
theta=0:angleStep:360;
elementNumber=16;
spacingMIN=0.5; % minimum space between two elements
spacingMAX=1.6;   % maximum space between two elements
mainbeam=10;     % width of the main lobe
thetaM=45;      % location of the main lobe

%% random spacing
% d=[0, spacingMIN + (spacingMAX-spacingMIN).*rand(1,elementNumber-1)];
% % d=[spacingMIN/2,spacingMIN*ones(1,elementNumber/2-1)];
%
% for nn=2:length(d)
%     d(nn)=d(nn-1)+d(nn);
% end

% figure(1);
% plot(d,zeros(1,length(d)),'x');
% axis([0,d(length(d)),-1,1]);

%% check results
load('resultd.mat');
load('resultd.mat');
load('SLL.mat')
[r,c]=size(resultd);

for kk=1:r
    %kk=21;
    d=resultd(kk,:);
    % w=resultw(nn,:)';
    
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
        resultd_verify(count,:)=d;
        resultw_verify(count,:)=w';
        resultSLL_verify(count)=max(abs(A_SL*w));
        resultSLL_main(count)=resultSLL(kk);
        save('resultd_verify.mat', 'resultd_verify');
        save('resultw_verify.mat', 'resultw_verify');
        save('SLL_verify.mat', 'resultSLL_verify');
        save('SLL_main.mat', 'resultSLL_main');
        count=count+1;
        %plot(theta,20*log10(abs(A*x(2:9))));
        %axis([0,180,-70,0]);
        %hold on;
    end
    
end
