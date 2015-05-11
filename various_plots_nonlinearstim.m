%% run sim and plot results with fit weights        
Estartpos=-10
init(:,1)=HalfWaveRectify((r0(1:n)+Estartpos*slopes(1:n)));
% make a perturbation   
numperturb=0;
 toperturb=ceil(rand(1,numperturb)*n);
        init(toperturb,1)=zeros(1,numperturb);
        
 timemult=1;
ymax=100;        
inits=sinfE(init(:,1));
[R,INT,U,S]=nonlinear_network_simulation(n,dub,tonic,slopes,init(:,1),init(:,1),inits,0:delt:(duration*timemult)); 
figure(2)
subplot(1,3,1)
plot(R(1:1000,:))
title('real network activity')
ylim([0 ymax])
[RC,INT,U,S]=nonlinear_network_simulation(n,betac,tonic,slopes,init(:,1),init(:,1),inits,0:delt:(duration*timemult)); 
figure(2)
subplot(1,3,2)
plot(RC(1:1000,:))
title('W*cp fit network activity')
ylim([0 ymax])
[RC0,INT,U,S]=nonlinear_network_simulation(n,betac0,tonic,slopes,init(:,1),init(:,1),inits,0:delt:(duration*timemult)); 
figure(2)
subplot(1,3,3)
plot(RC0(1:1000,:))
title('W*c fit network activity')
ylim([0 ymax])

figure(5)
subplot(1,2,1)
plot(R(1:1000,:)-RC(1:1000,:))
title('Actual vs. W*cp simulated activity')
subplot(1,2,2)
plot(R(1:1000,:)-RC0(1:1000,:))
title('Actual vs. W*c simulated activity')
