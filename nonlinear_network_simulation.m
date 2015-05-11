% nonlinear network simulation function
% updated 3-23 to use correlated noise, along the direction of eye movement
function [allr,int,allu,alls]=nonlinear_network_simulation(n,w,tonic,slopes,initr,initu,inits,t)

noiseamp=0;
% dynamic parameters.  in milliseconds (from p20-21 of supplementary)
    tauR = 10;
    tauS = 1000;
    tauU = 100;
    dt = t(2)-t(1); 
% high synaptic threshold
   RfI = 44; RfE = 44; thetaI = 6; thetaE = 6;
% high recruitment threshold
%   RfI = 40; RfE = 40; thetaI = 22; thetaE = 22; 
    
    
    
% functions
    % synaptic activation function in steady-state
    HalfWaveRectify = @(x) max(x,0);
    sinfI_helper = @(r) 1/(1-1/(1+exp(RfI/thetaI)))*( 1./(1+exp((RfI-r)/thetaI)) - 1/(1+exp(RfI/thetaI)) );
    sinfE_helper = @(r) 1/(1-1/(1+exp(RfE/thetaE)))*( 1./(1+exp((RfE-r)/thetaE)) - 1/(1+exp(RfE/thetaE)) );
    sinfI = @(r) sinfI_helper(HalfWaveRectify(r));
    sinfE = @(r) sinfE_helper(HalfWaveRectify(r));
    F = @(x) 80/0.8*HalfWaveRectify(x); % units of Hz/nA:  input nA, returns Hz
    FInverse = @(x) 0.8/80*HalfWaveRectify(x); % units of nA/Hz: so, if you input Hz, returns nA


%% do simulation
    allr=zeros(length(t),n);
    allu=zeros(length(t),n);
    alls=zeros(length(t),n);
    allr(1,:)=initr;
    allu(1,:)=initu;
    alls(1,:)=inits;
    zerospots=find(allr(1,:)==0);

    for step=2:length(t)
        allr(step,:) = HalfWaveRectify(allr(step-1,:) + dt/tauR*(-allr(step-1,:)+F(tonic'+(alls(step-1,:))*w)  ) + noiseamp*sqrt(dt)*randn*slopes(1:n)');
        allu(step,:) = allu(step-1,:) + dt/tauU*(-allu(step-1,:)+allr(step-1,:)) + 0*noiseamp*sqrt(dt)*randn(1,n);
        alls(step,:) = alls(step-1,:) + dt/tauS*(-alls(step-1,:)+sinfE(allu(step-1,:))) + 0*noiseamp/100*sqrt(dt)*randn(1,n);
    end 


    
%% calculate r to output
int=tauR*diff(allr)/dt+allr(1:end-1,:);
allr=allr(1:end-1,:);

end
