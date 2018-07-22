
tic 

% number of Monte Carlo replications: 
repetitions = 2000;


% initialize matrices for coefficients
QMLE_alpha_plus= NaN(1,repetitions);
QMLE_alpha_minus = NaN(1,repetitions);
QMLE_beta= NaN(1,repetitions);
QMLE_sigma = NaN(1,repetitions);

 
 
for w = 1: repetitions
 
%% sample size
N=2000;

%% innovation of the GJR-GARCH process:
t5_distr= trnd(5, [N,1] )';
epsilon =   t5_distr./ (sqrt(5/3)) ;
 
%% parameter values
alfa0_minus= 0.15 ; 
alfa0_plus= 0.075 ;
beta0= 0.8;


%% generate data:
h=zeros(1,N);
y=zeros(1,N);
h(1)= 1./ (1- 0.5.*alfa0_minus - 0.5.*alfa0_plus    -beta0);
y(1)= ((h(1)).^0.5) * (epsilon(1)); 
 
for  i=1: (N-1);
     h(i+1)= 1+ beta0 .* h(i) +  alfa0_plus .* ((y(i)).^2) .* ((y(i))>0) + alfa0_minus .* ((y(i)).^2) .* ((y(i))<0)   ;
        y(i+1)=  ( (h(i+1)).^0.5) * (epsilon(i+1))  ;
end
 



 
%% QML estimation

% initial parameter values:
startingvalues = [0.05;0.05; 0.05;     0.05  ];
lowerbound = [  0 ;  0  ; 0 ; 0.00000000001  ];
upperbound = [  3; 3 ;  0.999999999 ;   3  ];        

 
%% likelihood function: 

likelihood = @(x) MLE_normal_tgarch(x(1),x(2),x(3),x(4), y,h );

 
options    = optimset('fmincon');
options    = optimset(options, 'TolFun', 1e-006);
options    = optimset(options, 'LargeScale', 'off');
options    = optimset(options, 'MaxFunEvals', 1000);
options    = optimset(options, 'MaxIter', 400);
  
 
 
[PARMLE_QMLE, fval] = fmincon(likelihood,startingvalues,[],[],[],[],lowerbound ,upperbound , @mycon_threshold_garch ,options);
 
QMLE_alpha_plus(1,w)= PARMLE_QMLE(1,1);
QMLE_alpha_minus(1,w)= PARMLE_QMLE(2,1);
QMLE_beta(1,w)= PARMLE_QMLE(3,1);
QMLE_sigma(1,w)= PARMLE_QMLE(4,1);
 


w
 
end
 



% Parameters for QMLE estimator: 
 
Mean_alpha_PLUS_QMLE= mean(QMLE_alpha_plus,2)
Mean_alpha_MINUS_QMLE= mean(QMLE_alpha_minus,2)
Mean_betaQMLE= mean(QMLE_beta,2)
Mean_sigmaQMLE= mean(QMLE_sigma,2)

Std_dev_alpha_PLUS_QMLE=std(QMLE_alpha_plus)
Std_dev_alpha_MINUS_QMLE=std( QMLE_alpha_minus )
Std_dev_betaQMLE=std(QMLE_beta)
Std_dev_sigmaQMLE=std(QMLE_sigma) 
 






toc



