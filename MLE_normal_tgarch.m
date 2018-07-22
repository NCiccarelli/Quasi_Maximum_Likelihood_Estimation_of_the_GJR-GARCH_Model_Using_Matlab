function [ sum_lik  ] = MLE_normal_tgarch( alfa_plus_hat,  alfa_minus_hat  ,betahat, sigma, y,h )
% 

likelihoods= zeros(1,1999); 
T=2000;
for  j=2: T;    

     h(j)= 1+ betahat .* h(j-1) +  alfa_plus_hat .* ((y(j-1)).^2) .* ((y(j-1))>0) + alfa_minus_hat .* ((y(j-1)).^2) .* ((y(j-1))<0)   ;

    likelihoods(j) =   - 0.5 .*(log(( 2.*pi ))) -0.5 .* (  (y(j)^2 ) ./ ((sigma.^2) .* h(j)   ) ) - log(sigma ) - 0.5.* log(h(j))  ; 

end
sum_lik= -sum(likelihoods) ;
 
end

