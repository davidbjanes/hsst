function Loglikelihood = Loglikelihood_v4(X_t,nu,mu,H_z_idx,N_Hz)
%
%
%   INPUTS
%   ==================================
%   X_t : data transposed
%   nu  : precision matrix
%   mu  : mean
%   H_z_idx : cell array of indices of X which belong to each
%             group
%   N_Hz : # of indices in each H_z_idx array element
%   

%OLD CODE
%=================================
% % % datacc=H_z;
% % %       for mm=1:H_K
% % %          idx=find(datacc==mm);
% % %          KM=length(idx);
% % %          RR=chol(nu{mm});
% % %          xRinv=(X(:,idx)-repmat(mu{mm},1,KM)).'*RR.';
% % %          midval2=sum(log(diag(RR)));  
% % %          quadform = sum(sum(xRinv.^2, 2));
% % %          PH_z(mm)=midval2*KM-0.5*quadform; 
% % %       end
% % %    Loglikelihood=(sum(PH_z));

Loglikelihood = 0;
for mm = find(N_Hz ~= 0)
    RR       = chol(nu{mm});
    midval2  = sum(log(diag(RR)));
    
    xRinv    = bsxfun(@minus,X_t(H_z_idx{mm},:),mu{mm}')*RR';
    quadform = sum(xRinv(:).^2);
    
    Loglikelihood = Loglikelihood + midval2*N_Hz(mm) - 0.5*quadform; 
end




