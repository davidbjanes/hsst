function spl = Testing_DPBPSVD_Time_depedent_model_v4(X,K,H_K,TuningParameter1,burnin,num,space,debug)
%
%
%   INPUTS
%   =======================================================================
%   X       : (samples x snippet_example)
%   K       :
%   H_K     :
%   TuningParameter1
%   burnin  :
%   num     :
%   space   :
%   debug   :
%
%
%   IN CODE
%   H_z : array, value for each snippet, value indicates cluster of snippet
%
%
%%% Hierarchical Dirichelet Process and Beta Process Pseudo-SVD via Markov Chain Monte Carlo (HDPBPSVD_MCMC)
%%% written by Bo Chen, 4.19
%%%%% Please cite: Bo Chen, David Carlson and Lawrence Carin,
%%%%% 'On the Analysis of Multi-Channel Neural Spike Data', NIPS 2011.
%----------parameter setting----------------------
%G_w0=10^2*eye(K); d0 = 10^-4*ones(p,1); good result on HC-1 data
%-------------------------------------------------
% CRS20120914 original name Testing_DPBPSVD_Time_depedent_model
% renamed to DPBPFA_TIME_v6.m
rng('default')
rng(0)


%JAH: Really????
n=size(X,2);
Len=n;

%JAH FIX: Use bsxfun instead
X = bsxfun(@minus,X,mean(X,2));
%X=X-mean(X,2)*ones(1,n);
p=size(X,1);
K=min(K,p);

c0 = 10^-0*ones(p,1);
d0 = 10^-6*ones(p,1);
e0 = 1/K*ones(K,1);
f0 = 1*(K-1)/K*ones(K,1);
g0 = 1e-0*ones(p,K);
h0 = 1*10^-6*ones(p,K);
w0 = 0.01;
w0_inv = 1/w0;
eta0 = 0;
w0_inv_x_eta0 = w0_inv*eta0;
pai = betarnd(e0,f0);

[U0,H0,V0] = svd((X),'econ');
A    = U0(:,1:K);
S    = H0(1:K,1:K)*V0(:,1:K)';
RSS  = S*S';
phi  = TuningParameter1*ones(p,1);
PC_z = rand(Len,H_K);
PC_z = PC_z./repmat(sum(PC_z,2),1,H_K);


% s = rng; %Record random number state
%Generate random number weighted by PC_z
H_z = zeros(1,Len);
for i=1:Len
    H_z(i)=randsample(H_K,1,true,PC_z(i,:));
end

%This is from randsample()
% %CODE ACTUALLY BEING USED
% %===================================
% %p = w(:)' / sumw;
% %edges = min([0 cumsum(p)],1);
% %[~, y] = histc(rand(k,1),edges);


% rng(s) %reapply random number state
% %MY CODE
% %==================================
% p_rand_temp = [zeros(Len,1) cumsum(PC_z,2)];
% p_rand_temp(:,end) = 1;
% 
% rand_n_temp = rand(Len,1);
% 
%TODO: Finish this, 
%H_z2 = asdfasdfasdf
%
% if ~isequal(H_z,H_z2)
%   error('Nope')
% end






w = rand(K,1);
z = rand(K,1)<pai;
alpha = gamrnd(g0,1./(h0));
G_v0  = K;
G_w0  = eye(K);
G_w0_D = cholcov(G_w0,1); %NOTE: technically this equals eye(K)
G_beta0 = 1;

%NOTE: These should really be 3d
G_lamda(1:H_K) = {eye(K)};
G_mu(1:H_K)    = {zeros(K,1)};


%%%%the parameters in the prior Beta of v
H_a = 1;   %1e-0;
H_b = 0.1; %1e-1;
HV_lamda = gamrnd(H_a,1/H_b);

%TODO: When done with changes, uncomment and remove code below
%Preferred Approach, but it throws off random number generator
% H_v = betarnd(ones(1,H_K),HV_lamda);
% H_v(H_K)=1;

%OLD APPROACH - remove when done with all other changes
%==================================
H_v = zeros(1,H_K);
for k=1:H_K-1
    H_v(k)=betarnd(1,HV_lamda);
end
H_v(H_K)=1;



%%%%the parameters in the prior Gamma of lamda
maxit = burnin + num*space;

acc=0;

mse  = zeros(1,maxit);
numC = zeros(1,maxit);
numZ = zeros(1,maxit);

for iter = 1:maxit
    
    tic
    G = A'.*(ones(K,1)*phi');
    E = G*A;
    F = G*X;
    F_S = diag(F*S'); %It would be nice if I didn't have to compute
    %everything then take the diagonal ...
    
    %Update z - Would be great to throw into a function
    %-----------------------------------------------------------
    %
    %   INPUTS
    %   ==================================
    %   F : 32 x 4454
    %   S : 32 x 4454
    %   E : 32 x 32
    %  RSS: 32 x 32

    tmprr_p1 = log(pai+eps) - log(1-pai+eps)-0.5.*w.^2.*diag(E).*diag(RSS);
    for k=1:K
        zw      = z.*w;
        zw(k)   = 0;
        midval2 = F_S(k)-E(k,:)*(zw.*RSS(:,k));
        tmprr = tmprr_p1(k) + w(k)*midval2;
        if z(k)
            z(k)=binornd(1,1-min(1,exp(-tmprr)));
        else
            z(k)=binornd(1,min(1,1/exp(-tmprr)));
        end
    end
    
    pai = betarnd(e0 + z, f0 + 1 - z);
    
    I_z   = find(z');
    
    %Update w - Would be great to throw into a function
    %-------------------------------------------------------------
    w = update_w(w,K,w0_inv,E,RSS,z,F_S,w0_inv_x_eta0,I_z,w0,eta0);
    
    
    
    
    
    
    %%%%%%%%%%%%%%%%%%   DP Part %%%%%%%%%%%%%%%%
    %------------------------------------------------------------
    %Chinese Restaurant Formulation of Dirichilet Process (DP)
    %
    
    %DETERMINATION OF H_z - TODO: Make this a function
    %---------------------------------------------------------

    
    H_z = compute_Hz(S,Len,H_K,H_v,G_mu,G_lamda,I_z);

    %DETERMINATION OF H_v - TODO: Make this a function
    %---------------------------------------------------------
    Hz_k_idx_all  = cell(1,H_K);
    [uHz,uHz_idx] = unique2(H_z);
    %NOTE: This is to make sure we have indices for all values from
    %1 to H_K, not just the unique values that are present
    %the uHz are values from 1 to H_K, and their value can be used
    %as an index into the original cell array (Hz_k_idx_all)
    Hz_k_idx_all(uHz) = uHz_idx;
    
    n_Hz = cellfun('length',Hz_k_idx_all);
    n_Hz_greater = Len - cumsum(n_Hz);
    
    %Concatentation & indexing done to make
    %random # generators match for comparision with old
    %Would be better to make this two lines and not index
    %and then just replace last value with 1
    
    
    HV_lamda = gamrnd(H_a+H_K-1, 1/(H_b-sum(log(1-H_v(1:H_K-1)+realmin))) );
    %H_v = [betarnd(0.001 + n_Hz(1:end-1),HV_lamda + n_Hz_greater(1:end-1)) 1];
    
    %Keep this for now to be the same :/
    %This function calls random twice, which means that
    %this won't match the old code without the loop (hence the commented
    %out line above)
    H_v = ones(1,H_K); %NOTE: We throw in the 1 at end using the ones function
    for k = 1:H_K-1
        H_v(k) =  betarnd(0.001 + n_Hz(k),HV_lamda + n_Hz_greater(k));
    end
    

    
    [G_mu,G_lamda] = updateModel(G_mu,G_lamda,n_Hz,G_w0,G_v0,G_w0_D,S,Hz_k_idx_all,K,H_K,G_beta0);

    zw     = z.*w;
    
    %%%%%%%%%% Sampling S %%%%%%%%
    %%-------------------update weight------------------
    S = updateS(S,X,zw,A,phi,K,Len,Hz_k_idx_all,G_lamda,G_mu,n_Hz);
    

    RSS = S*S';
    
    
    %Computing A
    %--------------------------------------------------------
    tmp1 = randn(p,K);
    %NOTE: Since I_z is sparse, I didn't do much to pull things
    %out of the loop
    for k = I_z
        A(:,k) = 0;
        Xm     = zw(k)*phi.*(X*S(k,:)' - A*(zw.*RSS(:,k)));
        tmpA1  = phi*(zw(k)^2*RSS(k,k));
        tmpA3  = 1./(alpha(:,k) +tmpA1);
        tmpA2  = tmpA3.*Xm;
        A(:,k) = tmpA2 + tmpA3.^(0.5).*tmp1(:,k);
    end
    mask = ~z;
    A(:,mask) = tmp1(:,mask)./(alpha(:,mask).^0.5);
    
    
    % There might be numerical problems for small alpha.
    alpha     = gamrnd(g0 + 0.5, 1./(h0 + 0.5*A.^2));
    res       = X - A(:,I_z)*(diag(zw(I_z))*S(I_z,:));
    phi       = gamrnd(c0+0.5*n, 1./(d0+0.5*sum(res.^2,2)));
    mse(iter) = mean(sqrt(sum((X - A*diag(z.*w)*S).^2,1)));
    
    ColC  = find(n_Hz > 50);
    uniqC = find(n_Hz > 0);
    
    numC(iter) = length(uniqC);
    numZ(iter) = sum(z);
    

    %Log Stuff
    %----------------------------------------------------------------------
    if iter > burnin
        ndx = iter - burnin;
        if mod(ndx,space) == 0
            spl.Likelihood(ndx) = Loglikelihood_v4(S,G_lamda,G_mu,H_z,H_K);
            spl.H_z{ndx} = H_z;
            spl.numCluster(ndx)=length(ColC);
        end
    end

    %Print Things
    %----------------------------------------------------------------------
    if mod(iter,25) == 0,
        if (debug)
            % fprintf('iter %d: mse=%g: numDic=%g  Gclassnum=%d FirstOpenum: %d%n',iter,mse(iter),sum(z(:,1)),numC(iter),FirstOpenum);
            fprintf('iter %d: mse=%g: numDic=%g  Gclassnum=%d%n',iter,mse(iter),sum(z(:,1)),numC(iter));
            toc
        else
            fprintf('.')
        end
    end
    
end

%NOTE: Only used at end for logging
spl.H_v = H_v;
spl.HV_lamda = HV_lamda;
spl.A = A;
spl.S = S;
spl.w = w;
spl.z = z;
spl.phi = phi;
spl.alpha = alpha;
spl.G_mu  = G_mu;
spl.G_lamda = G_lamda;
spl.numC = numC;
spl.pai  = pai;
spl.acc=acc;

end

function w = update_w(w,K,w0_inv,E,RSS,z,F_S,w0_inv_x_eta0,I_z,w0,eta0)
%
%
%
%   INPUTS
%   ======================================
%   
%

    %NOTE: This code removed the call to TNrnd
    %The equations were also simplified assuming:
    %a = 1
    %b = inf
    
    
    r_local = rand(1,K);

    %NOTE: We only iterate over non-zero entries
    for k = I_z
        w(k) = 0;
        tmpb = 1./(w0_inv+E(k,k)*RSS(k,k));
        tmpa = tmpb*(F_S(k)-E(k,:)*(z.*w.*RSS(:,k)) + w0_inv_x_eta0);

        %Old code for reference
        %Avoiding call to normcdf
        %PHIl_old = normcdf(-tmpa/sqrt(tmpb)); 
        PHIl = 0.5 * erfc(tmpa/sqrt(2*tmpb));

        w(k) = tmpa + sqrt(tmpb)*(sqrt(2)*erfinv(2*(PHIl + (1-PHIl)*r_local(k)) - 1));
        
    end

    %This is for z(k) = 0
    %tmpb and tmpa simplify to the following
    %NOTE: We update w in this loop, so w would change
    %but we also multiply w by z, which doesn't change
    %which means we can take this stuff out of the loop
    %otherwise we would need to keep it in for calculating tmpa
    tmpb = w0;
    tmpa = eta0;
    sigma_local = sqrt(tmpb);
    PHIl = normcdf(-tmpa/sigma_local);
    w(~z) = tmpa + sigma_local.*(sqrt(2).*erfinv(2*(PHIl + (1-PHIl).*r_local(~z)) - 1));



end

function S = updateS(S,X,zw,A,phi,K,Len,Hz_k_idx_all,G_lamda,G_mu,n_Hz)


    zw_A   = bsxfun(@times,zw,A');
    %zwa    = zw_A.*repmat(phi,1,K)';
    zwa    = bsxfun(@times,zw_A,phi');
    
    
    
    midval = zwa*zw_A';  %Local
    AX     = zwa*X;      %Local
    

    temp0  = randn(K,Len);
    for clusndx = find(n_Hz ~= 0)
        ndx         = Hz_k_idx_all{clusndx};
        S_mupart    = bsxfun(@plus,AX(:,ndx),G_lamda{clusndx}*G_mu{clusndx});
        
        %S_sigmapart = midval + G_lamda{clusndx};
        G           = chol(midval + G_lamda{clusndx});
        S(:,ndx)    = G\(temp0(:,ndx) + (G')\S_mupart);
    end


end

function [G_mu,G_lamda] = updateModel(G_mu,G_lamda,n_Hz,G_w0,G_v0,G_w0_D,S,Hz_k_idx_all,K,H_K,G_beta0)

    %---------------------------------------
    %NOTE: I didn't get to optimize too much in this section ...

    tmp0 = randn(K,H_K);


    cov_norm_factor = 1./(n_Hz - 1);
    for k=1:H_K
        
        kM     = n_Hz(k);
        if kM == 0
            G_mu{k}    = randn(K,1);
            G_lamda{k} = wishrnd(G_w0,G_v0,G_w0_D);
        else
            colS = S(:,Hz_k_idx_all{k}); %column of S
            
            averS = mean(colS,2);
            
            %Do covariance locally since we have the mean ...
            if any(size(colS) == 1)
               %I don't trust myself with this
               %I think TMW is wrong though ...
               %G_s = zeros(size(colS,1));
               G_s   = cov(colS');
            else
               %Copied from cov() function
               %m   = size(colS,2);
               xc  = bsxfun(@minus,colS,averS);
               G_s = cov_norm_factor(k)*(xc * xc');% / (m-1);
            end

            G_beta = G_beta0 + kM;
            G_v    = G_v0 + kM;
            
            invG_w = G_w0 + G_s*(kM-1) + G_beta0*kM/G_beta*averS*(averS)';
            
            G_w        = inv(invG_w);
            G_w        = 0.5*(G_w+G_w');
            G_lamda{k} = wishrnd(G_w,G_v);
            
            G_mumu    = kM*averS/G_beta;
            Sinvsigma = chol(inv(G_lamda{k}*G_beta))';
            G_mu{k}   = Sinvsigma*tmp0(:,k) + G_mumu;
        end
    end



end

function H_z = compute_Hz(S,Len,H_K,H_v,G_mu,G_lamda,I_z)

 PC_zz = zeros(Len,H_K);
    
    sum_log_Hv_plus_log_Hv = [0 cumsum(log(1 - H_v(1:end-1)))] + log(H_v);
    
    for k=1:H_K
        RR0        = chol(G_lamda{k}(I_z,I_z));
        xRinv0     = bsxfun(@minus,S(I_z,:),G_mu{k}(I_z))'*RR0;
        quadform0  = sum(xRinv0.^2, 2);
        PC_zz(:,k) = sum_log_Hv_plus_log_Hv(k) + sum(log(diag(RR0)))-0.5*quadform0;
    end

    
    
    
    PC_Z = exp(bsxfun(@plus,PC_zz,-max(PC_zz,[],2)));
    PC_Z = bsxfun(@rdivide,PC_Z,sum(PC_Z,2));
    
    H_z  = discreterndv3(PC_Z');



end