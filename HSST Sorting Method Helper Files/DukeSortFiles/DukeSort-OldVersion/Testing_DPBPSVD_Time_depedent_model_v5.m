function spl = Testing_DPBPSVD_Time_depedent_model_v5(X,K,H_K,TuningParameter1,burnin,num,space,debug)
%
%
%   INPUTS
%   =======================================================================
%   X       : (samples x snippet_example)
%   K       : samples of the waveform to use, this is only used if this is
%             less than the # of samples available, i.e. K = min(K,size(X,1))
%   H_K     :
%   TuningParameter1: Initialiation of phi, controls noise variance
%   burnin  :
%   num     :
%   space   :
%   debug   :
%
%
%   IN CODE
%   ========================================================================
%   H_z : array, value for each snippet, value indicates cluster of snippet
%
%
%   OUTPUT
%   ========================================================================
%     spl.H_v = H_v;
%     spl.HV_lamda = HV_lamda;
%     spl.A = A;
%     spl.S = S;
%     spl.w = w;
%     spl.z = z;
%     spl.phi = phi;
%     spl.alpha = alpha;
%     spl.G_mu  = G_mu;
%     spl.G_lamda = G_lamda;
%     spl.numC = numC;  %# of clusters for each iteration
%     spl.pai  = pai;
%     spl.acc  = acc;
%
%
%
%%% Hierarchical Dirichelet Process and Beta Process Pseudo-SVD via Markov Chain Monte Carlo (HDPBPSVD_MCMC)
%%% written by Bo Chen, 4.19
%%%%% Please cite: Bo Chen, David Carlson and Lawrence Carin,
%%%%% 'On the Analysis of Multi-Channel Neural Spike Data', NIPS 2011.
%----------parameter setting----------------------
%G_w0=10^2*eye(K);
%d0 = 10^-4*ones(p,1);
%good result on HC-1 data
%-------------------------------------------------
% CRS20120914 original name Testing_DPBPSVD_Time_depedent_model
% renamed to DPBPFA_TIME_v6.m
%
%   See Also:
%       discreterndv3

MIN_CLUSTER_SIZE = 50; %For counting # of clusters during iterations
%Clusters with less samples than this during iterations are not counted as clusters
%NOTE: This does not impact results, but is only for debugging

rng('default')
rng(0)


n   = size(X,2);
Len = n;

X = bsxfun(@minus,X,mean(X,2));

p = size(X,1);
K = min(K,p);

c0   = 10^-0*ones(p,1);
d0   = 10^-6*ones(p,1);

g0   = 1e-0*ones(p,K);
h0   = 1*10^-6*ones(p,K);
w0   = 0.01;
eta0 = 0;


%pai = pi_k,
%a   = 1 %These are currently hardcoded constants below
%b   = 1
%NOTE:
%e0 = a/K ...
%f0 = b*(K-1)/K ...
e0   = 1/K*ones(K,1);
f0   = 1*(K-1)/K*ones(K,1);
pai  = betarnd(e0,f0);

%precalulated constants
%----------------------------------------
g0_plus_0p5   = g0 + 0.5;
c0_plus_0p5n  = c0 + 0.5*n;

phi_beta_empty_Iz_const = 1./(d0 + 0.5*sum(X.^2,2));

w0_inv        = 1/w0;
w0_inv_x_eta0 = w0_inv*eta0;

mse_Iz_empty  = mean(sum(X.^2,1).^0.5);

PHIl_Iz_empty = 0.5 * erfc(eta0./sqrt(2*w0));
sqrt_tmpb_2   = sqrt(w0*2);



[U0,H0,V0] = svd((X),'econ');
%NOTE: X = U0*H0*V0
A    = U0(:,1:K);
S    = H0(1:K,1:K)*V0(:,1:K)';
RSS  = S*S';


phi  = TuningParameter1*ones(p,1);


%NOT ACTUALLY USED ...
%========================================
%PC_z = rand(Len,H_K);
%PC_z = PC_z./repmat(sum(PC_z,2),1,H_K);
%NOTE: It turns out this isn't used ...
%H_z = initalize_Hz(PC_z,Len);

%Jim Match: lines added to make results the same:
%These can eventually be deleted ...
rand(Len,H_K);
rand(1,Len);


w       = rand(K,1);

%z = b_k = Bernoulli(pai)
z       = rand(K,1) < pai;

alpha   = gamrnd(g0,1./(h0));
G_v0    = K;
G_w0    = eye(K); %NOTE: CODE ASSUMES THIS IS SYMETTRIC ...
G_beta0 = 1;

%precalculated value for speed, copied from inside wishrnd()
G_w0_D  = cholcov(G_w0,1); %NOTE: technically this equals eye(K)

%NOTE: These should really be 3d
%Might be able to do optimizations then ...

%From normal-Wishart distribution
G_mu(1:H_K)    = {zeros(K,1)}; %Mean matrix
G_lamda(1:H_K) = {eye(K)};     %Precision matrix



%%%%the parameters in the prior Beta of v
H_a = 1;   %1e-0;
H_b = 0.1; %1e-1;

Ha_HK_minus_1 = H_a + H_K - 1;
HV_lamda = gamrnd(H_a,1/H_b);



%TODO: When done with changes, uncomment and remove code below
%Preferred Approach, but it throws off random number generator
%for comparision with old code when testing for equality ...
% H_v = betarnd(ones(1,H_K),HV_lamda);
% H_v(H_K)=1;

%OLD APPROACH - remove when done with all other changes
%=========================================================
H_v = zeros(1,H_K);
for k=1:H_K-1
    H_v(k) = betarnd(1,HV_lamda);
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
    
    G   = bsxfun(@times,A,phi)';
    E   = G*A;
    F   = G*X;
    F_S = diag(F*S');
    
    z   = update_z(z,pai,E,RSS,F_S,K,w);
    
    pai = betarnd(e0 + z, f0 + 1 - z);
    
    I_z = find(z'); %z is a boolean vector
    %NOTE: I_z represents indices where z is not false
    %It is also a row vector for use in for loops
    
    w = update_w(w,K,w0_inv,E,RSS,z,F_S,w0_inv_x_eta0,I_z,PHIl_Iz_empty,sqrt_tmpb_2,eta0);
    
    %Dirichilet Process (DP)
    %------------------------------------------------------------
    H_z      = compute_Hz(S,Len,H_K,H_v,G_mu,G_lamda,I_z);
    
    HV_lamda = gamrnd( Ha_HK_minus_1,   1/(H_b-sum(  log(1 - H_v(1:H_K-1) + realmin) ))   );
    
    [Hz_k_idx_all,n_Hz] = compute_Hz_stats(H_K,H_z);
    
    H_v = update_H_v(n_Hz,H_K,Len,HV_lamda);
    
    [G_mu,G_lamda] = update_Model(G_mu,G_lamda,n_Hz,G_w0,G_v0,G_w0_D,S,Hz_k_idx_all,K,H_K,G_beta0);
    
    zw = w;
    zw(~z) = 0;
    
    S   = update_S(S,X,zw,A,phi,K,Len,Hz_k_idx_all,G_lamda,G_mu,n_Hz);
    
    RSS = S*S';
    
    A   = update_A(A,I_z,z,p,K,X,S,zw,RSS,alpha,phi);
    
    
    % There might be numerical problems for small alpha.
    alpha   = gamrnd(g0_plus_0p5, 1./(h0 + 0.5*A.^2));
    if isempty(I_z)
        %When I_z is empty, res = X;
        %The beta constant is precomputed using that fact
        phi = gamrnd(c0_plus_0p5n,phi_beta_empty_Iz_const);
        mse(iter)  = mse_Iz_empty;
    else
        res =  X  -  A(:,I_z) * diag(zw(I_z))*S(I_z,:);
        phi = gamrnd(c0_plus_0p5n, 1./(d0 + 0.5*sum(res.^2,2)));
        mse(iter) = mean(sum((X - A*diag(zw)*S).^2,1).^0.5);
    end
    
    numC(iter) = sum(n_Hz > 0);
    numZ(iter) = length(I_z);
    
    %Log Stuff
    %----------------------------------------------------------------------
    if iter > burnin
        ndx = iter - burnin;
        if mod(ndx,space) == 0
            spl.Likelihood(ndx) = Loglikelihood_v4(S',G_lamda,G_mu,Hz_k_idx_all,n_Hz);
            spl.H_z{ndx}        = H_z;
            spl.numCluster(ndx) = sum(n_Hz > MIN_CLUSTER_SIZE);
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
%There is no need to populate these during the simulation
%========================================================
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
spl.acc  = acc;

end

function H_v = update_H_v(n_Hz,H_K,Len,HV_lamda)
%
%   NOTE: Code could run faster but has been
%   replaced with slower version to allow matching of previous code
%
%

%OLD CODE
%=======================================
% HV_lamda=gamrnd(H_a+H_K-1,1/(H_b-sum(log(1-H_v(1:H_K-1)+realmin))));
% for k=1:H_K-1
%     H_v(k)=betarnd( 10^-3+sum(H_z(:)==k),HV_lamda+sum(H_z(:)>k));
% end
% H_v(H_K)=1;


n_Hz_greater = Len - cumsum(n_Hz);

%Concatentation & indexing done to make
%random # generators match for comparision with old
%Would be better to make this two lines and not index
%and then just replace last value with 1
%BEST CODE (lines below)
%===============================

% H_v = betarnd(0.001 + n_Hz,HV_lamda + n_Hz_greater);
% H_v(end) = 1;

%Keep this for now to be the same :/
%This function calls random twice, which means that
%this won't match the old code without the loop (hence the commented
%out line above)

%CODE THAT MATCHES V3
%=========================================================================
H_v = ones(1,H_K); %NOTE: We throw in the 1 at end using the ones function
for k = 1:H_K-1
    H_v(k) =  betarnd(0.001 + n_Hz(k),HV_lamda + n_Hz_greater(k));
end

end

function A = update_A(A,I_z,z,p,K,X,S,zw,RSS,alpha,phi)
tmp1 = randn(p,K);
%NOTE: Since I_z is sparse, I didn't do much to pull things
%out of the loop

%OLD CODE
%=======================
% % % for k = 1:K
% % %     if z(k)==1
% % %         A(:,k) = 0;
% % %         zw= z.*w;
% % %         Xm = zw(k)*phi.*(X*S(k,:)' - A*(zw.*RSS(:,k)));
% % %         tmpA1=phi*(zw(k)^2*RSS(k,k));
% % %         tmpA3 = 1./(alpha(:,k) +tmpA1);
% % %         tmpA2 = tmpA3.*Xm;
% % %         A(:,k) = tmpA2 + tmpA3.^(0.5).*tmp1(:,k);
% % %     else  % Draw from base
% % %         A(:,k) = alpha(:,k).^(-0.5).*tmp1(:,k);
% % %     end
% % % end

%NOTE: A gets updated in this loop but zw does not
%values for which z(k) = 0 will have the value
%Xm equal to zero, meaning tmpA2 will be zero
%In addition, tmpA1 will also be zero, leaving
%no dependency of z(k) = 0 values on other values,
%and no dependency of z(k) = 1 values on A(:,k) values
%so they can be updated out of the loop
%NOTE: I_z is sparse, so pulling much more of the loop
%can actually be detrimental
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

end

function w = update_w(w,K,w0_inv,E,RSS,z,F_S,w0_inv_x_eta0,I_z,PHIl_Iz_empty,sqrt_tmpb_2,eta0)
%
%
%
%   INPUTS
%   ======================================
%

%OLD CODE
%======================================================
%     pai = betarnd(e0 + z, f0 + 1 - z);
%     for k = 1:K
%         w(k) = 0;
%         tmpb =  (w0^-1+z(k)^2*E(k,k)*RSS(k,k))^-1;
%         tmpa = tmpb*(z(k)*(F(k,:)*S(k,:)'-E(k,:)*(z.*w.*RSS(:,k)))+w0^-1*eta0);
%         w(k) = TNrnd(0,inf,tmpa,sqrt(tmpb),1);
%     end


%NOTE: This code removed the call to TNrnd
%The equations were also simplified ASSUMING:
%a = 1
%b = inf


r_local = rand(1,K);

%NOTE: We only iterate over non-zero entries
%Since in general I_z is sparse, this can
%lead to inefficiences if we take things out of the loop

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

%This is the calculation for the zero entries (i.e. z(k) = 0)
%=============================================================
%tmpb and tmpa simplify to the following
%NOTE: We update w in this loop, so w would change
%but we also multiply w by z, which doesn't change
%which means we can take this stuff out of the loop
%otherwise we would need to keep it in for calculating tmpa
% tmpb = w0;
% tmpa = eta0;
% sigma_local = sqrt(tmpb);
% %PHIl = normcdf(-tmpa/sigma_local);
% PHIl = 0.5 * erfc(PHIl_IZ_consttmpa/sqrt(2*tmpb));

w(~z) = eta0 + sqrt_tmpb_2 * erfinv(2*(PHIl_Iz_empty + (1-PHIl_Iz_empty).*r_local(~z)) - 1);



end

function [Hz_k_idx_all,n_Hz] = compute_Hz_stats(H_K,H_z)
Hz_k_idx_all  = cell(1,H_K);
[uHz,uHz_idx] = unique2(H_z);
%NOTE: This is to make sure we have indices for all values from
%1 to H_K, not just the unique values that are present
%the uHz are values from 1 to H_K, and their value can be used
%as an index into the original cell array (Hz_k_idx_all)
Hz_k_idx_all(uHz) = uHz_idx;

n_Hz = cellfun('length',Hz_k_idx_all);

end

function S = update_S(S,X,zw,A,phi,K,Len,Hz_k_idx_all,G_lamda,G_mu,n_Hz)
%
%
%
%   JAH STATUS: Optmized, needs documentation


zw_A   = bsxfun(@times,zw,A');
zwa    = bsxfun(@times,zw_A,phi');
midval = zwa*zw_A';
AX     = zwa*X;

temp0  = randn(K,Len); %VERY SLOW LINE
for clusndx = find(n_Hz ~= 0)
    ndx         = Hz_k_idx_all{clusndx};
    S_mupart    = bsxfun(@plus,AX(:,ndx),G_lamda{clusndx}*G_mu{clusndx});
    
    %NOTE: Combined two lines to save time
    %S_sigmapart = midval + G_lamda{clusndx};
    G           = chol(midval + G_lamda{clusndx});
    S(:,ndx)    = G\(temp0(:,ndx) + (G')\S_mupart); %VERY SLOW LINE
end


end

function [G_mu,G_lamda] = update_Model(G_mu,G_lamda,n_Hz,G_w0,G_v0,G_w0_D,S,Hz_k_idx_all,K,H_K,G_beta0)
%
%
%   JAH STATUS: Optimized
%   Needs documentation
%

%OLD CODE
%==============================================================
% % % %  tmp0=randn(K,H_K);
% % % %         for k=1:H_K
% % % %             h_pos=find(H_z==k);
% % % %             colS=S(:,h_pos);
% % % %             kM=size(colS,2);
% % % %             if kM==0
% % % %                 G_mu{k}= randn(K,1);
% % % %                 G_lamda{k}= wishrnd(G_w0,G_v0);
% % % %             else
% % % %                 averS=mean(colS,2);
% % % %                 G_s=cov(colS');
% % % %                 G_beta=G_beta0+kM;
% % % %                 G_v=G_v0+kM;
% % % %
% % % %                 invG_w=G_w0+G_s*(kM-1)+G_beta0*kM/G_beta*averS*(averS)';
% % % %
% % % %                 G_w=inv(invG_w);
% % % %                 G_w=(G_w+G_w')/2;
% % % %                 G_lamda{k}= wishrnd(G_w,G_v);
% % % %
% % % %                 G_mumu=kM*averS/G_beta;
% % % %                 Sinvsigma = chol(inv(G_lamda{k}*G_beta))';
% % % %                 %             G_mu{k} = Sinvsigma*randn(K,1)+G_mumu;
% % % %                 G_mu{k} = Sinvsigma*tmp0(:,k)+G_mumu;
% % % %             end
% % % %         end



tmp0            = randn(K,H_K);
%cov_norm_factor = 1./(n_Hz - 1); No longer used, since invG_w multiplied by this ...
G_beta_all      = G_beta0 + n_Hz;
G_v_all         = G_v0 + n_Hz;
kM_over_Gbeta   = n_Hz./G_beta_all;
sz_inv          = 1./n_Hz;

for k = 1:H_K
    kM = n_Hz(k);
    if kM == 0
        G_mu{k}    = randn(K,1);
        %NOTE: By passing in the third argument we save time
        G_lamda{k} = wishrnd(G_w0,G_v0,G_w0_D);
    else
        colS  = S(:,Hz_k_idx_all{k}); %column of S
        averS = sum(colS,2)*sz_inv(k);
        
        %Do covariance locally since we have the mean ...
        %and we can't pass the precomputed mean into cov()
        if any(size(colS) == 1)
            %I don't trust myself with this
            %I think TMW is wrong though ...
            %G_s = zeros(size(colS,1));
            G_s   = cov(colS');
            G_s_scaled   = G_s*(kM-1);
        else
            %Copied from cov() function
            %-------------------------------------------
            %m   = size(colS,2); %Don't need to know this because of bsxfun
            %and cov_norm_factor which was pulled out of the loop
            xc  = bsxfun(@minus,colS,averS);
            G_s_scaled = xc * xc'; %NOTE: We don't divided by (kM-1),
            %because we would multiply it later, so we skip both steps
        end
        
        G_mumu    = kM_over_Gbeta(k)*averS;
        invG_w    = G_w0 + G_s_scaled + G_beta0*G_mumu*(averS)';
        G_w       = inv(invG_w);
        %G_w      = 0.5*(G_w+G_w'); %Does nothing for symmetric matrices ...
        
        G_lamda{k} = wishrnd(G_w,G_v_all(k));
        
        Sinvsigma = chol(inv(G_lamda{k}*G_beta_all(k)))';
        G_mu{k}   = Sinvsigma*tmp0(:,k) + G_mumu;
    end
end






end


%NOTE: This function is no longer being called
function H_z = initalize_Hz(PC_z,Len) %#ok<DEFNU>
%
%
%   ASSUMPTIONS
%   =============================
%   1) randsample uses weighted sampling (PC_z)
%   2) single value requested from randsample



%ORIGINAL CODE
%==============================================
% %Generate random number weighted by PC_z
% H_z = zeros(1,Len);
% for i=1:Len
%     H_z(i)=randsample(H_K,1,true,PC_z(i,:));
% end


%This is from randsample()
%I've implemented it below in vector format
%for the old code above
%===================================
%randsample code:
%--------------------------
%p = w(:)' / sumw;
%edges = min([0 cumsum(p)],1);
%[~, y] = histc(rand(k,1),edges);


%MY CODE
%==================================
%NOTE: I played around with this a little bit
%and it seemed to be faster in this dimension
%then working in the other dimension
p_rand_temp = [zeros(1,Len); cumsum(PC_z')];
p_rand_temp(end,:) = 1;

rand_n_temp = rand(1,Len);

%NOTE: histc could operate on matrix
%If we could translate data by rand_n_temp
%to be in the same format, we would be all set
%i.e. if they all had the same cutoff
%Will leave for now
H_z = sum(bsxfun(@le,p_rand_temp,rand_n_temp));


end

function z = update_z(z,pai,E,RSS,F_S,K,w)
%
%
%   NOTE: (refer to Chen 2011)
%   z   = b (binary vector which is multiplied by lambda, see Chen 2011)
%   pai = pi_k
%
%   INPUTS
%   ==================================
%   F : 32 x 4454
%   S : 32 x 4454
%   E : 32 x 32
%  RSS: 32 x 32
%
%   See Also:
%       binornd_l_quick
%
%   JAH STATUS: Finished optimizing, needs documenting


%This was taken out of the loop
tmprr_p1 = log(pai+eps) - log(1-pai+eps)-0.5.*w.^2.*diag(E).*diag(RSS);

%OLD CODE
%     for k=1:K
%         signZ=z(k);
%         z(k) = 0;
%         zw=z.*w;
%         midval1=w(k)^2*E(k,k)*RSS(k,k);
%         midval2=F(k,:)*S(k,:)'-E(k,:)*(zw.*RSS(:,k));
%         tmprr = log(pai(k)+eps) - log(1-pai(k)+eps)-1/2*midval1 + w(k)*midval2;
%         %        z(k) = binornd(1,1/(1+exp(-tmprr)));
%         if signZ==0
%             z(k)=binornd(1,min(1,1/exp(-tmprr)));
%         else
%             z(k)=binornd(1,1-min(1,exp(-tmprr)));
%         end
%     end


%NOTE: z(k) is updated in the loop, and midval2 depends on it

r = rand(1,K);
n_log_1_over_r = -log(1./r);

zw     = w;
zw(~z) = 0;
for k = 1:K
    zw(k) = 0;
    
    midval2 = F_S(k) - E(k,:)*(zw.*RSS(:,k));
    tmprr   = tmprr_p1(k) + w(k)*midval2;
    
    if ~z(k) %This is the MUCH more common case
        %min(1,1/exp(-tmprr)
        if tmprr < 0
            %OLD
            %z(k) = r(k) < (1/exp(-tmprr));
            z(k) = tmprr > n_log_1_over_r(k); %Precompute scaling ...
        else
            z(k) = true;
        end
    else
        %1-min(1,exp(-tmprr)
        if tmprr > 0
            z(k) = r(k) < (1 - exp(-tmprr));   %Could precompute here but not common occurence
        else
            z(k) = false;
        end
    end
    
    %NOTE: We set if to zero above
    %so if z(k) is false, zw(k) is already valid
    if z(k)
        zw(k) = w(k);
    end
    
end

end

% % % function r = binornd_l_quick(p)
% % % %ASSUMPTIONS
% % % %=======================
% % % %1) Written for n = 1
% % % %2) Written for length(n) = 1
% % % %3) Written for logical output (0 or 1)
% % %
% % % r = rand(1) < p;
% % %
% % % %Code from binornd
% % % %r = cast(sum(rand([sizeOut,n]) < p, ndimsOut+1), outType);
% % % end

function H_z = compute_Hz(S,Len,H_K,H_v,G_mu,G_lamda,I_z)
%
%
%
%
%   See Also:
%       discreterndv3

%OLD CODE
%==================================
% % % % zz=logical(double(z)*double(z)');
% % % % PC_zz=zeros(Len,1);
% % % % for k=1:H_K
% % % %     RR0 = chol(reshape(G_lamda{k}(zz(:)),sum(z),sum(z)));
% % % %     xRinv0 = (S(z,:) - repmat(G_mu{k}(z),1,n))' * RR0;
% % % %     quadform0 = sum(xRinv0.^2, 2);
% % % %     PC_zz(:,k)=sum(log(1-H_v(1:k-1))) + log(H_v(k))+sum(log(diag(RR0)))-0.5*quadform0;
% % % %
% % % % end
% % % % PC_Z=zeros(Len,H_K);
% % % % PC_Z= PC_zz;
% % % % PC_Z = exp(PC_Z + repmat(-max(PC_Z,[],2),[1 H_K])) ;
% % % %
% % % % PC_Z = PC_Z./repmat(sum(PC_Z,2),[1 H_K]);
% % % %
% % % % H_z  = discreterndv2(PC_Z');



%NOTE: We could pull this out of the loop since this is constant given H_v
sum_log_Hv_plus_log_Hv = [0 cumsum(log(1 - H_v(1:end-1)))] + log(H_v); %Sum of all previous plus log(H_v)


if isempty(I_z)
    %NOTE: For emtpy I_z this simplifies significantly
    temp  = exp(sum_log_Hv_plus_log_Hv - max(sum_log_Hv_plus_log_Hv));
    
    %Since everything is the same
    %avoid replicating the costly cumsum
    %Since all divisions are the same for 1:Len, we can use histc
    %to determine boundaries, instead of the a sum for each column (see
    %discreterndv3 for comparison)
    p_cum = cumsum(temp);
    p_cum = p_cum./p_cum(end);
    r     = rand(1,Len);
    [~,H_z] = histc(r,p_cum);
    H_z = H_z+1;
    
else
    
    S_iz = S(I_z,:)'; %No need for redundant indexing, also transpose now instead of in loop

    sld_r         = zeros(1,H_K);
    quadform0_all = zeros(Len,H_K);
    for k = 1:H_K
        RR0      = chol(G_lamda{k}(I_z,I_z));
        sld_r(k) = sum(log(diag(RR0)));
        
        xRinv0   = bsxfun(@minus,S_iz,G_mu{k}(I_z)')*RR0; %NOTE: Lamda is precision, so multiply instead of dividing
        quadform0_all(:,k) = sum(xRinv0.^2, 2);
    end
    
    PC_Z = exp(bsxfun(@plus,-0.5*quadform0_all,sld_r+sum_log_Hv_plus_log_Hv));

    %NOTE: I don't think that the max subtraction is needed
    %Shouldn't overflow on the probabilities and the cumsum
    %PC_Z = exp(bsxfun(@minus,PC_zz,max(PC_zz,[],2)));
    %Not needed, performed in function below
    %PC_Z = bsxfun(@rdivide,PC_Z,sum(PC_Z,2));
    
    H_z  = discreterndv3(PC_Z');
    
end





end

function a = wishrnd(sigma,df,d)
%Local copy of wishrnd
%
%   Use chol instead of cholcov

n = size(sigma,1);
if (nargin <3)
    d = chol(sigma);
else
    n = size(d,2);
end

% For small degrees of freedom, generate the matrix using the definition
% of the Wishart distribution; see Krzanowski for example
if (df <= 81+n) && (df==round(df))
    x = randn(df,size(d,1)) * d;
    
    % Otherwise use the Smith & Hocking procedure
else
    % Load diagonal elements with square root of chi-square variates
    a = diag(sqrt(chi2rnd(df-(0:n-1))));
    
    % Load upper triangle with independent normal (0, 1) variates
    I = ones(n*(n-1)/2,1);
    I(1+cumsum(0:n-2))= n+1:-1:3;
    I = cumsum(I);
    
    a(I) = randn(n*(n-1)/2,1);
    
    % Desired matrix is D'(A'A)D
    x = a(:,1:size(d,1))*d;
end

a = x' * x;

end