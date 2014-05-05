function [Y,Bmatrix]= alignment_function(wf,Minposition,varargin)
% alignment at the point of the minimal value
% input of spike is a cell
varargin = sanitizeVarargin(varargin);
DEFINE_CONSTANTS
extrap_method = 'zero';
END_DEFINE_CONSTANTS

[N,nSpikes]   = size(wf);
[min_wf,min_idx]  = min(wf);
Y       = zeros(N,nSpikes);
Bmatrix = zeros(N,nSpikes);
% calculate waveform derivative for the filling in post alignment data
dwf = dif2(wf);
for iiSpike = 1:nSpikes
    
    J = min_idx(iiSpike);
    
    if length(J) > 1;
        J(2:end)=[];
    end
    
    if J < Minposition
        Delta                   = Minposition-J;
        
%         Y(1:Delta,iiSpike)      = 0;
                Y(1:Delta,iiSpike)          = wf(Delta+1,iiSpike)-dwf(Delta+1,iiSpike)*[Delta:-1:1];
        Y(Delta+1:N,iiSpike)        = wf(1:N-Delta,iiSpike);
        Bmatrix(Delta+1:N,iiSpike)  = 1;
    elseif J==Minposition
        Y(:,iiSpike)         = wf(:,iiSpike);
        Bmatrix(1:N,iiSpike) = 1;
    elseif J>Minposition
        Delta = J - Minposition;
        Y(1:N-Delta,iiSpike)      = wf(Delta+1:N,iiSpike);
                Y(N-Delta+1:N,iiSpike)    = wf(N,iiSpike)+dwf(N,iiSpike)*[1:Delta];
%         Y(N-Delta+1:N,iiSpike) = 0;
        Bmatrix(1:N-Delta,iiSpike)= 1;
    end
end
