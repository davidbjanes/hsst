function [spikes index Rate_ASU duration] =MSorter_detection (Signal,L,Fr,mode,Thr,dur,N_1)

import sqm.sortMethod.MSorter.*

% assign memory for 200Hz firing rate
ASU_detect = zeros(floor(length(Signal)/24000*200),dur); 
Signal=Signal-mean(Signal);
Wid=[0.5  1.1]; % ms [1 2.5]

%define relevant scales for detection
W1=Wid(1);
W2=Wid(2);
W1=round(W1*Fr)+1;
W2=round(W2*Fr)+1;
W1=W1-rem(W1,2);    %make sure it is even
W=W1:2:W2;          %filters should be of even length
CWDM_no=0;
index=[];
duration = [];
lap = floor(Fr/2);
for k=1:floor(length(Signal)/N_1)
  
        final_detection_duration = zeros(1,N_1);    
        idx=1+(k-1)*N_1;
        if idx+N_1-1+lap<= length(Signal)
   
            data=Signal(idx:idx+N_1-1+lap);
            %get all coefficients
            c=cwt(data,W,'coif5');


            % generating the mask
            [m n]=size(c);
            detection_result=zeros(m-(L-1),n);
            corr = ones(m-(L-1),n);
            new_corr = zeros(m-(L-1),n);
            com_corr = zeros(m-(L-1),n);
            for i=1:m-(L-1)

                for loop1=1:L 
                    corr(i,:)=c(i+loop1-1,:).*corr(i,:);
                end


                Pcorr=sum(corr(i,:).^2);
                PW=sum(c(i,:).^2);
                new_corr(i,:)=corr(i,:).*sqrt(PW/Pcorr)*3;
                com_corr(i,:) = abs(new_corr(i,:)) - abs(c(i,:));
                detection_result(i,com_corr(i,:)>0) = 1;  %% change 1

%                 Sigmaj=median(abs(c(i,1:W(i):end)-mean(c(i,:))))/0.6745;
%                 Thj=Sigmaj*sqrt(2*log(n));     %hard threshold
%                 new_corr(i,:)=corr(i,:).*sqrt(PW/Pcorr)*3;
%                 com_corr(i,:) = abs(new_corr(i,:)) - abs(Thj);

                detection_result(i,com_corr(i,:)>0) = 1;  %% change 1

            end
% 
%                 Pcorr1=sum(corr.^2,2);
%                 PW1=sum(c(1:m-(L-1),:).^2,2);
%                 new_corr1=corr.*(sqrt(PW1./Pcorr1)*ones(1,N_1));
%                 com_corr1 = abs(new_corr1) - abs(c(1:m-(L-1),:));
%                 mask(com_corr1>0) = 1;  %% change 1

            % detection by the policy
            % for each colomum, find the mask's elements
            % whose values are 1 on consective levels  (at least 2 levles)

            %detection_result=MulScaleCombine(mask);
%             clear Pcorr PW corr new_corr com_corr% release memory
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            [detection_flag detection_level detection_value]...
                =PickUpPeak(detection_result,c); %% change 2
%             clear c
                [final_detection_flag]=...
                    FinalDetection_modified(detection_flag,detection_value);
                final_detection_flag = final_detection_flag(1:N_1);
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%

            %  detection_flag=sum(detection_result);
            %  final_detection_flag=zeros(size(detection_flag));
            %  temp=find(detection_flag~=0);
            %  final_detection_flag(temp)=1;

            %%%%%%%%%%%%%%%%%%55
            detection_d = detection_flag==1; %% change 3
            final_detection_dur = detection_level(detection_d);
            final_detection_duration(detection_flag==1) = W(final_detection_dur);
%             final_detection_flag =...
%                 MergeMinReso_modified(final_detection_flag, final_detection_level, final_detection_value,W);
            %%%%%%%%%%%%%%%%%%%%%%% %%change 3 Yuan Sept 8 2011
            
            %% threshold
            x = ones(1,N_1)*Thr;
            if strcmp(mode, 'neg')
                test = find(data<Thr);
                % look for local minimum
                x(test) = data(test);
                y =sort([find( x(2:end) < x(1:end-1))+1, find( x(2:end) > x(1:end-1))]);
                diffe = y(diff(y) == 0);
            elseif strcmp(mode, 'pos')
                test = find(data>Thr);
                % look for local maximum
                x(test) = data(test);
                y =sort([find( x(2:end) > x(1:end-1))+1, find( x(2:end) < x(1:end-1))]);
                diffe = y(diff(y) == 0);
            else
                test = find(abs(data)>Thr);
                x(test) = abs(data(test));
                y =sort([find( x(2:end) > x(1:end-1))+1, find( x(2:end) < x(1:end-1))]);
                diffe = y(diff(y) == 0);
            end
            if ~isempty(test)
                T = diffe; % spike time
            end
            
            A=find(final_detection_flag~=0);
            final_spk = [];

            
            if ~isempty(diffe) && ~isempty(A)
                for th_p = 1:length(A)
                    [B, IX] = sort([A(th_p),T]);
                    B1 = diff(B);
                    if find(IX==1) <= length(B1) && find(IX==1) > 1
                        [l1 l2] = min([B1(find(IX==1)-1),B1(IX==1)]);
                        if l1 < lap && (B(find(IX==1)+(l2-1.5)*2) > max([final_spk,-lap])+lap)
                            final_spk = [final_spk, B(find(IX==1)+(l2-1.5)*2)];
                            duration = [duration,final_detection_duration(A(th_p))];
                        end
                    elseif find(IX==1) == length(B) && find(IX==1) ~= 1
                        if B1(find(IX==1)-1) && (B(find(IX==1)-1)> max([final_spk,-lap])+lap)
                            final_spk = [final_spk, B(find(IX==1)-1)];
                            duration = [duration,final_detection_duration(A(th_p))];
                        end
                    elseif find(IX==1) == 1
                        if numel(B1)~=0
                            if B1(1) < lap && (B(2) > max([final_spk,-lap])+lap)
                                final_spk = [final_spk, B(2)];
                                duration = [duration,final_detection_duration(A(th_p))];
                            end
                        end
                    end
                end

            end
            
            ASU_TT = final_spk+idx-1;
            ASU_TT = ASU_TT(diff([0,ASU_TT])~=0);
            temp=length(ASU_TT);
            index(1+length(index):length(index)+temp)=ASU_TT;
            duration1 = floor(dur/2)-1;
            duration2 = dur-duration1-1;
            if max(ASU_TT) + duration2 > length(Signal)
                ASU_TT = ASU_TT(1:end-1);
                duration = duration(1:end-1);
            end
            if min(ASU_TT) - duration1 < 0
                ASU_TT = ASU_TT(2:end);
                duration = duration(2:end);
            end
            if ~isempty(ASU_TT) && min(ASU_TT) > 17
                ASU_detect(CWDM_no+1:CWDM_no+length(ASU_TT),:) = Signal(ones(length(ASU_TT),1)*[-duration1:duration2] + ASU_TT'*ones(1,duration1+duration2+1));
                CWDM_no = CWDM_no + length(ASU_TT);
            end
        end
end

spikes=ASU_detect(sum(abs(ASU_detect'))>0,:);
clear ASU_detect
Rate_ASU=(length(index)-1)/length(Signal)*Fr*1000;

