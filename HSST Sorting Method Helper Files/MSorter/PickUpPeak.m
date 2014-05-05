function [final_detection_flag  final_detection_level final_detection_value]=PickUpPeak(detection_result,c)

%c : the matrix of wavelet transform coefficients


[c_m c_n] = size(c);
final_detection_flag=zeros(1,c_n);
final_detection_level=zeros(1,c_n);
final_detection_value=zeros(1,c_n);

% final_detection_result=0;
%% replace Yuan Sept 8 2011
% 
% for i=1: c_n % 3200 times
%     
%     non_zero=find(detection_result(:,i)>0, 1);
% 
%     if  ~isempty(non_zero)
% 
%         final_detection_result=zeros(c_m,1);
%         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         %final_detection_result(non_zero,1)=abs(c(non_zero,i)); % old
%         final_detection_result(:,1)=abs(c(:,i)); %new
%         
%         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         %temp=max(final_detection_result(:,i));
%         %temp= abs(c(non_zero,i)) ;
%         max_temp=max(final_detection_result);
%         
%         index = find(final_detection_result==max_temp);
%         final_detection_level(1,i)=index(1);
%         final_detection_flag(1,i)=1;
%         final_detection_value(1,i)=max_temp;
%         
%     end
% 
% end



t = find(detection_result > 0);
t = ceil(t/size(detection_result,1));
[max_temp1 index1] = max(abs(c(:,t)));
final_detection_level(1,t) = index1;
final_detection_flag(1,t)=1;
final_detection_value(1,t)=max_temp1;



    


