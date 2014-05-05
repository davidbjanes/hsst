function             [final_detection_flag]=...
    FinalDetection_modified(detection_flag,detection_value)



% x = detection_value;
% y =sort([find( x(2:end-1) > x(1:end-2)), find(x(2:end-1) > x(3:end))]);
% detection_index = 1 + y(diff(y) == 0);
% final_detection_flag(detection_index)=1;
% final_detection_level(detection_index)=detection_level(detection_index);
% final_detection_value(detection_index)=detection_value(detection_index);

n=length(detection_flag);
detection_flag(n+1)=0;
final_detection_flag=zeros(1,n);
t = find(detection_flag == 1);
if ~isempty(t)
    tt = diff(t);
    ad = t((tt~=1)) + 1;
    b =  t(find(tt~=1)+1);
    FirstPos = [t(1) b];
    LastPos = [ad t(end) + 1]-1;
    for index=1:length(FirstPos)
        if LastPos(index)-FirstPos(index)>0
            [~, temp_index]=max(detection_value(FirstPos(index):LastPos(index)));
            detection_index=temp_index-1+FirstPos(index);               
            final_detection_flag(detection_index)=1;
        end
    end
end

