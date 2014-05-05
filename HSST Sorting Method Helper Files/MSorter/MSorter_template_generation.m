 function [template] = MSorter_template_generation(filename, Fr,spikes, tempclass, MaxV)
% get rid of small and big spike noise
spikes = spikes(max(abs(spikes),[],2)<MaxV,:);
dur = size(spikes,2);
 %% kmeans
k = 10;
[ans1,c]=kmeans(spikes,k);
% close all
for i = 1:k
    if sum(ans1==i) > 1
        pcaresult(i,:) = mean(spikes(ans1==i,:));
    else
        pcaresult(i,:) = spikes(ans1==i,:);
    end
end
pcaresult = c;
dim = size(pcaresult,1);

for i=1:dim
    [~,y]=min(pcaresult(i,:));
    dur1 = length(find((pcaresult(i,1:y))<((max(pcaresult(i,1:y)))*0.7+min(pcaresult(i,1:y))*0.3)));
    dur2 = length(find((pcaresult(i,y:dur))<((max(pcaresult(i,y:dur)))*0.7+min(pcaresult(i,y:dur))*0.3)));
    duration(i) = dur1+dur2;
%     hold on
end

    kmeansresult = kmeans(pcaresult, tempclass,'distance','sqEuclidean');
    template = zeros(tempclass,dur);
%     figure;
    for i = 1:tempclass
%         subplot(2,ceil(tempclass/2),i)
%         plot([1:dur]/Fr,(pcaresult((kmeansresult==i),:)'))
%         axis([0 dur/Fr -1.5e-4 1e-4])
%         hold on
        if length(mean(pcaresult((kmeansresult==i),:)))~=1
%             plot([1:dur]/Fr,mean(pcaresult((kmeansresult==i),:))','linewidth',3)
            template(i,:) = mean(pcaresult((kmeansresult==i),:));
        else
%             plot([1:dur]/Fr,(pcaresult((kmeansresult==i),:)'),'linewidth',3)
            template(i,:) = pcaresult((kmeansresult==i),:)';
        end
    end

%% save figure
% save figure
% fname2 = strcat(filename,'_template');
% saveas(gcf,fname2,'jpg')