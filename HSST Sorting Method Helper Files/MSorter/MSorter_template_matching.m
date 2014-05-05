function sortCode = MSorter_template_matching(filename, template,Fr, spikes, index,duration,Thr_peak,dist_sig,corr_sig)
% close all
dur = size(template,2);
tempclass = size(template,1);
corrmatrix = zeros(size(spikes,1),tempclass);
idx = zeros(tempclass,size(spikes,1));
distmaxr = ones(tempclass,size(spikes,1));

c = [];
for i = 1:size(spikes,1)
    duration(i) = dur;
    if max(spikes(i,:))-min(spikes(i,:)) > Thr_peak 
        r1 = corrcoef([template(:,dur/2-duration(i)/2+1:dur/2+duration(i)/2)',spikes(i,dur/2-duration(i)/2+1:dur/2+duration(i)/2)'*ones(1,tempclass)]);
        r2 = r1(tempclass+1,1:tempclass);
        dist1 = (template(:,dur/2-duration(i)/2+1:dur/2+duration(i)/2) - ones(tempclass,1)*spikes(i,dur/2-duration(i)/2+1:dur/2+duration(i)/2))*(template(:,dur/2-duration(i)/2+1:dur/2+duration(i)/2) - ones(tempclass,1)*spikes(i,dur/2-duration(i)/2+1:dur/2+duration(i)/2))';
        dist = diag(dist1);
        v = max(r2(dist < dist_sig));
        distmaxr(dist < dist_sig,i) = dist(dist < dist_sig);
        corrmatrix(i,dist < dist_sig) = r2(dist < dist_sig);
        [c(i) d] = min(distmaxr(:,i));
        if ~isempty(v)
            if v>corr_sig && c(i)<1
                idx(d,i) = 1;
            end
        end
    end
end

% figure
% subplot(2,1,1)
%     plot(template')
% subplot(2,1,2)
%     hist(c,[0:0.01:10])

% figure
% for i = 1 : tempclass
%     subplot(2,ceil(tempclass/2),i)
%         if sum(idx(i,:)==1)~=0
%         plot([1:dur]/Fr,spikes(idx(i,:)==1,:)','b')
%         hold on
%         end
%     axis([0 dur/Fr -2e-4 1e-4])
% end

% fname4 = strcat(filename,'_spike_cluster');
% saveas(gcf,fname4,'jpg')


% figure(10)
% for i = 1:tempclass
%     subplot(2,ceil(tempclass/2),i)
%     plot([1:dur]/Fr,template(i,:),'linewidth',3)
%     axis([0 dur/Fr -1.5e-4 1e-4])
%     text(0.2,0.5e-4,num2str(size(find(idx(i,:)==1),2)/max(index)*Fr*1000));
%     hold on
% end
% title(['detection rate',num2str(1-size(find(sum(idx)==0),2)/size(idx,2))])

% save figure
% saveas(gcf,filename,'jpg')
% fname = strcat(filename,'.mat');
% save(fname,'spikes','template','index','idx')
% close all
sortCode = idx;

