classdef WaveClus_HSST < hsst.sorter
    
    properties (Constant)
        
        sortMethodLabel = 'WaveClus-HSST';
        defaultParameters = [2:5:22];
    end 
        
    methods (Static)
        
        function [sort_ids, wf, ts, property] = sort(wf, ts, Fs, thresh, ...
                                     ~, ~, raw_data, ...
                                     sortparameter)
                        
            % Create mock figure object: handles
            handles = struct;
            wave_clus_figure = figure('visible', 'off');
            set(wave_clus_figure, 'userdata', {});
            handles.wave_clus_figure = wave_clus_figure;
            
            % (wave_features_wc: Line 62) -> plotting
            handles.projections = axes;
            
            % (wave_features_wc: Line 9) -> toggle radio button
            handles.spike_shapes_button = uicontrol(wave_clus_figure, 'value', 1); 
            
            % (amp_detect: Line 37) -> updates status on main GUI
            handles.file_name = text('String', ''); 
            
            % (amp_detect_wc: Line 39) informs spike detection to threshold
            handles.datatype = 'CSC data';
           
                         
            % Are we sorting raw_data or already extracted snips?
            if ~isempty(raw_data),
                
                % Load parameters
                filename = '';
                handles.par = set_parameters_ascii(filename,handles);   
                handles.par.sr = Fs;

                % Ensures amp_detect does not plot continous data
                % (amp_detect_wc: Line 106)
                handles.flag = 0;
                handles.par.tmax = 'all';
            
                index_all=[];
                spikes_all=[];

                % that's for cutting the data into pieces
                for j=1:handles.par.segments  

                    % LOAD CONTINUOUS DATA
                    data = raw_data;
                    tsmin = (j-1)*floor(length(data)/handles.par.segments)+1;
                    tsmax = j*floor(length(data)/handles.par.segments);
                    x=data(tsmin:tsmax); 
                    clear data; 

                    % SPIKE DETECTION + HIGH PASS FILTERING WITH AMPLITUDE THRESHOLDING
                    [spikes,thresh,index]  = amp_detect_wc(x,handles); %detection with amp. thresh.
                    index = index + tsmin - 1;

                    index_all = [index_all index];
                    spikes_all = [spikes_all; spikes];
                end

                %spike times in ms.
                index = index_all *1e3/handles.par.sr; 
                spikes = spikes_all;
            
            else
                % Load parameters
                filename = '';
                handles.par = set_parameters_ascii_spikes(filename,handles);   
                handles.par.sr = Fs;
                
                num_spikes = length(ts);
                if (num_spikes == size(wf,1))
                    spikes = wf;
                else
                    spikes = wf';
                end
                
                % Time Stamp info
                index = ts;
                
            
                num_samples = size(spikes, 2);

                handles.par.detection = 'pos';
    %                 handles.par.w_pre = 29;
    %                 handles.par.w_post = 61;

                handles.par.w_pre = 10;  %number of pre-event data points stored
                handles.par.w_post = num_samples - 10; %number of post-event data points stored
                
                [spikes] = spike_alignment(spikes,handles);
            end
            
            
            USER_DATA = get(handles.wave_clus_figure,'userdata');
            USER_DATA{2} = spikes;
            USER_DATA{3} = index;
            set(handles.wave_clus_figure,'userdata',USER_DATA);            
            
            % Extract spike features.
            [inspk] = wave_features_wc(spikes,handles); 
            
            % change axis (wave_features_wc: Line 62) to hidden
            set(handles.wave_clus_figure, 'visible', 'off');

            if handles.par.permut == 'y'
                if handles.par.match == 'y';
                    naux = min(handles.par.max_spk,size(inspk,1));
                    ipermut = randperm(length(inspk));
                    ipermut(naux+1:end) = [];
                    inspk_aux = inspk(ipermut,:);
                else
                    ipermut = randperm(length(inspk));
                    inspk_aux = inspk(ipermut,:);
                end
            else
                if handles.par.match == 'y';
                    naux = min(handles.par.max_spk,size(inspk,1));
                    inspk_aux = inspk(1:naux,:);
                else
                    inspk_aux = inspk;
                end
            end

            % Path Management
            current_path = pwd;
            path_to_current_function = which('wave_features_wc');
            [pathstr, ~, ~] = fileparts(path_to_current_function);
            cd(fullfile(pathstr, 'SPC'));
            
            %Interaction with SPC
            handles.par.fname_in = 'tmp_data';
            fname_in = handles.par.fname_in;
            
            %Input file for SPC
            save([fname_in],'inspk_aux','-ascii'); 
            
            %filename if "save clusters" button is pressed
            handles.par.fnamesave = [handles.par.fname '_' filename(1:end-4)];   
            
            %Output filename of SPC
            handles.par.fname = [handles.par.fname '_wc'];   
            
            handles.par.fnamespc = handles.par.fname;

            [clu,tree] = run_cluster(handles);
            USER_DATA = get(handles.wave_clus_figure,'userdata');
            
            % Path Management
            cd(current_path);

            if exist('ipermut')
                clu_aux = zeros(size(clu,1),length(index)) + 1000;
                for i=1:length(ipermut)
                    clu_aux(:,ipermut(i)+2) = clu(:,i+2);
                end
                clu_aux(:,1:2) = clu(:,1:2);
                clu = clu_aux; 
                clear clu_aux;
                USER_DATA{12} = ipermut;
            end

            USER_DATA{4} = clu;
            USER_DATA{5} = tree;
            USER_DATA{7} = inspk;
            set(handles.wave_clus_figure,'userdata',USER_DATA);
            
            %Selects temperature. 
            if isempty(sortparameter) || isnan(sortparameter),
                temp = find_temp(tree,handles); 
            else
                temp = sortparameter;
            end

            if size(clu,2)-2 < size(spikes,1);
                classes = clu(temp(end),3:end)+1;
                if ~exist('ipermut')
                    classes = [classes(:)' zeros(1,size(spikes,1)-handles.par.max_spk)];
                end
            else
                classes = clu(temp(end),3:end)+1;
            end
            
            % Defines nclusters
            cluster_sizes=[];
            ifixflag=zeros(1,handles.par.max_clus);
            for i=1:handles.par.max_clus                                    
                eval(['cluster_sizes = [cluster_sizes length(find(classes==' num2str(i) '))];'])
            end

            i=1;
            while i<=min(max(classes),handles.par.max_clus);
                if isempty(classes(classes==i))
                    for k=i+1:handles.par.max_clus
                        classes(classes==k)=k-1;
                    end
                else
                    i=i+1;
                end
            end
            
            sizemin_clus = handles.par.min_clus;
            nclusters = length(find(cluster_sizes(:) >= sizemin_clus));

            
            % Defines classes
            clustered = [];
            cont=0;  
            for i=1:nclusters       
                eval(['class_temp = find(classes==' num2str(i) ');'])
                if ((ifixflag(i)==1) & (~isempty(class_temp)))
                    ifixflagc = 1;
                else
                    ifixflagc = 0;
                end    
                if ((length(class_temp) >= sizemin_clus) || (ifixflagc == 1))
                    cont=cont+1;        
                    eval(['class' num2str(cont) '= class_temp;'])
                    eval(['clustered = [clustered class' num2str(cont) '];'])
                end        
            end
            nclusters = cont;
            class0 = setdiff( 1:size(spikes,1), sort(clustered) );

            % Redefines classes
            classes = zeros(size(spikes,1),1);
            for i = 1:nclusters+1
                if ~ (isempty(class0) & i==1)
                    eval(['classes(class' num2str(i-1) ') = ' num2str(i-1) ';']);
                end
            end

            USER_DATA = get(handles.wave_clus_figure,'userdata');
            USER_DATA{6} = classes(:)';
            USER_DATA{8} = temp(end);
            
            
            % OUTPUTS
            sort_ids            = USER_DATA{6};
            wf                  = USER_DATA{2};
            ts                  = USER_DATA{3};
            property.USER_DATA  = USER_DATA;
            property.thresh     = thresh;
            
            cd(current_path);
            
        end
        
    end
    
end

