function signal = kc2onlineSmoothing(kc, initial) 
%online smoothing - how to smooth a ftion if you don't know how it'll
%evolve? -> lets try to hardcode it from the initial sample (I'll suppose a periodic ftion)

%% getting initial data
thresh = 0.95;
Res_Signal = XDGRASP_res_sig(kc,thresh); 
Sample = Res_Signal(1:initial);
Res_Signal_sorts = zeros(length(Sample),1);
rising_falling = zeros(length(Sample),1);

%sig = 0.0488;
%sort the data into groups of similar signal values
sorts = 20; %I am supposing roughly evenly spaced data
sample_max = max(Sample);
sample_min = min(Sample);
gap = (sample_max-sample_min)/sorts;
borders = [sample_min:gap:sample_max];%evenly spaced data sorts
follow_average = zeros(length(borders),2); %for each data sort

% sorting all the init data
avrg_1 = zeros(length(borders), 2);
avrg_2 = zeros(length(borders), 2);
avrg_3 = zeros(length(borders), 2);
avrg_4 = zeros(length(borders), 2);
avrg_5 = zeros(length(borders), 2);
for ii = 1:(length(borders)-1)
    Indices = find(borders(ii)<Sample(1:end-5) & Sample(1:end-5)<borders(ii+1));
    avrg_r = 0;
    avrg_f = 0;
    
    kk = 0;
    for jj = 1:length(Indices)
        Res_Signal_sorts(Indices(jj)) = ii; 
       
        if length(Sample) >= Indices(jj)+3
            if Indices(jj)>1 & Sample(Indices(jj)-1)>Sample(Indices(jj))
                rising_falling(Indices(jj)) = 0; %it is falling
                avrg_f = avrg_f + (Sample(Indices(jj)+1) + Sample(Indices(jj)+2) + Sample(Indices(jj)+3))/3;
                
                avrg_1(ii,2) = avrg_1(ii,2) + Sample(Indices(jj)+1);
                avrg_2(ii,2) = avrg_2(ii,2) + Sample(Indices(jj)+2);
                avrg_3(ii,2) = avrg_3(ii,2) + Sample(Indices(jj)+3);
                avrg_4(ii,2) = avrg_4(ii,2) + Sample(Indices(jj)+4);
                avrg_5(ii,2) = avrg_5(ii,2) + Sample(Indices(jj)+5);
                kk = kk +1;
            else
                rising_falling(Indices(jj)) = 1; %else: it must be rising
                avrg_r = avrg_r + (Sample(Indices(jj)+1) + Sample(Indices(jj)+2) + Sample(Indices(jj)+3))/3;
                
                avrg_1(ii,1) = avrg_1(ii,1) + Sample(Indices(jj)+1);
                avrg_2(ii,1) = avrg_2(ii,1) + Sample(Indices(jj)+2);
                avrg_3(ii,1) = avrg_3(ii,1) + Sample(Indices(jj)+3);
                avrg_4(ii,1) = avrg_4(ii,1) + Sample(Indices(jj)+4);
                avrg_5(ii,1) = avrg_5(ii,1) + Sample(Indices(jj)+5);
            end  
        end
        
    end
    if jj-kk ~= 0
        follow_average(ii,1) = avrg_r/(jj-kk);   %first column is follow average for rising number of this sort
        avrg_1(ii,1) = avrg_1(ii,1)/(jj-kk);
        avrg_2(ii,1) = avrg_2(ii,1)/(jj-kk);
        avrg_3(ii,1) = avrg_3(ii,1)/(jj-kk);
        avrg_4(ii,1) = avrg_4(ii,1)/(jj-kk);
        avrg_5(ii,1) = avrg_5(ii,1)/(jj-kk);
    end
    if kk ~= 0
        follow_average(ii,2) = avrg_f/kk ;       %2nd column is follow average for a falling signal 
        avrg_1(ii,2) = avrg_1(ii,2)/kk;
        avrg_2(ii,2) = avrg_2(ii,2)/kk;
        avrg_3(ii,2) = avrg_3(ii,2)/kk;
        avrg_4(ii,2) = avrg_4(ii,2)/kk;
        avrg_5(ii,2) = avrg_5(ii,2)/kk;
    end
end

%classify the new data point into correct sort and rising/falling
Res_Signal_smooth = zeros(length(Res_Signal),1);
Res_Signal_smooth(1:initial) = Res_Signal(1:initial);
for rr = initial:length(Res_Signal)
    sig = Res_Signal(rr);
    %unorder = [sig borders]; %sig is new data point
    %[~,order] = sort(unorder);
    %sig_sort = find(order==1)%we've classified our input into the corresponding sort
    
    sig_sort = floor(sig/gap);
    if sig_sort <1
        sig_sort = 1;
    elseif sig_sort>length(borders)
        sig_sort = length(borders);
    end    
    
    if Res_Signal(rr-1)>sig
        rising_falling(rr) = 0;%it is falling
        signal_segment = [Res_Signal(rr-5:rr)', avrg_1(sig_sort,2), avrg_2(sig_sort,2),avrg_3(sig_sort,2), avrg_4(sig_sort,2), avrg_5(sig_sort,2)];
    else
        rising_falling(rr) = 1;%it must be rising
        signal_segment = [Res_Signal(rr-5:rr)', avrg_1(sig_sort,1), avrg_2(sig_sort,1),avrg_3(sig_sort,1), avrg_4(sig_sort,1), avrg_5(sig_sort,1)];
    end    

    %output the smoothing correction to the new number
    signal_segment_smooth = smooth(signal_segment, 6, 'lowess');  
    corr_sig = signal_segment_smooth(6); %select the middle number
    
    Res_Signal_smooth(rr) = corr_sig;
    
end
signal = Res_Signal_smooth;
end