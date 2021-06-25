%Correct (smoothed signal) vs. Online signal
load kc125.mat
thresh = 0.95;
correct_sig = XDGRASP_res_sig(kc,thresh);
online_sig = kc2Res_Sig(kc,150,thresh);

figure, plot(1-correct_sig); %notice the shift of the sign
hold on, plot(online_sig), title('Comparison of the correct res. signal with online res. signal'), xlabel('Spokes'), ylabel('respiratory signal'),
legend('correct res. signal','online res. signal')


%Correct (smoothed signal) vs. Online signal smoothed after the full
%acquisition (offline smoothing)
smooth_sig = online_sig;
smooth_sig(150:end) = smooth(online_sig(150:end), 6, 'lowess');%smoothing the non-smooth part

figure, plot(1-correct_sig), hold on ,plot(smooth_sig), title('Comparision of the correct res. signal with the online res. signal smoothed retrospectively (offline smoothing)'),
xlabel('Spokes'), ylabel('respiratory signal'), legend('correct res. signal','online res. signal smoothed')


%Correct signal vs. online smoothed online signal (online smoothing)
online_smooth_sig = kc2onlineSmoothing(kc, 150);
figure, plot(correct_sig), hold on ,plot(online_smooth_sig), title('Comparision of the correct res. signal with the online res. signal with online smoothing'),
xlabel('Spokes'), ylabel('respiratory signal'), legend('correct res. signal','online res. signal smoothed')

%Changing threshold
compare_thresh = [0.05:0.3:0.95];     %different values of threshold to compare
signal_legend = strsplit(num2str(compare_thresh),' ');
signals_thresh = zeros(length(online_sig), length(compare_thresh));
figure
signal1 = zeros(length(online_sig), 1);
for ii = 1:length(compare_thresh)
   signals_thresh(:,ii) =  XDGRASP_res_sig(kc,compare_thresh(ii));
   if ii==1
       signal1 = signals_thresh(:,ii);
   end
   if norm(signals_thresh(:,ii) + signal1 -1) > norm(signals_thresh(:,ii)-signal1)     %finding the orientation/sign of the res. signal
       Sign = 1;
   else 
       Sign = -1;
   end
   
   plot(Sign*signals_thresh(:,ii)+ 1*ii + (1-Sign)/2), hold on       %plotting each signal with an offset: +1*ii
   
end
title('Comparision of the res. signal with different threshold values'),
xlabel('Spokes'), ylabel('respiratory signal'), legend(signal_legend)



