
prob_flag = p_MgII>0.85;
flagged_N = map_N_MgIIL2(prob_flag);
flagged_z = map_z_MgIIL2(prob_flag);
flagged_sig = map_sigma_MgIIL2(prob_flag);

figure(1)
histogram(flagged_z,BinWidth=0.25)
xlabel('Absorber Redshift')
ylabel('Number Detected')
title('Histogram of Detected MgII Absorber Redshifts (probability >85%)')

figure(2)
histogram(flagged_N,BinWidth=0.5,colormap='')
xlabel('Column Density [$cm^{-2}$]', 'Interpreter', 'latex');
ylabel('Number Detected')
title('Histogram of Detected MgII Absorber Column Densities (probability >85%)')

figure(3)
histogram(flagged_sig,BinWidth=1E6)
xlabel('Sigma', 'Interpreter', 'latex');
ylabel('Number Detected')
title('Histogram of Detected MgII Absorber Sigma (probability >85%)')



% [a,b] = size(p_MgII);
% for i=1:a
% 	sliced_p_MgII = p_MgII(i,:);
% 	if sum(sliced_p_MgII>0.85) >= 1
% 		qso_flag(i) = 1;
% 	else
% 		qso_flag(i) = 0;
% 	end
% end
% qso_flag = qso_flag';
% 
% z_qso = all_zqso_dr1(logical(~filter_flags));
% 
% numel(qso_flag)