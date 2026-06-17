data = importdata("apjabbb34t4_mrt.txt");
sigma_pixels = ones(1,numel(wave))*0.66;
% converts relative velocity in km s^-1 to redshift difference
speed_of_light = 299792458;                   % speed of light                     m s⁻¹
kms_to_z = @(kms) (kms * 1000) / speed_of_light;
%%
wave_unmatched = 3600:0.01:9800;
wave = wave_unmatched(1:end-6);
abs_1_data = data.data(1:8,:);

figure(1)
flux_L2 = ones(length(wave),1);
flux_L1 = ones(length(wave),1);
for i = 1:8
flux_L1 = flux_L1.*voigt_iP(wave_unmatched,abs_1_data(i,1) + kms_to_z(abs_1_data(i,3)),10^(abs_1_data(i,4)),1,...
    (abs_1_data(i,6))*(10^5)/sqrt(2),sigma_pixels);
flux_L2 = flux_L2.*voigt_iP(wave_unmatched,abs_1_data(i,1) + kms_to_z(abs_1_data(i,3)),10^(abs_1_data(i,4)),2,...
    (abs_1_data(i,6))*(10^5)/sqrt(2),sigma_pixels);

end
plot(wave,flux_L2)
hold on 
% plot(wave,voigt_iP(wave_unmatched,0.8708,10^13.5,2,30*(10^5),sigma_pixels));
hold off


axis([5220,5250,-0.01,1.05])

REW_2796 = trapz(wave, 1-flux_L1)/(1+abs_1_data(1,1));
REW_2803 = trapz(wave, 1-flux_L2./flux_L1)/(1+abs_1_data(1,1));
REW_tot =  trapz(wave, 1-flux_L2)/(1+abs_1_data(1,1));
figure(2)
plot(wave, 1-flux_L1,'r')
hold on
plot(wave, 1-flux_L2./flux_L1,'b')
hold off
axis([5220,5250,-0.01, 1.05])


%%
abs_data = data.data;

% qsos = strings(1,length(data.data));
% 
% for i = 1:length(data.data)
%     qsos(i) = string(data.textdata{i+40});
%     % disp(qso_id(i))
% end
% 
% unqiue_qsos = sort(unqiue(qsos),"ascend");

unique_abs = nnz(abs_data(:,2) == 1);
unique_indicies = cell(1,unique_abs);
i = 1;
start_ind = 1;
for j = 2:length(abs_data)-1
    if abs_data(j,2) == 1
        unique_indicies{i} = start_ind:j-1;
        start_ind = j;
        i = i+1;
    end
end