close all force
clear all
clc
format compact
format long

desired_relative_radial_shift=["0p5"]; % "1p0" "1p5" "2p0" "2p5" "3p0" "3p5" "4p0" "4p5" "5p0" "5p5" "6p0" "6p5" "7p0" "7p5" "8p0" "8p5" "9p0" "9p5" "10p0" "10p5" "11p0" "11p5" "12p0" "12p5" "13p0" "13p5" "14p0" "14p5" "15p0"];

import com.comsol.model.*
import com.comsol.model.util.*

model=mphload('step_2_compute_Etot_achieved.mph');

wb=waitbar(0,'Please wait...');
for slucaj=1:length(desired_relative_radial_shift)
    apap=load(strcat("./SLM_8bit_from_Arturo/SLM_Einc_lens_R_4mu_rad_shift_",desired_relative_radial_shift(slucaj),"lam0_theta_0deg.txt"));
    aaaa=load(strcat("./Einc/abs_Einc_lens_R_4mu_rad_shift_",desired_relative_radial_shift(slucaj), "lam0_theta_0deg.txt"));
    % aaaa=load(strcat("./Einc/abs_Einc_DIGITIZED_lens_R_4mu_rad_shift_",desired_relative_radial_shift(slucaj), "lam0_theta_0deg.txt"));
    % pppp=load(strcat("./Einc/pha_Einc_DIGITIZED_lens_R_4mu_rad_shift_",desired_relative_radial_shift(slucaj), "lam0_theta_0deg.txt"));
    xs=aaaa(:,1);
    clear aaaa;
    amp=apap(:,2); % amp=aaaa(:,2);
    pha=apap(:,3); % pppp(:,2);
    model.func('int15').set('table', cellfun( @(x) num2str(x, '%0.10g'), num2cell( [xs, amp] ), 'UniformOutput', false));
    model.func('int16').set('table', cellfun( @(x) num2str(x, '%0.10g'), num2cell( [xs, pha] ), 'UniformOutput', false));
    model.sol('sol1').runAll;
    pd= mpheval(model,'emw.normE');
    save(strcat("./Etot_achieved/abs_Etot_achieved_lens_R_4mu_desired_rad_shift_",desired_relative_radial_shift(slucaj),"lam0_desired_theta_0deg_SLM_Einc.mat"),'-struct','pd');
    % E_tot_near(slucaj,:)=pd.d1;
    % savefile=strcat("./Etot_achieved/abs_Etot_achieved_lens_R_4mu_desired_rad_shift_",desired_relative_radial_shift(slucaj),"lam0_desired_theta_0deg_DIGITIZED_Einc.txt");
    % snimiti=pd; % E_tot_near(slucaj,:);
    % save(savefile,'snimiti'); % ,'-ascii');
    waitbar(slucaj/length(desired_relative_radial_shift),wb);
end
close(wb);