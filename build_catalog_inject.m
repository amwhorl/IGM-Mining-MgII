% Build catalogs usable for spectra from DESI dr1 and from SALSA simulated
% absorbers

DESI_dr1q = ...
fitsread('/home/cosmic/desi-env/sim_abs_mgii/QSO_catalog.fits', 'binarytable');
all_tid_dr1               = DESI_dr1q{1};
all_tid_dr1               = all_tid_dr1(null_search_filter);
all_tid_dr1               = all_tid_dr1(null_QSO_ind);
all_zqso_dr1              = DESI_dr1q{2};
all_zqso_dr1              = all_zqso_dr1(null_search_filter);
all_zqso_dr1              = all_zqso_dr1(null_QSO_ind);
all_RA_dr1                = DESI_dr1q{3};
all_RA_dr1                = all_RA_dr1(null_search_filter);
all_RA_dr1                = all_RA_dr1(null_QSO_ind);
all_DEC_dr1               = DESI_dr1q{4};
all_DEC_dr1               = all_DEC_dr1(null_search_filter);
all_DEC_dr1               = all_DEC_dr1(null_QSO_ind);

% Load all 5 simulated absorber catalogs 
sim_cat_05 = ... % z = 0.5
fitsread('/home/cosmic/desi-env/sim_abs_mgii/MgII_Simulated_Catalog_z0.5.fits', 'binarytable');
sim_id_05                 = sim_cat_05{1};
sim_z_05                  = sim_cat_05{2};
sim_EW2796_05             = sim_cat_05{3};
sim_EW2803_05             = sim_cat_05{4};
sim_N2796_05              = sim_cat_05{5};
sim_N2803_05              = sim_cat_05{6};

sim_cat_07 = ... % z = 0.7
fitsread('/home/cosmic/desi-env/sim_abs_mgii/MgII_Simulated_Catalog_z0.7.fits', 'binarytable'); 
sim_id_07                 = sim_cat_07{1};
sim_z_07                  = sim_cat_07{2};
sim_EW2796_07             = sim_cat_07{3};
sim_EW2803_07             = sim_cat_07{4};
sim_N2796_07              = sim_cat_07{5};
sim_N2803_07              = sim_cat_07{6};

sim_cat_10 = ... % z = 1.0
fitsread('/home/cosmic/desi-env/sim_abs_mgii/MgII_Simulated_Catalog_z1.0.fits', 'binarytable'); 
sim_id_10                 = sim_cat_10{1};
sim_z_10                  = sim_cat_10{2};
sim_EW2796_10             = sim_cat_10{3};
sim_EW2803_10             = sim_cat_10{4};
sim_N2796_10              = sim_cat_10{5};
sim_N2803_10              = sim_cat_10{6};

sim_cat_15 = ... % z = 1.5
fitsread('/home/cosmic/desi-env/sim_abs_mgii/MgII_Simulated_Catalog_z1.5.fits', 'binarytable'); 
sim_id_15                 = sim_cat_15{1};
sim_z_15                  = sim_cat_15{2};
sim_EW2796_15             = sim_cat_15{3};
sim_EW2803_15             = sim_cat_15{4};
sim_N2796_15              = sim_cat_15{5};
sim_N2803_15              = sim_cat_15{6};

sim_cat_20 = ... % z = 2.0
fitsread('/home/cosmic/desi-env/sim_abs_mgii/MgII_Simulated_Catalog_z2.0.fits', 'binarytable'); 
sim_id_20                 = sim_cat_20{1};
sim_z_20                  = sim_cat_20{2};
sim_EW2796_20             = sim_cat_20{3};
sim_EW2803_20             = sim_cat_20{4};
sim_N2796_20              = sim_cat_20{5};
sim_N2803_20              = sim_cat_20{6};





% save catalog 
variables_to_save = {'all_tid_dr1', 'all_RA_dr1', 'all_DEC_dr1', 'all_zqso_dr1'};
save(sprintf('%s/catalog', processed_directory(releaseTest)), ...
    variables_to_save{:}, '-v7.3');

variables_to_save = {
    'sim_id_05','sim_z_05','sim_EW2796_05','sim_EW2803_05','sim_N2796_05','sim_N2803_05',...
    'sim_id_07','sim_z_07','sim_EW2796_07','sim_EW2803_07','sim_N2796_07','sim_N2803_07',...
    'sim_id_10','sim_z_10','sim_EW2796_10','sim_EW2803_10','sim_N2796_10','sim_N2803_10',...
    'sim_id_15','sim_z_15','sim_EW2796_15','sim_EW2803_15','sim_N2796_15','sim_N2803_15',...
    'sim_id_20','sim_z_20','sim_EW2796_20','sim_EW2803_20','sim_N2796_20','sim_N2803_20'};
save(sprintf('%s/SALSA_catalog', processed_directory(releaseTest)), ...
    variables_to_save{:}, '-v7.3')