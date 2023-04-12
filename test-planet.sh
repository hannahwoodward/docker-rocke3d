cd /home/app/modelE2_planet_1.0/decks
make rundeck RUNSRC=P1SoM40 RUN=P1SoM40_Test
../exec/get_input_data -w P1SoM40_Test /home/app/ModelE_Support/prod_input_files

sed -i "s|!\{0,1\} \{0,1\}solar_spec_dir=.*|solar_spec_dir='/home/app/ModelE_Support/socrates/stellar_spectra'|" P1SoM40_Test.R
sed -i "s|!\{0,1\} \{0,1\}spectral_dir=.*|spectral_dir='/home/app/ModelE_Support/socrates/spectral_files'|" P1SoM40_Test.R
sed -i "s|!\{0,1\} \{0,1\}solar_spec=.*|solar_spec='sun'|" P1SoM40_Test.R
sed -i "s|!\{0,1\} \{0,1\}spectral_file_lw=.*|spectral_file_lw='sp_lw_ga7/sp_lw_ga7_dsa'|" P1SoM40_Test.R
sed -i "s|!\{0,1\} \{0,1\}spectral_file_sw=.*|spectral_file_sw='sp_sw_ga7/sp_sw_ga7_dsa'|" P1SoM40_Test.R

make clean; make -j setup RUN=P1SoM40_Test
../exec/runE P1SoM40_Test -cold-restart -np 2
