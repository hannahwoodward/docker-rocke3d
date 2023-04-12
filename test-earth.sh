cd /home/app/modelE2_planet_1.0/decks
make rundeck RUNSRC=E1oM20 RUN=E1oM20_Test
../exec/get_input_data -w E1oM20_Test /home/app/ModelE_Support/prod_input_files
make clean; make -j setup RUN=E1oM20_Test
../exec/runE E1oM20_Test -cold-restart -np 2
