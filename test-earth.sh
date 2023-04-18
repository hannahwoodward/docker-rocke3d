# Create Earth test rundeck
cd $MODELDIR/decks
make rundeck RUNSRC=E1oM20 RUN=E1oM20_Test

# Download input data
../exec/get_input_data -w E1oM20_Test /home/app/ModelE_Support/prod_input_files

# Compile & run
make clean; make -j setup RUN=E1oM20_Test
../exec/runE E1oM20_Test -cold-restart -np 2

# Create readable netcdf outputs
cd /home/app/ModelE_Support/huge_space/E1oM20_Test
scaleacc PARTIAL.accE1oM20_Test.nc aij
