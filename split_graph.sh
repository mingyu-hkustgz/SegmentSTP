rm -r ./cmake-build-debug
mkdir cmake-build-debug
cd cmake-build-debug
cmake ..
make -j 40
cd ..


graph=NY-d.gr
save=NY-d.gr.5
./cmake-build-debug/test/test_range_stp -g ${graph} -p 5 -s ${save}

