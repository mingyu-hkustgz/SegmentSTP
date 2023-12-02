rm -r ./cmake-build-debug
mkdir cmake-build-debug
cd cmake-build-debug
cmake ..
make -j 40
cd ..


graph=NY-d.gr
index=NY.index
method=2
length=10
dataset=NY
res_path=./res.log
./cmake-build-debug/test/test_build_range_index -g ${graph} -n ${dataset} -i ${index} -l ${length} -m ${method} -s ${res_path}

