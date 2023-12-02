rm -r ./cmake-build-debug
mkdir cmake-build-debug
cd cmake-build-debug
cmake ..
make -j 40
cd ..

length=10

graph=NY-d.gr
index=NY.index
method=1
length=10
dataset=NY
res_path=./res.log

./cmake-build-debug/test/test_index_search_time -g ${graph} -n ${dataset} -i ${index} -l ${length} -m ${method} -s ${res_path}