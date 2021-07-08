rm -rf print_testres.txt
./int2-test-cfmm2d
./int2-test-cfmm2d-vec
./int2-test-lfmm2d
./int2-test-lfmm2d-vec
./int2-test-rfmm2d
./int2-test-rfmm2d-vec
mv print_testres.txt ../../print_testreslap.txt
rm -rf fort.13
