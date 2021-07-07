rm -rf print_testres.txt
./int2-test-hfmm2d
./int2-test-hfmm2d-vec
mv print_testres.txt ../../print_testreshelm.txt
rm -rf fort.13
