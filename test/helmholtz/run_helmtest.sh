rm -rf print_testres.txt
./int2-test-hfmm2d
./int2-test-hfmm2d-hf
./int2-test-hfmm2d-vec
./int2-test-hfmm2d-mps
mv print_testres.txt ../../print_testreshelm.txt
rm -rf fort.13
