ConvDiff2d: ConvDiff2dDense.cpp
	emcc ConvDiff2dDense.cpp -O3 \
	-I /home/ryan/Downloads/eigen-3.3.7 \
	-o ./out/ConvDiff2d.html

