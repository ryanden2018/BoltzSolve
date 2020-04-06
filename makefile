ConvDiff2d: ConvDiff2d.cpp
	emcc ConvDiff2d.cpp -O3 \
	-I /home/ryan/Downloads/eigen-3.3.7 \
	-s ALLOW_MEMORY_GROWTH=1 \
	-o ./out/ConvDiff2d.html \
	--shell-file ./html_template/shell_minimal.html

