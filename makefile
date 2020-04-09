ConvDiff2d: ConvDiff2d.cpp
	mkdir -p out
	cp html_template/*.png out
	emcc ConvDiff2d.cpp -O3 \
	-I /home/ryan/Downloads/eigen-3.3.7 \
	-s ALLOW_MEMORY_GROWTH=1 \
	-o ./out/ConvDiff2d.html \
	--shell-file ./html_template/shell_minimal.html
	emcc ConvDiff2d.cpp -O3 \
	-I /home/ryan/Downloads/eigen-3.3.7 \
	-s ALLOW_MEMORY_GROWTH=1 \
	-s WASM=0 \
	-o ./out/ConvDiff2dJS.html \
	--shell-file ./html_template/shell_minimalJS.html

