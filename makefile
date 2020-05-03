BoltzSolve: BoltzSolve.cpp
	mkdir -p out
	emcc BoltzSolve.cpp -O3 \
	-s ALLOW_MEMORY_GROWTH=1 \
	-s NO_EXIT_RUNTIME=1  \
	-s "EXTRA_EXPORTED_RUNTIME_METHODS=['ccall']" \
	-o ./out/BoltzSolve.html \
	--shell-file ./html_template/shell_minimal.html
	