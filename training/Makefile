TARGET := rf-train

build:

	@clang++-17 *.cpp -o $(TARGET) \
	-Iinclude \
	-O2 \
	-fopenmp \
	-fplugin=include/ClangEnzyme-17.so \
	
run:
	./$(TARGET)

clean:
	rm $(TARGET)
