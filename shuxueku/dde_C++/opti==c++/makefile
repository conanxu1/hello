# dir= /d/++codelibrary++/myddev1
# eigendir=   /d/++codelibrary++/myddev1/eigen-3.3.8

# flags = -O3 -lm -lnlopt -Wall
flags = -O3 -lm -lnlopt  -W



all: example1

.PHONY: all

dir:
	@echo "ok"
	
example1:

	@g++ $(flags)  -o example1 main.cpp pro.cpp
	@echo "example1"
clean:
	@rm example1