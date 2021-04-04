dir= /d/++codelibrary++/myddev1
eigendir=   /d/++codelibrary++/myddev1/eigen-3.3.8
flags = -O3 -std=c++17 -mavx


all: example1
.PHONY: all

dir:
	@echo "ok"
	
example1:
	@g++ $(flags)  -o example1 -I$(dir)   -I$(eigendir) mydde_tool.cpp ddetool.cpp
	@echo "example1"
clean:
	@rm example1