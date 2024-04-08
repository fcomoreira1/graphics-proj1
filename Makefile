CXX=g++
CXXFLAGS=-std=c++17 -O3

EXE=render
ODIR=build

SRCS = $(wildcard *.cpp)
DEPS = $(wildcard *.h) $(wildcard ./scenes/*.h) 
_OBJS = $(subst .cpp,.o,$(SRCS))
OBJS = $(patsubst %,$(ODIR)/%,$(_OBJS))

$(ODIR)/%.o: %.cpp $(DEPS)
	$(CXX) -c -o $@ $< $(CXXFLAGS)

render: $(SRCS) $(DEPS)
	$(CXX) -o render $(SRCS) $(CXXFLAGS) $(LIBS)

.PHONY: incremental
incremental: $(OBJS)
	$(CXX) -o render $^ $(CXXFLAGS) $(LIBS)

.PHONY: clean
clean:
	rm -r $(OBJS) $(EXE)
