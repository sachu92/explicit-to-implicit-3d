# Declaration of variables
CC = g++

SRCDIR = src
BUILDDIR = build
EXEC = bin/run

SRCEXT = cpp

DEBUG_FLAG = # -g
OPTI_FLAGS = -O3 
WARN_FLAGS = -Wall 
CC_FLAGS = -I/usr/include/eigen3 -std=c++0x $(DEBUG_FLAG) $(OPTI_FLAGS) $(WARN_FLAGS)  

# File names
SOURCES = $(shell find $(SRCDIR) -type f -name *.$(SRCEXT))
OBJECTS = $(patsubst $(SRCDIR)/%,$(BUILDDIR)/%,$(SOURCES:.$(SRCEXT)=.o))

HEADERS = $(wildcard src/*.hpp)
LIBFLAGS = #-lprofiler

all: $(EXEC)

# Main target
$(EXEC): $(OBJECTS)
	@echo " Linking...";
	$(CC) $^ -o $(EXEC) $(LIBFLAGS) $(CC_FLAGS)

# To obtain object files
$(BUILDDIR)/%.o: $(SRCDIR)/%.$(SRCEXT)
	@mkdir -p $(BUILDDIR);
	$(CC) $(CC_FLAGS) -c -o $@ $< 

# To remove generated files
clean:
	@echo "Cleaning...";
	rm -r $(EXEC) $(BUILDDIR)

.PHONY: clean

