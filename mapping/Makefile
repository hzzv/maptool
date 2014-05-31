# executable name
TARGET=maptool

# directories
SRCDIR=./
INCLUDE=include
OBJDIR=obj
BINDIR=./

# sources
SOURCES = $(wildcard $(SRCDIR)/*.cpp)
INCLUDES= $(wildcard $(INCLUDE)/*.h)
OBJECTS = $(SOURCES:$(SRCDIR)/%.cpp=$(OBJDIR)/%.o)

# compiler and flags

CXX=g++
CXXFLAGS= -O2 -std=gnu++11 -static
RM=rm
WFLAGS=-Wall -Wextra -Wno-unused-result 
#-g -pg
MYLIBS=-lbgzf -lz -L../../ocaml-bgzf

# commands

all: $(BINDIR)/$(TARGET)

$(BINDIR)/$(TARGET): $(OBJECTS)
	$(CXX) -o $@ $(WFLAGS) $(OBJECTS) $(MYLIBS)
	@echo "Done."

$(OBJECTS): $(OBJDIR)/%.o : $(SRCDIR)/%.cpp
	$(CXX) $(CXXFLAGS) $(WFLAGS) -c $< -o $@
	@echo "Compiled "$<" successfully."

clean:
	$(RM) $(OBJECTS)
	@echo "Clean complete."

remove: clean
	$(RM) $(BINDIR)/$(TARGET)
	@echo "Remove complete."

$(shell   mkdir -p $(OBJDIR)) 