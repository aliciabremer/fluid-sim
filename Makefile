CC = gcc
CXX = g++

LIBS = -L$(lib_dir)
CXXFLAGS = -g $(INCLUDES)
INCLUDES = -I $(inc_dir) 
LDFLAGS = $(LIBS) -lglfw3 -lGL -lX11 -lpthread -lXrandr -lXi -ldl -lGLEW -lboost_system

BIN = bin
SRC_DIR = src
inc_dir = include
lib_dir = lib

TARGET = fluid-sim
CPP_SRC = $(shell find src -type f -name "*.cpp")
SRC_OBJS = $(CPP_SRC:.cpp=.o)
OBJS = $(patsubst $(SRC_DIR)/%,$(BIN)/%,$(SRC_OBJS))

all: $(TARGET)

$(BIN)/%.o: $(SRC_DIR)/%.cpp
	@mkdir -p $(@D)
	$(CXX) -c $(CXXFLAGS) -o $@ $<

$(TARGET): $(OBJS)
	$(CXX) -o $@ $^ $(LDFLAGS)

.PHONY : clean
clean :
	$(RM) $(OBJS)
	$(RM) $(TARGET)