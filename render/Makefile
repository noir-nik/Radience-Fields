CC := g++

BUILD_DIR := build
TARGET := rf-render

RENDER_DIR := output

INCLUDES := -Iinclude -Iinclude/glfw/include

CXXFLAGS := -MP -MMD -fopenmp -O2 $(INCLUDES)

SRCS := $(wildcard *.cpp)
GLAD_SRC := include/glad/glad.c
OBJS := $(patsubst %.cpp, $(BUILD_DIR)/%.o, $(notdir $(SRCS))) $(BUILD_DIR)/glad.o

SRC_GLFW := include/glfw/src
OBJS_GLFW  := $(patsubst %.c, $(BUILD_DIR)/glfw/%.o, $(notdir $(wildcard $(SRC_GLFW)/*.c)))

LDFLAGS :=
ifeq ($(OS), Windows_NT)
LDFLAGS += -lgdi32
GLFW_PLATFORM := _GLFW_WIN32
else
LDFLAGS += -lGL -lm
GLFW_PLATFORM := _GLFW_X11
endif

all: build

build: create_dir $(TARGET)

create_dir:
	@mkdir -p $(BUILD_DIR)
	@mkdir -p $(BUILD_DIR)/glfw
	@mkdir -p $(RENDER_DIR)


$(TARGET): $(OBJS_GLFW) $(OBJS)
	@echo "Linking $(TARGET)"
	@$(CC) $(CXXFLAGS)  $(INCLUDES) $(OBJS) $(OBJS_GLFW) $(LDFLAGS) -o $(TARGET)
	@echo "Build complete for $(TARGET)"

$(BUILD_DIR)/glfw/%.o: $(SRC_GLFW)/%.c
	@echo "Compiling $(notdir $<)"
	@$(CC) $(CXXFLAGS) -x c $ -c $< -o$@ -D$(GLFW_PLATFORM)

$(BUILD_DIR)/%.o: %.cpp
	@echo "Compiling $(notdir $<)"
	@$(CC) $(CXXFLAGS) -o $@ -c $<

$(BUILD_DIR)/%.o: include/glad/%.c
	@echo "Compiling $(notdir $<)"
	@$(CC) $(CXXFLAGS) -x c -o $@ -c $<

run:
	@$(TARGET) -m data/model.dat -s 128

clean:
	@rm -f $(TARGET) $(wildcard $(BUILD_DIR)/*.o) $(wildcard $(BUILD_DIR)/glfw/*.o) $(DEPFILES)