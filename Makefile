CXX = g++
CXXFLAGS += -Wall -Wextra -MMD -std=c++11 -O3
LDFLAGS += -lm
TARGET = gsmf
SRCIGNORE =
SRC_DIR = src
SRCS = $(filter-out $(SRC_DIR)/$(SRCIGNORE), $(wildcard $(SRC_DIR)/*.cc))
OBJS = $(SRCS:.cc=.o)
DEPS = $(OBJS:.o=.d)
INI_DIR = ~/template
INIT = for.cc
RM = rm -fv
MODE = release
ifeq ($(MODE),debug)
	CXXFLAGS += -ggdb -DNDEBUG
	CXXFLAGS := $(filter-out -O3, $(CXXFLAGS))
endif
ifeq ($(MODE),profile)
	CXXFLAGS += -pg -fno-inline-functions-called-once -fno-optimize-sibling-calls
endif

.PHONY:all information
all: information $(TARGET)
information:
ifneq ($(MODE),release)
ifneq ($(MODE),debug)
ifneq ($(MODE),profile)
	@echo "Invalid build mode."
	@exit 1
endif
endif
endif
	@echo "Building $(TARGET) on "$(MODE)" MODE"
	@echo "Flags:" $(CXXFLAGS)
	@echo "............................."

$(TARGET): $(OBJS) Makefile
	@echo linking $(OBJS)
	@$(CXX) $(CXXFLAGS) -o $@ $(OBJS) $(LDFLAGS)

-include $(DEPS)
%.o:%.cc Makefile
	@echo compiling $<
	@$(CXX) -c $(CXXFLAGS) -o $@ $<

.SILENT:debug profile clean init d prof p c i
debug d:
	$(MAKE) MODE=debug   --always-make --no-print-directory
profile prof p:
	$(MAKE) MODE=profile --always-make --no-print-directory
clean c:
	$(RM) $(TARGET) $(OBJS) $(DEPS) gmon.out .*.swp .*.swo *~
init i:
	-@if [ -e $(SRC_DIR)/main.cc ] ;then \
	echo main.cc already exists. && echo remove main.cc and try again. ;\
	else \
	if [ ! -e $(SRC_DIR) ] ; then mkdir $(SRC_DIR) ; fi && \
	cp -v $(INI_DIR)/$(INIT) $(SRC_DIR)/main.cc; fi
