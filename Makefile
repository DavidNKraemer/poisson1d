# ######################## CONFIGURATION ########################################

# C++ compiler and compling flags
CC := g++
CFLAGS := -Wall -Wextra -g -std=c++11 -fPIC

# Directories
BUILDDIR := build
SRCDIR := solvers
DATADIR := data
TARGET := bin/main
MISC := *.txt

SRCEXT := cpp
SOURCES := $(shell find $(SRCDIR) -type f -name *.$(SRCEXT))
OBJECTS := $(patsubst $(SRCDIR)/%,$(BUILDDIR)/%,$(SOURCES:.$(SRCEXT)=.o))
DATA := $(shell find $(DATADIR) -type f -name *.txt)
INC := -I include

# linking flags, "-L/path/to/your/libs"
LDFLAGS :=
# libraries, either "/path/to/static/libs.a" or just -llibs
LIBS := -llapacke -lblas

# dependencies, your headers that are depended by objs
# note that in general, you need to group objs for different
# deps
DEPS := Complex.h

# Formatting
# Colors #
COM_COLOR := \033[0;34m
OBJ_COLOR := \033[0;36m
OK_COLOR := \033[0;32m
ERROR_COLOR := \033[0;31m
WARN_COLOR := \033[0;33m
NOCOLOR := \033[m

# Formatting #
BF := \033[1m
UL := \033[4m
NM := \033[0m
DM := \033[2m


# Strings #
OK_STRING := "[OK]"
ERROR_STRING := "[ERROR]"
WARN_STRING := "[WARNING]"
COM_STRING := "Compiling"
LNK_STRING := "Linking:"
CLN_STRING := "Cleaning:"



# suffix rule
.SUFFIXES: .cpp .cc .c++ .h .hpp .hh .h++ .o

.cpp.o:
	$(CC) $(CFLAGS) -c $< -o $@

############## RULES ##################

# Commands #
default: $(TARGET) $(OBJECTS)
	@echo "$(COM_COLOR)$(BF)Running:$(NM)$(NOCOLOR) $(TARGET)"
	@$(TARGET)

plots: $(DATA)
	@echo "plotting data"
	$(python3 scripts/plot_results.py)

$(TARGET): $(OBJECTS)
	@echo "$(COM_COLOR)$(BF)$(LNK_STRING)$(NM)$(NOCOLOR) $(OBJECTS)"
	@echo "    $(UL)$(CC) $^ $(LIBS) -o $(TARGET) $(NM)"
	@$(CC) $^ $(LIBS) -o $(TARGET) 

$(BUILDDIR)/%.o: $(SRCDIR)/%.$(SRCEXT)
	@echo "$(COM_COLOR)$(BF)$(COM_STRING)$(NM)$(NOCOLOR)"
	@mkdir -p $(BUILDDIR)
	@echo "    $(UL)$(CC) $(CFLAGS) $(INC) -c -o $@ $<$(NM)"
	@$(CC) $(CFLAGS) $(INC) -c -o $@ $<

clean:
	@echo "$(COM_COLOR)$(BF)$(CLN_STRING)$(NM)$(NOCOLOR)"
	@echo "    $(RM) -r $(BUILDDIR) $(TARGET) $(MISC)";
	@$(RM) -r $(BUILDDIR)/ $(TARGET) $(MISC);

.PHONY: clean plots
