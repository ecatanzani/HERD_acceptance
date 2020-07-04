MKDIR_P       ?=mkdir -p
COMPILER      ?=g++
TESTS         ?=true
TOP           ?=$(shell pwd)
HAVE_TERM     :=$(shell echo $$TERM)
#undef to none (linux)
ifndef HAVE_TERM
	HAVE_TERM = none
endif
#dump to none (macOS)
ifeq ($(HAVE_TERM),dumb)
	HAVE_TERM = none
endif

#dependencie
DIPS_INCLUDE = $(shell root-config --incdir)
DIPS_LIBS = $(shell root-config --ldflags) $(shell root-config --libs) 

S_DIR  = $(TOP)/source/
S_INC  = $(TOP)/include/

ANYOPT_SRC = $(TOP)/deps/AnyOption/anyoption.cpp
ANYOPT_INC = $(TOP)/deps/AnyOption

HEALPIX_CPP_INCLUDE = /HEALPiX/Healpix_3.60/include/healpix_cxx
HEALPIX_C_INCLUDE = /HEALPiX/Healpix_3.60/include
HEALPIX_CPP_LIBS_PATH = /HEALPiX/Healpix_3.60/lib
CFITSIO_LIBS_PATH = /HEALPiX/packages/cfitsio/lib

HEALPIX_CPP_LIBS = -lhealpix_cxx -lchealpix -lhealpix -lsharp

DEBUG_DIR    = Debug/obj
RELEASE_DIR  = Release/obj
DEBUG_PROG   = Debug/mapsBuilder
RELEASE_PROG = Release/mapsBuilder

# Don't compile specific files
FILTER := 

# Subdirs
SUB_DIRS := $(wildcard $(S_DIR)/**/.)\
            $(wildcard $(S_DIR)/**/**/.)
SUB_DIRS := $(subst $(S_DIR)/,,$(SUB_DIRS))

####################################################
# C FLAGS
C_FLAGS = -fPIC -D_FORCE_INLINES
# CPP FLAGS
CC_FLAGS = -std=c++14 -I$(S_INC)
# RELEASE_FLAGS
RELEASE_FLAGS = -O3
# DEBUG_FLAGS
DEBUG_FLAGS = -g -D_DEBUG -Wall -Wno-unknown-pragmas
# Linker
LDFLAGS :=
#add dips
ifneq ($(DIPS_INCLUDE),)
	CC_FLAGS += -I$(DIPS_INCLUDE) -I$(ANYOPT_INC) -I$(HEALPIX_CPP_INCLUDE) -I$(HEALPIX_C_INCLUDE)
	LDFLAGS += $(DIPS_LIBS) -lTreePlayer -L$(HEALPIX_CPP_LIBS_PATH) $(HEALPIX_CPP_LIBS) -L$(CFITSIO_LIBS_PATH) -lcfitsio
endif
####################################################
# Flags by OS
ifeq ($(shell uname -s),Linux) # LINUX
# threads
C_FLAGS += -pthread
LDFLAGS += -lpthread
endif

####################################################
# All source files
ALL_SOURCE_FILES := $(wildcard $(S_DIR)/*.cpp)\
					$(wildcard $(S_DIR)/**/*.cpp)\
					$(wildcard $(S_DIR)/**/**/*.cpp)
					
####################################################
# Object files
SOURCE_FILES = $(filter-out $(FILTER), $(ALL_SOURCE_FILES))
SOURCE_DEBUG_OBJS = $(subst $(S_DIR),$(DEBUG_DIR),$(subst .cpp,.o,$(SOURCE_FILES)))
SOURCE_RELEASE_OBJS = $(subst $(S_DIR),$(RELEASE_DIR),$(subst .cpp,.o,$(SOURCE_FILES)))

SOURCE_DEPENDENCY_FILES = $(filter-out $(FILTER), $(ANYOPT_SRC))
SOURCE_DEBUG_DEPENDENCY_OBJS = $(subst $(S_DIR),$(DEBUG_DIR),$(subst .cpp,.o,$(SOURCE_DEPENDENCY_FILES)))
SOURCE_RELEASE_DEPENDENCY_OBJS = $(subst $(S_DIR),$(RELEASE_DIR),$(subst .cpp,.o,$(SOURCE_DEPENDENCY_FILES)))

####################################################
# Output dirs
O_DEBUG_DIR    = $(TOP)/$(DEBUG_DIR)
O_RELEASE_DIR  = $(TOP)/$(RELEASE_DIR)
O_DEBUG_PROG   = $(TOP)/$(DEBUG_PROG)
O_RELEASE_PROG = $(TOP)/$(RELEASE_PROG)

##
# Support function for colored output
# Args:
#     - $(1) = Color Type
#     - $(2) = String to print
ifneq ($(HAVE_TERM),none)
define colorecho
	@tput setaf $(1)
	@echo $(2)
	@tput sgr0
endef
else
define colorecho
	@echo $(2)
endef
endif

# Color Types
COLOR_BLACK = 0
COLOR_RED = 1
COLOR_GREEN = 2
COLOR_YELLOW = 3
COLOR_BLUE = 4
COLOR_MAGENTA = 5
COLOR_CYAN = 6
COLOR_WHITE = 7

all: directories show_debug_flags debug release clean

directories: debug_make_dirs release_make_dirs

rebuild: clean directories debug release

rebuild_debug: clean_all_debug debug

rebuild_release: clean_all_release release

debug: directories show_debug_flags $(SOURCE_DEBUG_OBJS) $(SOURCE_DEBUG_DEPENDENCY_OBJS) $(TEST_SOURCE_DEBUG_OBJS)
	$(COMPILER) $(C_FLAGS) $(CC_FLAGS) $(SOURCE_DEBUG_OBJS) $(SOURCE_DEBUG_DEPENDENCY_OBJS) $(LDFLAGS) -o $(O_DEBUG_PROG)

release: directories show_release_flags $(SOURCE_RELEASE_OBJS) $(SOURCE_RELEASE_DEPENDENCY_OBJS) $(TEST_SOURCE_RELEASE_OBJS)
	$(COMPILER) $(C_FLAGS) $(CC_FLAGS) $(SOURCE_RELEASE_OBJS) $(SOURCE_RELEASE_DEPENDENCY_OBJS) $(LDFLAGS) -o $(O_RELEASE_PROG)

# makedir
debug_make_dirs:
	@${MKDIR_P} $(DEBUG_DIR)
	@for dir in $(SUB_DIRS); do \
		${MKDIR_P} $(DEBUG_DIR)/$$dir; \
	done

# makedir
release_make_dirs:
	@${MKDIR_P} $(RELEASE_DIR)
	@for dir in $(SUB_DIRS); do \
		${MKDIR_P} $(RELEASE_DIR)/$$dir; \
	done

show_debug_flags:
	$(call colorecho,$(COLOR_YELLOW),"[ Debug flags: $(C_FLAGS) $(CC_FLAGS) $(DEBUG_FLAGS) ]")

show_release_flags:
	$(call colorecho,$(COLOR_YELLOW),"[ Release flags: $(C_FLAGS) $(CC_FLAGS) $(RELEASE_FLAGS) ]")

##################################################################################################################
# DEBUG
$(SOURCE_DEBUG_OBJS):
	$(call colorecho,$(COLOR_GREEN),"[ Make debug object: $(subst $(DEBUG_DIR),,$(@:.o=.cpp)) => $(subst $(TOP)/,,$(@)) ]")
	@$(COMPILER) $(C_FLAGS) $(CC_FLAGS) $(DEBUG_FLAGS) -c $(subst $(DEBUG_DIR),$(S_DIR),$(@:.o=.cpp)) -o $@

##################################################################################################################
# RELEASE
$(SOURCE_RELEASE_OBJS):
	$(call colorecho,$(COLOR_GREEN),"[ Make release object: $(subst $(RELEASE_DIR),,$(@:.o=.cpp)) => $(subst $(TOP)/,,$(@)) ]")
	@$(COMPILER) $(C_FLAGS) $(CC_FLAGS) $(RELEASE_FLAGS) -c $(subst $(RELEASE_DIR),$(S_DIR),$(@:.o=.cpp)) -o $@

# Clean
clean: clean_debug clean_release
clean_all: clean_all_debug clean_all_release

clean_debug:
	$(call colorecho,$(COLOR_MAGENTA),"[ Delete debug obj files ]")
	@rm -f -R $(O_DEBUG_DIR)
	@rm -f $(ANYOPT_INC)/anyoption.o

clean_release:
	$(call colorecho,$(COLOR_MAGENTA),"[ Delete release obj files ]")
	@rm -f -R $(O_RELEASE_DIR)
	@rm -f $(ANYOPT_INC)/anyoption.o

clean_all_debug:
	$(call colorecho,$(COLOR_MAGENTA),"[ Delete debug obj files ]")
	@rm -f -R $(O_DEBUG_DIR)
	@rm -f $(ANYOPT_INC)/anyoption.o
	$(call colorecho,$(COLOR_MAGENTA),"[ Delete debug executable files ]")
	@rm -f -R $(O_DEBUG_PROG)

clean_all_release:
	$(call colorecho,$(COLOR_MAGENTA),"[ Delete release obj files ]")
	@rm -f -R $(O_RELEASE_DIR)
	@rm -f $(ANYOPT_INC)/anyoption.o
	$(call colorecho,$(COLOR_MAGENTA),"[ Delete release executable files ]")
	@rm -f -R $(O_RELEASE_PROG)
