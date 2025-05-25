include ${FSLCONFDIR}/default.mk

PROJNAME = fabber_pet
XFILES   = fabber_pet
SOFILES  = libfsl-fabber_models_pet.so
AFILES   = libfabber_models_pet.a

# The FSL build system changed
# substantially in FSL 6.0.6
# FSL >= 6.0.6
ifeq (${FSL_GE_606}, true)
  LIBS = -lfsl-fabberexec -lfsl-fabbercore -lfsl-newimage \
         -lfsl-miscmaths -lfsl-cprob -lfsl-utils \
         -lfsl-NewNifti -lfsl-znz -lz -ldl
# FSL <= 6.0.5
else
  ifeq ($(shell uname -s), Linux)
    MATLIB := -lopenblas
  endif
  USRINCFLAGS = -I${INC_NEWMAT} -I${INC_CPROB} -I${INC_BOOST} \
                -I.. -I${FSLDIR}/extras/include/armawrap
  USRLDFLAGS  = -L${LIB_NEWMAT} -L${LIB_CPROB} -L../fabber_core  \
                -lfabberexec -lfabbercore -lutils -lnewimage     \
                -lmiscmaths -lcprob ${MATLIB} -lNewNifti -lznz -lz -ldl
endif

# Forward models
OBJS =  fwdmodel_pet.o fwdmodel_pet_1TCM.o fwdmodel_pet_2TCM.o fwdmodel_pet_2TCM_IR.o

# For debugging:
#OPTFLAGS = -ggdb

# Pass Git revision details
GIT_SHA1 := $(shell git describe --match=NeVeRmAtCh --always --abbrev=40 --dirty)
GIT_DATE := $(shell git log -1 --format=%ad --date=local)
CXXFLAGS += -DGIT_SHA1=\"${GIT_SHA1}\" -DGIT_DATE="\"${GIT_DATE}\""

#
# Build
#

# FSL >=606 uses dynamic linking
ifeq (${FSL_GE_606}, true)
all: ${XFILES} ${SOFILES}

# models in a library
libfsl-fabber_models_pet.so : ${OBJS}
	${CXX} ${CXXFLAGS} -shared -o $@ $^ ${LDFLAGS}

# fabber built from the FSL fabbercore library including the models specifieid in this project
fabber_pet : fabber_client.o | libfsl-fabber_models_pet.so
	${CXX} ${CXXFLAGS} -o $@ $< -lfsl-fabber_models_pet ${LDFLAGS}

# FSL <=605 uses static linking
else
all: ${XFILES} ${AFILES}

libfabber_models_pet.a : ${OBJS}
	${AR} -r $@ $^

fabber_pet : fabber_client.o ${OBJS}
	${CXX} ${CXXFLAGS} -o $@ $^ ${LDFLAGS}
endif
# DO NOT DELETE
