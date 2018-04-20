LIBDIR=../lib
CXXFLAGS= -O -I.

target=$(LIBDIR)/libmodel_sm_gqlnullgq.a

all: $(target)

objects=HelAmps_sm_gqlnullgq.o Parameters_sm_gqlnullgq.o rambo.o ../../PROC_SA_CPP_sm_4/src/read_slha.o

$(target): $(objects)
	$(AR) cru $@ $(objects)
	ranlib $(target)

.PHONY: clean

clean:
	rm -f $(target)
	rm -f $(objects)

