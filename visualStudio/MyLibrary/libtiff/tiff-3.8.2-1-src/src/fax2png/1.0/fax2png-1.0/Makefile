CFLAGS=-D__GW32__ -D_LARGEFILE_SOURCE=1 -D_LARGEFILE64_SOURCE=1 -D_FILE_OFFSET_BITS=64 \
	-Wall -O3 -fms-extensions -mms-bitfields  -fno-exceptions -fomit-frame-pointer \
	-march=i386 -mcpu=i686 \
	-idirafter D:/progra~1/gnuwin32/include -idirafter G:/glibc/include
LDFLAGS = -Wl,-s -Wl,--force-exe-suffix -Wl,--enable-auto-import -Wl,--enable-runtime-pseudo-reloc \
	-LD:/progra~1/gnuwin32/lib -lgw32c
EXEEXT = .exe
	
.SUFFIXES:
.SUFFIXES: .c .o .exe .rc

# Rule to make compiled resource (Windows)
.rc.o:
	windres --include-dir $(resdir) -i $< -o $@

fax2png$(EXEEXT): fax2png.o fax2png-res.o
	gcc fax2png.o  fax2png-res.o -o $@ $(LDFLAGS) -lpng -lz -ltiff 

install: fax2png$(EXEEXT)
	/bin/install.exe -cpD $< c:/progra~1/Fax2Png/bin/$<

clean:
	rm -f fax2png$(EXEEXT) fax2png.o

