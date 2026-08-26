SHELL := /bin/sh
ROOT := $(abspath .)
VENV := $(ROOT)/.venv
PYTHON := $(VENV)/bin/python
BIN := $(ROOT)/bin/
EXTERN_LIBS := $(ROOT)/extern/libs
BUILD := $(ROOT)/.build

FFTW_LIBS := \
	$(EXTERN_LIBS)/libfftw3.a \
	$(EXTERN_LIBS)/libfftw3_threads.a \
	$(EXTERN_LIBS)/libfftw3f.a \
	$(EXTERN_LIBS)/libfftw3f_threads.a
LBFGS_LIB := $(EXTERN_LIBS)/liblbfgs.a

.PHONY: setup venv extern ensure-fftw ensure-lbfgs build-fftw build-lbfgs \
	fuda native spinunidec spintools catia portable-runtime bytecode test doctor \
	clean distclean rebuild-extern app dmg installer

# Deliberately express the build as a dependency graph.  This makes `make -jN`
# safe: doctor cannot run until native compilation and runtime bundling
# have completed. Python application launchers are permanent files in bin/.
setup: doctor

# Assemble a double-clickable macOS application from the completed runtime.
# External packages (NMRPipe, MDDNMR, csh) are deliberately not bundled.
app: doctor
	$(PYTHON) scripts/assemble_macos_app.py

# Create a drag-to-Applications disk image for distribution.
# `installer` is an alias for the same target.
dmg: app
	$(PYTHON) scripts/create_macos_dmg.py

installer: dmg

doctor: portable-runtime bytecode
	$(PYTHON) scripts/doctor.py

# Precompile application bytecode as part of every completed build. Shipping
# checked-hash .pyc files avoids a first-launch compile penalty and remains
# valid when application packaging changes source file timestamps.
bytecode: venv
	$(PYTHON) -m compileall -q -f -j 0 --invalidation-mode checked-hash applications scripts

portable-runtime: native
	$(PYTHON) scripts/bundle_macos_runtime.py

native: fuda spinunidec spintools catia

venv:
	./scripts/bootstrap.py

extern: ensure-fftw ensure-lbfgs

ensure-fftw: | $(EXTERN_LIBS)
	@if [ -f "$(EXTERN_LIBS)/libfftw3.a" ] && \
	    [ -f "$(EXTERN_LIBS)/libfftw3_threads.a" ] && \
	    [ -f "$(EXTERN_LIBS)/libfftw3f.a" ] && \
	    [ -f "$(EXTERN_LIBS)/libfftw3f_threads.a" ]; then \
		echo "FFTW: already installed in $(EXTERN_LIBS) (4/4 libraries)"; \
	else \
		$(MAKE) --no-print-directory build-fftw; \
	fi

ensure-lbfgs: | $(EXTERN_LIBS)
	@if [ -f "$(LBFGS_LIB)" ]; then \
		echo "libLBFGS: already installed in $(EXTERN_LIBS)"; \
	else \
		$(MAKE) --no-print-directory build-lbfgs; \
	fi

# Build only the four FFTW archives used by the suite.  Do not build FFTW's
# tests, benchmark programs, MPI helpers, documentation or wisdom tools.
build-fftw: | $(EXTERN_LIBS)
	@echo "Building FFTW static libraries only"
	@if [ -f "$(ROOT)/extern/fftw-3.3.10/config.status" ]; then \
		$(MAKE) -C "$(ROOT)/extern/fftw-3.3.10" distclean >/dev/null 2>&1 || true; \
	fi
	rm -rf "$(BUILD)/fftw-double" "$(BUILD)/fftw-float"
	mkdir -p "$(BUILD)/fftw-double" "$(BUILD)/fftw-float"
	cd "$(BUILD)/fftw-double" && "$(ROOT)/extern/fftw-3.3.10/configure" --disable-shared --enable-static --enable-threads
	$(MAKE) -C "$(BUILD)/fftw-double/kernel" all
	$(MAKE) -C "$(BUILD)/fftw-double/simd-support" all
	$(MAKE) -C "$(BUILD)/fftw-double/dft" all
	$(MAKE) -C "$(BUILD)/fftw-double/rdft" all
	$(MAKE) -C "$(BUILD)/fftw-double/reodft" all
	$(MAKE) -C "$(BUILD)/fftw-double/api" all
	$(MAKE) -C "$(BUILD)/fftw-double" libfftw3.la
	$(MAKE) -C "$(BUILD)/fftw-double/threads" libfftw3_threads.la
	cp "$(BUILD)/fftw-double/.libs/libfftw3.a" "$(BUILD)/fftw-double/threads/.libs/libfftw3_threads.a" "$(EXTERN_LIBS)/"
	cd "$(BUILD)/fftw-float" && "$(ROOT)/extern/fftw-3.3.10/configure" --disable-shared --enable-static --enable-threads --enable-float
	$(MAKE) -C "$(BUILD)/fftw-float/kernel" all
	$(MAKE) -C "$(BUILD)/fftw-float/simd-support" all
	$(MAKE) -C "$(BUILD)/fftw-float/dft" all
	$(MAKE) -C "$(BUILD)/fftw-float/rdft" all
	$(MAKE) -C "$(BUILD)/fftw-float/reodft" all
	$(MAKE) -C "$(BUILD)/fftw-float/api" all
	$(MAKE) -C "$(BUILD)/fftw-float" libfftw3f.la
	$(MAKE) -C "$(BUILD)/fftw-float/threads" libfftw3f_threads.la
	cp "$(BUILD)/fftw-float/.libs/libfftw3f.a" "$(BUILD)/fftw-float/threads/.libs/libfftw3f_threads.a" "$(EXTERN_LIBS)/"

# Keep libLBFGS's configure-time optimisations.  SSE2 is enabled on x86_64;
# Apple Silicon uses native ARM optimisation and compiler vectorisation.
build-lbfgs: | $(EXTERN_LIBS)
	@echo "Building performance-optimised libLBFGS"
	@if [ -f "$(ROOT)/extern/liblbfgs-1.10/config.status" ]; then \
		$(MAKE) -C "$(ROOT)/extern/liblbfgs-1.10" distclean >/dev/null 2>&1 || true; \
	fi
	rm -rf "$(BUILD)/lbfgs"
	mkdir -p "$(BUILD)/lbfgs"
	@arch=`uname -m`; \
	case "$$arch" in \
		x86_64|amd64) cflags='-O3 -DNDEBUG -march=native'; cfg='--enable-sse2' ;; \
		arm64|aarch64) cflags='-O3 -DNDEBUG -mcpu=native'; cfg='' ;; \
		*) cflags='-O3 -DNDEBUG'; cfg='' ;; \
	esac; \
	echo "libLBFGS architecture: $$arch"; \
	echo "libLBFGS CFLAGS: $$cflags"; \
	cd "$(BUILD)/lbfgs" && CFLAGS="$$cflags" "$(ROOT)/extern/liblbfgs-1.10/configure" --disable-shared --enable-static $$cfg
	$(MAKE) -C "$(BUILD)/lbfgs"
	cp "$(BUILD)/lbfgs/lib/.libs/liblbfgs.a" "$(LBFGS_LIB)"

$(EXTERN_LIBS):
	mkdir -p "$@"

$(BIN):
	mkdir -p "$@"

fuda: venv | $(BIN)
	$(MAKE) -C native/fuda_py3 PYTHON_EXE=$(PYTHON) BIN=$(BIN) install
	$(PYTHON) -c 'import sysconfig,shutil,pathlib; dst=pathlib.Path(sysconfig.get_paths()["purelib"]); src=pathlib.Path("native/fuda_py3/py"); [shutil.copy2(p,dst/p.name) for p in src.iterdir() if p.name in {"fuda.py","fudaIO.py","fudalib.so","dataIO.so"}]'
	$(PYTHON) -c 'import fuda, fudaIO, fudalib, dataIO'

spinunidec: extern | $(BIN)
	$(MAKE) -C native/spinUnidec BIN=$(BIN) EXTERN_LIBS=$(EXTERN_LIBS) install

spintools: | $(BIN)
	$(MAKE) -C native/spinTools BIN=$(BIN) install

catia: extern | $(BIN)
	$(MAKE) -C native/CATIA BIN=$(BIN) LIB=$(EXTERN_LIBS)/ install

test: venv
	PYTHONPATH=$(ROOT)/applications $(PYTHON) -m pytest applications/spinDecon/tests

clean:
	$(MAKE) -C native/fuda_py3 clean
	$(MAKE) -C native/spinUnidec clean
	$(MAKE) -C native/spinTools clean
	$(MAKE) -C native/CATIA clean
	rm -rf "$(BUILD)" "$(BIN)lib" "$(ROOT)/dist"

rebuild-extern:
	rm -f $(FFTW_LIBS) $(LBFGS_LIB)
	$(MAKE) extern

distclean: clean
	rm -rf "$(VENV)"
	find "$(BIN)" -mindepth 1 ! -name .gitkeep ! -name spinDecon ! -name spinHub -delete
	find "$(EXTERN_LIBS)" -mindepth 1 ! -name .gitkeep -delete
