HTSLIB_static_LDFLAGS = -Wl,-O2 -Wl,--sort-common -Wl,--as-needed -Wl,-z,relro -Wl,-z,now -Wl,--disable-new-dtags -Wl,--gc-sections -Wl,--allow-shlib-undefined -Wl,-rpath,/home/gww/micromamba/envs/GWW/lib -Wl,-rpath-link,/home/gww/micromamba/envs/GWW/lib -L/home/gww/micromamba/envs/GWW/lib
HTSLIB_static_LIBS = -lpthread -lz -lm -lbz2 -llzma
