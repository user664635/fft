CC = clang
CFLAGS = -Wall -std=gnu2y -O3 -ffast-math -march=native -fuse-ld=lld -lm -lfftw3f -lSDL3
all: main
	./main
