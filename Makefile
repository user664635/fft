CC = clang
CFLAGS = -Wall -std=gnu2y -O3 -march=native -fuse-ld=lld -lm -lfftw3f -lSDL3
all: main
	./main
