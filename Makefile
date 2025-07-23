CC = clang
CFLAGS = -std=gnu2y -O3 -march=native -fuse-ld=lld -lm -lSDL3
all: main
	./main
