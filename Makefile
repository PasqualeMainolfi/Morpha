CC = clang++
FLAGS = -std=c++2b -Wextra -Wall
LIBS = -lfftw3 -lsndfile
SRC_INLCUDE = -I/opt/homebrew/opt/fftw/include -I/opt/homebrew/opt/libsndfile/include
SRC_LIBS = -L/opt/homebrew/opt/fftw/lib -L/opt/homebrew/opt/libsndfile/lib

SOURCE = main.cpp
TARGET = build

.PHONY = all
all = compile run

compile:
	$(CC) $(SOURCE) -o $(TARGET) $(FLAGS) $(LIBS) $(SRC_LIBS) $(SRC_INLCUDE)

run:
	./$(TARGET) && rm -f $(TARGET)
