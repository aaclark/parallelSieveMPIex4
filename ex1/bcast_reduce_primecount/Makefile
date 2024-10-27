# Makefile for seives program

# Compiler
CC = mpicc

# Compiler flags
CFLAGS = -lm

# Adding additional flags from the command line or the Makefile
CFLAGS += $(OUTPUT_FLAG)

# Target executable
TARGET = sieve_prime_count_reduce

# Source files
SRC = sieve_prime_count_reduce.c

# Default target
all: $(TARGET)

# Build target
$(TARGET): $(SRC)
	$(CC) $(SRC) -o $(TARGET) $(CFLAGS)

# Clean target
clean:
	rm -f $(TARGET)

# Run target with arguments
run: $(TARGET)
	@echo "Running $(TARGET) with arguments:"
	@read -p "Enter processes: " p; \
        read -p "Enter Max: " max; \
	mpiexec -np $$p ./$(TARGET) $$max

.PHONY: all clean run
