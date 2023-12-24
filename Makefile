CFLAGS = -g -Wall -Wpedantic -Wunused-const-variable -Werror -Wno-attributes -fsigned-zeros
LDFLAGS = -lraylib -lm
TARGET = tile

all:
	$(CC) $(CFLAGS) $(LDFLAGS) $(TARGET).c -o $(TARGET)
debug:
	$(CC) -DDEBUG $(CFLAGS) $(LDFLAGS) $(TARGET).c -o $(TARGET)

.PHONY: all
