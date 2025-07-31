CFLAGS=-I/usr/local/include -g -Wall
LDFLAGS=-L/usr/local/lib -lgsl -lgslcblas -lm
CC=gcc

TARGET= feautrier

all: $(TARGET)

$(TARGET): $(TARGET).c
	$(CC) -o $(TARGET) $(TARGET).c $(CFLAGS) $(LDFLAGS)

clean:
	$(RM) $(TARGET)
