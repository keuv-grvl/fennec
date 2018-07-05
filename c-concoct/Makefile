CFLAGS=-I${CONDA_PREFIX}/include/python3.6m -pthread -fPIC -Wall -Wextra -pedantic
LDFLAGS=-L${CONDA_PREFIX}/lib -lpython3.6m -lpthread -lm -lgsl -lgslcblas

SRC=vbgmmmodule.c
OBJ=$(SRC:.c=.o)

SOFILE=vbgmmmodule.so
TARGET=../vbgmm.cpython-36m-x86_64-linux-gnu.so

all: $(TARGET)

$(TARGET): $(SOFILE)
	cp -f $^ $@

$(SOFILE): $(OBJ)
	$(CC) $(CFLAGS) -shared -o $@ $(OBJ) $(LDFLAGS)

%.o: %.c
	$(CC) $(CFLAGS) -o $@ -c $^
