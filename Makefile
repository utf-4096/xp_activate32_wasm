# compiler options
CC     = clang
CFLAGS = \
	-target wasm32 \
	-nostdlib \
	-mbulk-memory \
	-funroll-loops \
	-Ofast

# xp_activate32 target
XPA_SOURCE_DIR = src/
XPA_SOURCES    = $(wildcard $(XPA_SOURCE_DIR)/*.c)
XPA_OBJECTS    = $(patsubst $(XPA_SOURCE_DIR)/%.c,$(XPA_SOURCE_DIR)/%.o,$(XPA_SOURCES))
XPA_BINARY     = web/assets/wasm/xpa.wasm

.PHONY: clean

$(XPA_BINARY): $(XPA_OBJECTS)
	$(CC) -Wl,--no-entry -Wl,--export-all $(CFLAGS) -o $@ $^

$(XPA_SOURCE_DIR)/%.o: $(XPA_SOURCE_DIR)/%.c
	$(CC) $(CFLAGS) -c -o $@ $<

clean:
	rm $(XPA_OBJECTS) $(XPA_BINARY)
