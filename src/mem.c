char host_buf[128];

/* there are better ways to pass data from js to wasm of course, but this is simple */
char* get_buffer() {
    return host_buf;
}
