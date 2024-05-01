# xp_activate32_wasm

This is a port of `xp_activate32.exe` to the browser using WebAssembly, built from
[this source](https://archive.org/details/xp_activate32_src).

## Requirements

- make
- a version of clang modern enough to support the `wasm32` target

## Building

```
git clone https://github.com/utf-4096/xp_activate32_wasm.git
make -j$(nproc)
```

This will compile and place the binary at `web/assets/wasm/xpa.wasm`. You can now upload
the `web/` folder to your webhost.

## Acknowledgements

- Famfamfam for the Silk icons at [`web/assets/icons/`](web/assets/icons/) and
  [`web/favicon.ico`](web/favicon.ico).
- [maniekx86/xpactivate32_php](https://github.com/maniekx86/xpactivate32_php)
  for the idea of running `xp_activate32` in the browser.
