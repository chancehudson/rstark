# rstark [![CircleCI](https://img.shields.io/circleci/build/github/vimwitch/rstark/main)](https://app.circleci.com/pipelines/github/vimwitch/rstark)

WIP: ZK-STARK prover written in rust.

## Building wasm

First install [`wasm-pack`](https://rustwasm.github.io/wasm-pack/installer/), then run `wasm-pack build` in the project root. A wasm file with some supporting JS will be generated in `pkg/`.

## Profiling

Install the flamegraph crate globally using `cargo flamegraph install`. Ensure that either `perf` or `dtrace` is present in your `PATH`.

To generate a flamegraph SVG:

```
cargo flamegraph --example squares
```

If you're on Mac add the `--root` flag to the above command.

Open the SVG file in a browser to interact with the stack elements.

## Browser support

See [here](https://rstark.net) and [here](https://github.com/vimwitch/benchstark) for a React based project that uses the wasm built in this repo.
