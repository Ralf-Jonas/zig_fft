# zig_fft
Fast Fourier Transformation in pure Zig

## Example

```zig
// Use ArrayList flavour
const fft = Fourier(f64).init(allocator);
const result = try fft.arraylist_dft(list);
defer result.deinit();

// Use Slice flavour
const c1 = Complex(f64).init(1, 0);
const c2 = Complex(f64).init(-1, 0);
const list: [4]Complex(f64) = .{ c1, c2, c1, c2 };

const fft = Fourier(f64).init(allocator);
const result = try fft.array_dft(4, &list);
defer allocator.free(result);
```

## Description

The FFT works according to the Cooley & Tukey method with arrays of complex numbers, whereby the arrays must have a length of 2^x, i.e. 2, 4, 8, 16, 32, etc. The analogue output signal is then entered into the real part of the complex numbers. To use the FFT 
the analogue output signal is entered in the real part of the complex numbers, the imaginary parts remain zero. The return of the FFT consists of an array of the same length, the first element of which contains the mean value of the input signal, which must be divided by the number of field elements before use. For example, if the input signal is constant 1.0 and contains 4 elements, then the first element of the FFT is a 4.0. This leaves an odd number of fields, with the centre field containing the highest frequency (i.e. the third element for 4 fields). The lowest mapped frequency (and therefore also the frequency spacing between two neighbouring fields) is derived from the reciprocal of the length of the input signal. For example, if the input signal is 20ms long, the first frequency is 50Hz, then 100Hz, 150Hz, etc. The highest frequency (the field in the centre of the array, see above) is determined by the duration from one input field to the next and is exactly half the frequency of the reciprocal value. If the time interval between two input fields is 1 ms, then the highest calculated frequency is not 1 kHz, but 500 Hz (see also Shannon's sampling theorem).

The FFT can work with floating point numbers of different precision, the data type can be f16, f32, f64 or f128. 

```zig
// Use ArrayList flavour
const fft = Fourier(f64).init(allocator);
                    ^^^
```

The FFT is implemented twice: With ArrayList(Complex(dataType)) and with [] Complex(dataType).

## Installation

Add zig_fft to your `build.zig.zon` with the following command:

`zig fetch --save https://github.com/Ralf-Jonas/zig_fft/archive/[commit_hash].tar.gz`

Replace [commit_hash] with the latest commit or tagged release.

Then add the following to your `build.zig`:

```zig
const fft = b.dependency("fft", .{
    .target = target,
    .optimize = optimize,
});
exe.root_module.addImport("fft", zgl.module("fft"));
```

Then import it with `const gl = @import("fft");`, and build as normal with `zig build`.

Examples can be found in the test cases.
