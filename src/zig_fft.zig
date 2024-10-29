const std = @import("std");
const testing = std.testing;
const List = std.ArrayList;
const Allocator = std.mem.Allocator;
const Complex = std.math.Complex;
const assert = std.debug.assert;
const expect = std.testing.expect;
const exp = std.math.complex.exp;

pub fn Fourier(comptime dataType: type) type {
    return struct {
        const Self = @This();
        allocator: Allocator,

        pub fn init(allocator: Allocator) Self {
            return Self{ .allocator = allocator };
        }

        // Cooley & Tukey algorithm. Length of data has to be a power of 2, i.e.
        // 2, 4, 8, 16, 32, 64, etc.
        pub fn arraylist_dft(self: Self, data: List(Complex(dataType))) !List(Complex(dataType)) {
            const n: u32 = @truncate(data.items.len);
            if (n == 1) {
                var g = List(Complex(dataType)).init(self.allocator);
                try g.append(data.items[0]);
                return g;
            } else {
                assert(n % 2 == 0);
                var g = try List(Complex(dataType)).initCapacity(self.allocator, n >> 1);
                var u = try List(Complex(dataType)).initCapacity(self.allocator, n >> 1);

                // Disperse data in g, u
                for (data.items, 0..) |item, idx| {
                    if (idx & 1 == 0) {
                        try g.append(item);
                    } else {
                        try u.append(item);
                    }
                }

                var g_ = try arraylist_dft(self, g);
                var u_ = try arraylist_dft(self, u);
                g.deinit();
                u.deinit();

                var rpos = List(Complex(dataType)).init(self.allocator);
                var rneg = List(Complex(dataType)).init(self.allocator);

                for (0..g_.items.len) |idx| {
                    const gk = g_.items[idx];
                    const uk = u_.items[idx];
                    const uk_term =
                        uk.mul(exp(Complex(dataType).init(0.0, -2.0 * std.math.pi *
                        @as(dataType, @floatFromInt(idx)) / @as(dataType, @floatFromInt(n)))));
                    try rpos.append(gk.add(uk_term));
                    try rneg.append(gk.sub(uk_term));
                }
                g_.deinit();
                u_.deinit();

                try rpos.appendSlice(rneg.items);
                rneg.deinit();
                return rpos;
            }
        }

        pub fn array_dft(self: Self, n: u32, data: []const Complex(dataType)) ![]Complex(dataType) {
            if (n == 1) {
                const ptr = try self.allocator.alloc(Complex(dataType), 1);
                ptr[0] = data[0];
                return ptr;
            } else {
                assert(n % 2 == 0);
                const g = try self.allocator.alloc(Complex(dataType), n >> 1);
                const u = try self.allocator.alloc(Complex(dataType), n >> 1);

                // Disperse data in g, u
                for (0..n) |idx| {
                    if (idx & 1 == 0) {
                        g[idx >> 1] = data[idx];
                    } else {
                        u[idx >> 1] = data[idx];
                    }
                }

                const g_ = try array_dft(self, n >> 1, g);
                const u_ = try array_dft(self, n >> 1, u);
                self.allocator.free(g);
                self.allocator.free(u);

                const result = try self.allocator.alloc(Complex(dataType), n);

                for (0..(n >> 1)) |idx| {
                    const gk = g_[idx];
                    const uk = u_[idx];
                    const uk_term =
                        uk.mul(exp(Complex(dataType).init(0.0, -2.0 * std.math.pi *
                        @as(dataType, @floatFromInt(idx)) / @as(dataType, @floatFromInt(n)))));
                    result[idx] = gk.add(uk_term);
                    result[idx + (n >> 1)] = gk.sub(uk_term);
                }
                self.allocator.free(g_);
                self.allocator.free(u_);
                return result;
            }
        }
    };
}

fn epsilonTest(val_1: f64, val_2: f64) bool {
    return @abs(val_1 - val_2) < 1e-10;
}

test "dft -0 Hz" {
    var list = List(Complex(f64)).init(std.testing.allocator);
    defer list.deinit();

    const c1 = Complex(f64).init(1, 0);
    try list.append(c1);
    try list.append(c1);
    try list.append(c1);
    try list.append(c1);

    const fft = Fourier(f64).init(std.testing.allocator);
    const result = try fft.arraylist_dft(list);

    try expect(result.items.len == 4);
    try expect(epsilonTest(result.items[0].re, 4));
    try expect(epsilonTest(result.items[0].im, 0));
    try expect(epsilonTest(result.items[1].re, 0));
    try expect(epsilonTest(result.items[1].im, 0));
    try expect(epsilonTest(result.items[2].re, 0));
    try expect(epsilonTest(result.items[2].im, 0));
    try expect(epsilonTest(result.items[3].re, 0));
    try expect(epsilonTest(result.items[3].im, 0));

    result.deinit();
}

test "dft +0 Hz" {
    var list = List(Complex(f64)).init(std.testing.allocator);
    defer list.deinit();

    const c2 = Complex(f64).init(-1, 0);
    try list.append(c2);
    try list.append(c2);
    try list.append(c2);
    try list.append(c2);

    const fft = Fourier(f64).init(std.testing.allocator);
    const result = try fft.arraylist_dft(list);

    try expect(result.items.len == 4);
    try expect(epsilonTest(result.items[0].re, -4));
    try expect(epsilonTest(result.items[0].im, 0));
    try expect(epsilonTest(result.items[1].re, 0));
    try expect(epsilonTest(result.items[1].im, 0));
    try expect(epsilonTest(result.items[2].re, 0));
    try expect(epsilonTest(result.items[2].im, 0));
    try expect(epsilonTest(result.items[3].re, 0));
    try expect(epsilonTest(result.items[3].im, 0));

    result.deinit();
}

test "dft 1/2 Hz" {
    var list = List(Complex(f64)).init(std.testing.allocator);
    defer list.deinit();

    const c0 = Complex(f64).init(0, 0);
    const c1 = Complex(f64).init(1, 0);
    const c2 = Complex(f64).init(-1, 0);
    try list.append(c0);
    try list.append(c1);
    try list.append(c0);
    try list.append(c2);

    const fft = Fourier(f64).init(std.testing.allocator);
    const result = try fft.arraylist_dft(list);

    try expect(result.items.len == 4);
    try expect(epsilonTest(result.items[0].re, 0));
    try expect(epsilonTest(result.items[0].im, 0));
    try expect(epsilonTest(result.items[1].re, 0));
    try expect(epsilonTest(result.items[1].im, -2));
    try expect(epsilonTest(result.items[2].re, 0));
    try expect(epsilonTest(result.items[2].im, 0));
    try expect(epsilonTest(result.items[3].re, 0));
    try expect(epsilonTest(result.items[3].im, 2));
    try expect(epsilonTest(result.items[0].magnitude() +
        result.items[1].magnitude() +
        result.items[2].magnitude() +
        result.items[3].magnitude(), 4));
    result.deinit();
}

test "dft 1 Hz" {
    var list = List(Complex(f64)).init(std.testing.allocator);
    defer list.deinit();

    const c1 = Complex(f64).init(1, 0);
    const c2 = Complex(f64).init(-1, 0);
    try list.append(c1);
    try list.append(c1);
    try list.append(c2);
    try list.append(c2);

    const fft = Fourier(f64).init(std.testing.allocator);
    const result = try fft.arraylist_dft(list);

    try expect(result.items.len == 4);
    try expect(epsilonTest(result.items[0].re, 0));
    try expect(epsilonTest(result.items[0].im, 0));
    try expect(epsilonTest(result.items[1].re, 2));
    try expect(epsilonTest(result.items[1].im, -2));
    try expect(epsilonTest(result.items[2].re, 0));
    try expect(epsilonTest(result.items[2].im, 0));
    try expect(epsilonTest(result.items[3].re, 2));
    try expect(epsilonTest(result.items[3].im, 2));

    result.deinit();
}

test "dft 2 Hz" {
    var list = List(Complex(f64)).init(std.testing.allocator);
    defer list.deinit();

    const c1 = Complex(f64).init(1, 0);
    const c2 = Complex(f64).init(-1, 0);
    try list.append(c1);
    try list.append(c2);
    try list.append(c1);
    try list.append(c2);

    const fft = Fourier(f64).init(std.testing.allocator);
    const result = try fft.arraylist_dft(list);

    try expect(result.items.len == 4);
    try expect(epsilonTest(result.items[0].re, 0));
    try expect(epsilonTest(result.items[0].im, 0));
    try expect(epsilonTest(result.items[1].re, 0));
    try expect(epsilonTest(result.items[1].im, 0));
    try expect(epsilonTest(result.items[2].re, 4));
    try expect(epsilonTest(result.items[2].im, 0));
    try expect(epsilonTest(result.items[3].re, 0));
    try expect(epsilonTest(result.items[3].im, 0));
    try expect(epsilonTest(result.items[0].magnitude() +
        result.items[1].magnitude() +
        result.items[2].magnitude() +
        result.items[3].magnitude(), 4));

    result.deinit();
}

test "array dft 2 Hz" {
    const c1 = Complex(f64).init(1, 0);
    const c2 = Complex(f64).init(-1, 0);
    const list: [4]Complex(f64) = .{ c1, c2, c1, c2 };

    const fft = Fourier(f64).init(std.testing.allocator);
    const result = try fft.array_dft(4, &list);

    try expect(epsilonTest(result[0].re, 0));
    try expect(epsilonTest(result[0].im, 0));
    try expect(epsilonTest(result[1].re, 0));
    try expect(epsilonTest(result[1].im, 0));
    try expect(epsilonTest(result[2].re, 4));
    try expect(epsilonTest(result[2].im, 0));
    try expect(epsilonTest(result[3].re, 0));
    try expect(epsilonTest(result[3].im, 0));
    try expect(epsilonTest(result[0].magnitude() +
        result[1].magnitude() +
        result[2].magnitude() +
        result[3].magnitude(), 4));

    std.testing.allocator.free(result);
}
