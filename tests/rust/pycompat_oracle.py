"""The CPython side of the pycompat differential test.

    python3 -I pycompat_oracle.py N SEED [OP ...]
    python3 -I pycompat_oracle.py ops        the operations, one line

Writes N random cases per operation (the hand-picked edge cases are among
the values drawn) to stdout, each line followed by the line this Python
returns for it. tests/rust/pycompat-diff.sh splits the two, feeds the cases
to `og pycompat-probe` and requires its output to equal the expected lines
byte for byte. Nothing is mounted, so the oracle runs in any image. The case format
is the one rust/og/src/probe.rs reads: tab-separated, floats as the hex of
their bits, text as the hex of its UTF-8.
"""

import csv
import html
import io
import json
import math
import random
import struct
import sys

SPECIAL = [
    0.0, -0.0, math.nan, math.inf, -math.inf, 0.125, 0.375, 2.675, 0.35, 0.5, 1.5, 2.5, -0.5, -2.5,
    1e-05, 1e-04, 1e16, 1e15, 1e17, 1e22, 1e23, 5e-324, 2.2250738585072014e-308, 2.225073858507201e-308,
    1.7976931348623157e308, 9.9999995, 999999.5, 0.1 + 0.2, 1 / 3, 2 / 3, 123456789012345678.0, 0.05, 0.15,
    0.25, 1.005, 1.115, 4.35, 8.345, 99.995, 1234.5, 2.0 ** 52 + 0.5, 2.0 ** 53, 2.0 ** 60, 1e300, -1e-300,
]
OPS = ["repr", "f", "e", "g", "pct", "grpf", "round0", "round", "grp", "zpad", "ljust", "rjust",
       "float", "int", "fsum", "sum", "json", "lines", "readline"]
WHITESPACE = " \t\n\x0b\x0c\r\x1c\x1d\x1e\x1f\x85\xa0\u1680\u2000\u2005\u200a\u2028\u2029\u202f\u205f\u3000"


def bits(x):
    return "%016x" % struct.unpack("<Q", struct.pack("<d", x))[0]


def hexs(s):
    return s.encode("utf-8").hex()


def escape(s):
    return (s.replace("\\", "\\\\").replace("\n", "\\n").replace("\r", "\\r").replace("|", "\\|")
            .replace("\0", "\\0"))


def a_float(r):
    k = r.randrange(9)
    if k == 0:
        return r.choice(SPECIAL) * r.choice([1, -1])
    if k == 1:  # any bit pattern: subnormals, huge, tiny, nan payloads
        return struct.unpack("<d", struct.pack("<Q", r.getrandbits(64)))[0]
    if k == 2:  # subnormal
        return struct.unpack("<d", struct.pack("<Q", r.getrandbits(52) | (r.getrandbits(1) << 63)))[0]
    if k == 3:  # a decimal written with few digits
        return round(r.uniform(-1e6, 1e6), r.randrange(8))
    if k == 4:  # an exact binary tie at some decimal place
        return r.randrange(-100000, 100000) / 2 ** r.randrange(1, 12)
    if k == 5:  # next to a power of ten
        x = 10.0 ** r.randrange(-30, 30)
        return math.nextafter(x, r.choice([0, math.inf])) if r.random() < 0.7 else x
    if k == 6:  # a tie at the last digit of a short decimal
        return float("%d.%d5" % (r.randrange(1000), r.randrange(100)))
    return r.random() * 10.0 ** r.randrange(-30, 30) * r.choice([1, -1])


def a_prec(r):
    return r.choice([r.randrange(18), r.randrange(18), r.randrange(18), 20, 25, 30, 50])


def an_int(r):
    return r.choice([1, -1]) * r.randrange(10 ** r.randrange(1, 38))


def a_char(r):
    k = r.randrange(10)
    if k < 5:
        return chr(r.randrange(32, 127))
    if k == 5:
        return chr(r.randrange(0, 32))
    if k == 6:
        return r.choice('"\\/\x7f|')
    if k == 7:
        return chr(r.randrange(0x80, 0x800))
    if k == 8:
        c = r.randrange(0x800, 0x10000)
        return chr(c) if not 0xD800 <= c < 0xE000 else "\ufffd"
    return chr(r.randrange(0x10000, 0x110000))


def a_str(r, n=12):
    return "".join(a_char(r) for _ in range(r.randrange(n)))


def a_float_text(r):
    k = r.randrange(6)
    if k == 0:
        s = r.choice(["inf", "infinity", "nan", "Inf", "INFINITY", "NaN", "iNf", "infinit", "nan1", "-nan"])
    elif k == 1:
        s = r.choice(["", ".", "1.", ".5", "1e", "e5", "1e+", "1..2", "0x10", "1_", "_1", "1__0", "1_.5", "1._5",
                      "1e1_0", "1_000.000_1", "+-1", "--1", "1 2", "١٢", "1e400", "-1e400", "1E5", "00012", "1d5"])
    elif k == 2:  # underscores between digits, valid or not
        s = repr(a_float(r))
        out = []
        for i, c in enumerate(s):
            out.append(c)
            if r.random() < 0.15 and i + 1 < len(s):
                out.append("_")
        s = "".join(out)
    else:
        x = a_float(r)
        s = r.choice([repr(x), "%.17g" % x, "%.3e" % x, "%f" % x, str(r.randrange(10 ** 6))])
        s = r.choice(["", "+"]) + s if not s.startswith("-") else s
        if r.random() < 0.3:
            s = s.upper()
    if r.random() < 0.3:
        s = "".join(r.choice(WHITESPACE) for _ in range(r.randrange(3))) + s
        s += "".join(r.choice(WHITESPACE) for _ in range(r.randrange(3)))
    return "".join(c for c in s if not ("\u0660" <= c <= "\u0669")) or "x"


def an_int_text(r):
    k = r.randrange(5)
    if k == 0:
        s = r.choice(["", "-", "+", "1.0", "1e3", "1_", "_1", "1__2", "0x1f", "--1", "1 2", "0_0", "000", "+_1"])
    else:
        s = str(an_int(r))
        if r.random() < 0.3:
            s = s[0] + "0" * r.randrange(3) + s[1:] if s[0] in "+-" else "0" * r.randrange(3) + s
        if r.random() < 0.3:
            out = []
            for i, c in enumerate(s):
                out.append(c)
                if r.random() < 0.2 and i + 1 < len(s):
                    out.append("_")
            s = "".join(out)
    if r.random() < 0.3:
        s = r.choice(WHITESPACE) + s + r.choice(WHITESPACE)
    return s


def a_list(r):
    n = r.randrange(1, 20)
    xs = [a_float(r) for _ in range(n)]
    k = r.randrange(4)
    if k == 0:  # cancellation
        xs += [-x for x in xs[: n // 2]] + [r.random() * 1e-20]
    elif k == 1:
        xs = [r.choice([0.1, 0.2, 0.3, 1e100, -1e100, 1.0, 1e-16, 1e16]) for _ in range(n)]
    elif k == 2:
        xs = [x for x in xs if math.isfinite(x)] or [1.0]
    r.shuffle(xs)
    return xs


def a_value(r, depth=0):
    k = r.randrange(10 if depth < 4 else 7)
    if k == 0:
        return None
    if k == 1:
        return r.random() < 0.5
    if k == 2:
        return r.choice([1, -1]) * r.randrange(2 ** r.randrange(1, 100))
    if k == 3:
        return a_float(r)
    if k in (4, 5, 6):
        return a_str(r)
    if k in (7, 8):
        return [a_value(r, depth + 1) for _ in range(r.randrange(5))]
    return {a_str(r, 6): a_value(r, depth + 1) for _ in range(r.randrange(5))}


BAD_JSON = ['', '[1,]', '{"a" 1}', '01', '1.', '[1] x', "{'a': 1}", '"\\x"', '[', '{"a":1,}', 'nul', 'Infinity1',
            '-', '"\x01"', '[1 2]', '{"a":}', '-NaN', '+1', '.5', '1e', '"\\u12"', 'tru', '[]]']


def a_json_text(r):
    k = r.randrange(8)
    if k == 0:
        return r.choice(BAD_JSON)
    if k == 1:  # duplicate keys, NaN and Infinity written by hand
        keys = [r.choice("abc") for _ in range(r.randrange(1, 6))]
        vals = [r.choice(["1", "2.5", "NaN", "-Infinity", "Infinity", "[]", "{}", '"x"', "null", "-0", "1E400"])
                for _ in keys]
        return "{" + ", ".join('"%s": %s' % kv for kv in zip(keys, vals)) + "}"
    v = a_value(r)
    text = json.dumps(v, ensure_ascii=r.random() < 0.5, indent=r.choice([None, None, 0, 1, 2, 4]),
                      separators=r.choice([None, (",", ":"), (", ", ": "), (" , ", " : ")]))
    if r.random() < 0.2:
        text = r.choice([" ", "\n", "\t", "\r\n"]) + text + r.choice([" ", "\n", ""])
    return text


def json_opts(r):
    indent = r.choice(["-", "-", "1", "2", "0", "4"])
    return ",".join([indent, r.choice("dcpi"), r.choice("10"), r.choice("10"), r.choice("1110")])


SEPS = {"d": None, "c": (",", ":"), "p": (", ", ": "), "i": (",", ": ")}


def a_lines_input(r):
    parts = [b"a", b"bc", b"\r", b"\n", b"\r\n", b"\xff", b"\xe2\x82", b"\xc3\xa9", b"\x00", b"|", b"\\", b"\xf0\x9f\x98"]
    return b"".join(r.choice(parts) for _ in range(r.randrange(12)))


def case(op, r):
    """One case line and its Python result."""
    if op in TEXT_CASES:
        return TEXT_CASES[op](r)
    if op == "repr":
        x = a_float(r)
        return [bits(x)], lambda: repr(x)
    if op in ("f", "e", "g", "pct", "grpf"):
        x, p = a_float(r), a_prec(r)
        spec = {"f": "f", "e": "e", "g": "g", "pct": "%", "grpf": ",f"}[op]
        return [bits(x), str(p)], lambda: format(x, "%s.%d%s" % ("," if op == "grpf" else "", p, spec[-1]))
    if op == "round0":
        x = a_float(r)
        return [bits(x)], lambda: str(round(x))
    if op == "round":
        x, n = a_float(r), r.choice([r.randrange(-20, 30), r.randrange(-20, 30), -310, -308, -309, 323, 324, 330])
        return [bits(x), str(n)], lambda: repr(round(x, n))
    if op == "grp":
        n = an_int(r)
        return [str(n)], lambda: format(n, ",")
    if op == "zpad":
        n, w = an_int(r) // 10 ** r.randrange(30), r.randrange(1, 25)
        return [str(n), str(w)], lambda: format(n, "0%dd" % w)
    if op in ("ljust", "rjust"):
        s, w = a_str(r), r.randrange(30)
        return [hexs(s), str(w)], lambda: ("%%-%ds" if op == "ljust" else "%%%ds") % w % s
    if op == "float":
        s = a_float_text(r)
        return [hexs(s)], lambda: repr(float(s))
    if op == "int":
        s = an_int_text(r)
        return [hexs(s)], lambda: str(int(s))
    if op in ("fsum", "sum"):
        xs = a_list(r)
        f = math.fsum if op == "fsum" else sum
        return [",".join(bits(x) for x in xs)], lambda: repr(f(xs))
    if op == "json":
        text, opts = a_json_text(r), json_opts(r)
        indent, sep, ascii_, sort, nan = opts.split(",")

        def run():
            v = json.loads(text)
            return json.dumps(v, indent=None if indent == "-" else int(indent), separators=SEPS[sep],
                              ensure_ascii=ascii_ == "1", sort_keys=sort == "1", allow_nan=nan == "1")
        return [hexs(text), opts], run
    if op == "lines":
        b = a_lines_input(r)

        def run():
            t = io.TextIOWrapper(io.BytesIO(b), encoding="utf-8", errors="replace", newline=None)
            return "|".join(escape(line[:-1] if line.endswith("\n") else line) for line in t)
        return [b.hex()], run
    if op == "readline":
        b = a_lines_input(r)

        def run():
            t = io.TextIOWrapper(io.BytesIO(b), encoding="utf-8", errors="replace", newline=None)
            return "|".join(escape(line) for line in iter(t.readline, ""))
        return [b.hex()], run
    raise ValueError(op)

# csv and text (pycompat::csv, pycompat::text). A list is "-" when empty,
# else its items' hex joined by ","; inside a result, pe() writes the
# characters that give it structure as %XX, so no two results look alike.
CSV_PARTS = ["a", "bc", " ", "\t", ",", ";", "|", '"', '""', "\n", "\r", "\r\n", "", "\xe9", "\U0001f600", "\0",
             "\\", "\x1c", "x" * 5]
DELIMS = ["\t", "\t", ",", ";", "|", " ", "\xe9"]
TEXT_PARTS = [" ", "  ", "\t", "\n", "\r", "\x0b", "\x0c", "\x1c", "\x1f", "\x85", "\xa0", "\u2028", "\u3000",
              "\u200b", "\u180e", "a", "Bc", "\xdf", "\u0130", "\u0131", "\u03a3", "\u0391\u03a3", "\u03c3",
              "\u03c2", "\u01c4", "\u01c5", "\ufb01", "\xb2", "\u0663", "\u2163", "\u0301", "'", "&", "<", ">",
              '"', ",", "--", "\U0001f600", "\U0001d7d8", "\0", "x"]
FIELD_LIMIT = 131072
TEXT_COUNT = {}


def pe(s):
    return "".join("%%%02X" % ord(c) if c in "%[],{}=~\\|" or ord(c) < 32 else c for c in s)


def hexlist(xs):
    return ",".join(x.encode("utf-8").hex() for x in xs) if xs else "-"


def a_text(r, parts=TEXT_PARTS):
    return "".join(r.choice(parts) if r.random() < 0.8 else a_char(r) for _ in range(r.randrange(10)))


def nth(op):
    TEXT_COUNT[op] = TEXT_COUNT.get(op, 0) + 1
    return TEXT_COUNT[op] - 1


def a_csv_text(r):
    return "".join(r.choice(CSV_PARTS) for _ in range(r.randrange(14)))


def csv_lines(r, i):
    """Lines as a file iterator yields them; some cut anywhere, and the
    field size limit's edges first."""
    special = [["x" * FIELD_LIMIT + "\n"], ["x" * (FIELD_LIMIT + 1) + "\n"], ['"' + "x" * FIELD_LIMIT + '"\n'],
               ['"' + "x" * (FIELD_LIMIT + 1)], ["a\t" + "é" * FIELD_LIMIT + "\tb\n"]]
    if i < len(special):
        return special[i]
    text = "".join(a_csv_text(r) + r.choice(["\n", "\n", "\r\n", "\r", ""]) for _ in range(r.randrange(5)))
    if r.random() < 0.3:  # cut anywhere: a \r inside a line, a line without its ending
        cuts = sorted(r.randrange(len(text) + 1) for _ in range(r.randrange(4)))
        return [text[a:b] for a, b in zip([0] + cuts, cuts + [len(text)])]
    return list(io.StringIO(text, newline=""))


def render_list(xs):
    return "[%d:%s]" % (len(xs), ",".join(pe(x) for x in xs))


def render_rows(rows):
    return "".join(render_list(row) for row in rows)


def render_dicts(reader, acc):
    """The header, then each row as it is read, so an error keeps the rows before it."""
    acc.append("~" if reader.fieldnames is None else render_list(reader.fieldnames))
    for d in reader:
        acc.append("{" + ",".join(
            ("~=" + render_list(val)) if key is None
            else pe(key) + "=" + ("~" if val is None else pe(val)) for key, val in d.items()) + "}")


def caught(f):
    """f's result, or what it got to and the class of what it raised."""
    acc = []
    try:
        f(acc)
    except Exception as e:
        acc.append("!" + type(e).__name__)
    return "".join(acc)


def case_csvw(r):
    d = r.choice(DELIMS)
    fields = [a_csv_text(r) for _ in range(r.choice([0, 1, 1, 2, 3, 5]))]

    def run():
        buf = io.StringIO()
        csv.writer(buf, delimiter=d, lineterminator="\n").writerow(fields)
        return buf.getvalue()
    return [hexs(d), hexlist(fields)], run


def case_csvr(r, op="csvr"):
    i, d, q = nth(op), r.choice(DELIMS), r.choice("mn")
    quoting = csv.QUOTE_MINIMAL if q == "m" else csv.QUOTE_NONE
    if op == "dictf":  # bytes read in Python's text mode, as a tool opens a file
        data = "".join(csv_lines(r, i)).encode("utf-8") + r.choice([b"", b"", b"\xff", b"\xe2\x82"])
        arg = data.hex()

        def lines():
            return io.TextIOWrapper(io.BytesIO(data), encoding="utf-8", errors="replace", newline=None)
    else:
        given = csv_lines(r, i)
        arg = hexlist(given)

        def lines():
            return iter(given)
    if op == "csvr":
        def run(acc):
            for row in csv.reader(lines(), delimiter=d, quoting=quoting):
                acc.append(render_rows([row]))
    else:
        def run(acc):
            render_dicts(csv.DictReader(lines(), delimiter=d, quoting=quoting), acc)
    return [hexs(d), q, arg], lambda: caught(run)


def case_split(r):
    s = a_text(r)
    sep = r.choice([None, None, ",", " ", "--", "Σ", "\U0001f600"])
    n = r.choice([-1, -1, 0, 1, 2, 5])
    return [hexs(s), "-" if sep is None else hexs(sep), str(n)], lambda: render_rows([s.split(sep, n)])


def case_strip(op):
    def make(r):
        s = a_text(r)
        chars = r.choice([None, None, "", "x", "a ", "\U0001f600Σ"])
        return [hexs(s), "-" if chars is None else hexs(chars)], lambda: getattr(s, op)(chars)
    return make


def case_str(op):
    def make(r):
        s = a_text(r)
        if op == "is":
            return [hexs(s)], lambda: "".join("1" if f() else "0" for f in (s.isspace, s.isdigit, s.isdecimal, s.isalpha))
        if op == "html":
            return [hexs(s)], lambda: html.escape(s)
        return [hexs(s)], getattr(s, op)
    return make


def case_chars(r):
    """Every code point in turn, a chunk per case, so any N covers them all
    (surrogates excluded: no UTF-8 text holds one)."""
    i = nth("chars")
    k = -(-0x110000 // N_CASES[0])
    start = i * k % 0x110000
    cps = [c for c in range(start, min(start + k, 0x110000)) if not 0xD800 <= c < 0xE000]

    def run():
        out = []
        for c in map(chr, cps):
            f = (c.isspace(), c.isdigit(), c.isdecimal(), c.isalpha())
            out.append("%X[%s,%s]" % (sum(b << n for n, b in enumerate(f)), pe(c.upper()), pe(c.lower())))
        return "".join(out)
    return [str(start), str(k)], run


def case_slice(r):
    s = a_text(r)
    a, b = (r.choice([None, r.randrange(-12, 12)]) for _ in range(2))
    return [hexs(s), "-" if a is None else str(a), "-" if b is None else str(b)], lambda: s[a:b]


N_CASES = [1]
TEXT_CASES = {
    "csvw": case_csvw, "csvr": case_csvr, "dictr": lambda r: case_csvr(r, "dictr"),
    "dictf": lambda r: case_csvr(r, "dictf"), "split": case_split, "strip": case_strip("strip"),
    "lstrip": case_strip("lstrip"), "rstrip": case_strip("rstrip"), "upper": case_str("upper"),
    "lower": case_str("lower"), "is": case_str("is"), "html": case_str("html"), "chars": case_chars,
    "slice": case_slice,
}
OPS += list(TEXT_CASES)


def main():
    if sys.argv[1:] == ["ops"]:
        print(" ".join(OPS))
        return
    n, seed = int(sys.argv[1]), int(sys.argv[2])
    N_CASES[0] = n
    ops = sys.argv[3:] or OPS
    r = random.Random(seed)
    out = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8", newline="\n")
    for op in ops:
        for _ in range(n):
            args, run = case(op, r)
            try:
                res = run()
            except Exception as e:  # the exception class is the result
                res = type(e).__name__
            out.write("\t".join([op] + args) + "\n")
            out.write((res if op in ("lines", "readline") else escape(res)) + "\n")
    out.flush()


if __name__ == "__main__":
    main()
