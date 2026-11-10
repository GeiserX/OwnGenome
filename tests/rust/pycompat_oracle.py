"""The CPython side of the pycompat differential test.

    python3 -I pycompat_oracle.py N SEED [OP ...]
    python3 -I pycompat_oracle.py ops        the operations, one line

Writes N random cases per operation (the hand-picked edge cases are among
the values drawn) to stdout, each line followed by the line this Python
returns for it. tests/rust/pycompat-diff.sh splits the two, feeds the cases
to `og pycompat-probe` and requires its output to equal the expected lines
byte for byte. The tree operations (glob, mtimes) work on a folder the
oracle fills under $PYCOMPAT_TREE and og reads after it; nothing else is
mounted, so the oracle runs in any image. The case format
is the one rust/og/src/probe.rs reads: tab-separated, floats as the hex of
their bits, text as the hex of its UTF-8.
"""

import glob
import io
import json
import math
import os
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
# pycompat::value and pycompat::mtime. The tree operations work on a folder
# the oracle fills under $PYCOMPAT_TREE, which og then reads in place.
VALUE_OPS = ["truthy", "vstr", "vrepr", "getor", "fnmatch", "escape", "sortnames"]
TREE_OPS = ["glob", "first", "newest", "mtime", "mtimens"]
OPS += VALUE_OPS + TREE_OPS
# A tree case walks folders, so it gets one case in TREE_SHARE.
TREE_SHARE = 20
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


FALSY = [None, False, 0, 0.0, -0.0, "", [], {}, math.nan, "0", [0], {"": None}, 1e-300]
# Characters a name or a pattern is drawn from: the ones fnmatch treats
# specially, a non-ASCII and an astral one, and a stray byte as its surrogate.
NAME_CHARS = "abc-]![^\\*?.&~|" + "\u00e9\U0001d400\udcf5"
NAMES = ["a", "b", "ab", "a.txt", "b.txt", "a.vcf.gz", ".h", ".h.txt", "_x", "[a]", "a*", "a?", "-", "!x",
         "\u00e9", "\u4e2d", "\U0001d400", "x y", "A", "z\udcff", "s_eh.vcf", "**", "[!a]"]
RAW_BYTES = [b"a", b"b", b"A", b"\x80", b"\xc3\xa9", b"\xf0\x9d\x90\x80", b"\xf5", b"\xed\xa0\x80", b"\xe4\xb8",
             b"\xff", b"\xee\x80\x80", b"/"]
NS_BASE = 1_700_000_000 * 10 ** 9


def fs_hex(s):
    return os.fsencode(s).hex()


def a_latin1_str(r, n=8):
    picks = ["a", "'", '"', "\\", "\t", "\n", "\r", "\x00", "\x1f", " ", "\x7f"]
    return "".join(r.choice(picks + [chr(r.randrange(0x80, 0x100))]) for _ in range(r.randrange(n)))


def a_py_value(r, depth=0):
    """A JSON value whose strings inside lists and dicts are mostly Latin-1,
    so most cases get past the known limit of value::py_repr."""
    k = r.randrange(10 if depth < 3 else 6)
    if k == 0:
        return r.choice(FALSY)
    if k == 1:
        return a_value(r, 4)
    if k in (2, 3):
        return a_str(r) if depth == 0 else a_latin1_str(r)
    if k in (4, 5):
        return a_float(r)
    if k in (6, 7):
        return [a_py_value(r, depth + 1) for _ in range(r.randrange(4))]
    return {a_latin1_str(r, 5): a_py_value(r, depth + 1) for _ in range(r.randrange(4))}


def needs_unicode_repr(v, top):
    """Whether str() (TOP true) or repr() of V reprs a string with a
    character above U+00FF: the known limit, where value::py_repr refuses."""
    if isinstance(v, str):
        return not top and any(ord(c) > 0xFF for c in v)
    if isinstance(v, list):
        return any(needs_unicode_repr(x, False) for x in v)
    if isinstance(v, dict):
        return any(needs_unicode_repr(k, False) or needs_unicode_repr(x, False) for k, x in v.items())
    return False


def a_name(r, chars, n=6):
    return "".join(r.choice(chars) for _ in range(r.randrange(n)))


def value_case(op, r):
    if op == "truthy":
        text = json.dumps(r.choice([r.choice(FALSY), a_value(r)]))
        return [hexs(text)], lambda: str(bool(json.loads(text)))
    if op in ("vstr", "vrepr"):
        text = json.dumps(a_py_value(r))
        f = str if op == "vstr" else repr

        def run():
            v = json.loads(text)
            return "Unsupported" if needs_unicode_repr(v, op == "vstr") else f(v)
        return [hexs(text)], run
    if op == "getor":
        d = r.choice([{a_str(r, 4): a_value(r, 2) for _ in range(r.randrange(5))}] * 4 + [a_value(r, 3)])
        keys = list(d) if isinstance(d, dict) else []
        key = r.choice(keys) if keys and r.random() < 0.8 else a_str(r, 4)
        dtext, default = json.dumps(d), json.dumps(r.choice([{}, [], None, "", 0, a_value(r, 3)]))
        return ([hexs(dtext), hexs(key), hexs(default)],
                lambda: json.dumps(json.loads(dtext).get(key) or json.loads(default)))
    if op == "fnmatch":
        name, pat = a_name(r, NAME_CHARS.replace("[", "").replace("*", "") + "ab"), a_name(r, NAME_CHARS, 9)
        return [fs_hex(name), fs_hex(pat)], lambda: str(fnmatch_case(name, pat))
    if op == "escape":
        s = a_name(r, "a*?[]!/.\u00e9", 8) + a_str(r, 4)
        return [hexs(s)], lambda: glob.escape(s)
    if op == "sortnames":
        names = [b"".join(r.choice(RAW_BYTES) for _ in range(r.randrange(1, 4))) for _ in range(r.randrange(6))]
        return ([",".join(b.hex() for b in names)],
                lambda: ",".join(os.fsencode(s).hex() for s in sorted(os.fsdecode(b) for b in names)))
    raise ValueError(op)


def fnmatch_case(name, pat):
    import fnmatch
    return fnmatch.fnmatchcase(name, pat)


def build_tree(root, r):
    """Fills ROOT with folders, files and symlinks (to a file, to a folder,
    dangling), each with an mtime drawn so that some differ by less than a
    float can tell apart. Returns every path made."""
    made = []

    def fill(d, depth):
        for name in r.sample(NAMES, r.randrange(2, 9)):
            p = os.path.join(d, name)
            try:
                if depth < 3 and r.random() < 0.35:
                    os.mkdir(p)
                    made.append(p)
                    fill(p, depth + 1)
                else:
                    open(p, "w").close()
                    made.append(p)
            except OSError:  # a name this file system refuses (a stray byte on macOS)
                pass
    fill(root, 0)
    dirs = [p for p in made if os.path.isdir(p)]
    files = [p for p in made if os.path.isfile(p)]
    links = [("ln_file", r.choice(files) if files else "/nonexistent"), ("ln_gone", os.path.join(root, "gone")),
             ("ln_dir", r.choice(dirs) if dirs else root + "/missing")]
    for name, target in links:
        os.symlink(target, os.path.join(root, name))
        made.append(os.path.join(root, name))
    # Deepest first, so setting a folder's time comes after its last change.
    for p in sorted(made + [root], key=lambda p: -p.count("/")):
        ns = NS_BASE + r.choice([0, 0, 1, 2, 100, 239, 240, 241, 499, 10 ** 9, -10 ** 9, r.randrange(10 ** 12)])
        if r.random() < 0.05:
            ns = -r.randrange(10 ** 18)
        os.utime(p, ns=(ns, ns), follow_symlinks=False)
    return made


def a_rel_pattern(r):
    parts = []
    for _ in range(r.randrange(1, 4)):
        k = r.randrange(9)
        if k == 0:
            parts.append("**")
        elif k == 1:
            parts.append(r.choice(NAMES))
        elif k == 2:
            parts.append(glob.escape(r.choice(NAMES)))
        elif k == 3:
            parts.append("")
        elif k == 4:
            parts.append(a_name(r, "ab*?.[]!-\u00e9", 5))
        else:
            parts.append(r.choice(["*", "?", "*.txt", "*.vcf.gz", "[ab]*", "[!a]*", ".*", "*_eh.vcf", "a*", "[a-c]",
                                   "[]]*", "[!]*", "ln_*", "*\udcff"]))
    # ".." would leave the tree for the folder the test writes its own files in.
    return "/".join(p if p != ".." else "." for p in parts) + ("/" if r.random() < 0.1 else "")


def tree_case(op, r, root, made):
    if op == "glob":
        pat, rec = root + "/" + a_rel_pattern(r), r.random() < 0.7
        return [fs_hex(pat), str(int(rec))], lambda: "|".join(fs_hex(p) for p in glob.glob(pat, recursive=rec))
    if op in ("first", "newest"):
        pats = [root + "/" + a_rel_pattern(r) for _ in range(r.randrange(1, 4))]

        def run():  # collect_summary's first_existing and newest
            if op == "first":
                for p in pats:
                    for f in sorted(glob.glob(p, recursive=True)):
                        if os.path.isfile(f):
                            return fs_hex(f)
                return "None"
            found = [f for p in pats for f in glob.glob(p, recursive=True) if os.path.isfile(f)]
            return fs_hex(max(found, key=os.path.getmtime)) if found else "None"
        return [",".join(fs_hex(p) for p in pats)], run
    if op in ("mtime", "mtimens"):
        p = r.choice(made + [root])
        p = r.choice([p, p, p, p + "/x", root + "/missing"])
        if op == "mtime":
            return [fs_hex(p)], lambda: repr(os.path.getmtime(p))
        return [fs_hex(p)], lambda: str(os.stat(p).st_mtime_ns)
    raise ValueError(op)


def case(op, r):
    """One case line and its Python result."""
    if op in VALUE_OPS:
        return value_case(op, r)
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


def main():
    if sys.argv[1:] == ["ops"]:
        print(" ".join(OPS))
        return
    n, seed = int(sys.argv[1]), int(sys.argv[2])
    ops = sys.argv[3:] or OPS
    r = random.Random(seed)
    out = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8", newline="\n")
    root = made = None
    if set(ops) & set(TREE_OPS):
        root = os.environ.get("PYCOMPAT_TREE") or sys.exit("PYCOMPAT_TREE: an empty folder for the tree operations")
        made = build_tree(root, random.Random(seed))
    for op in ops:
        for _ in range(max(1, n // TREE_SHARE) if op in TREE_OPS else n):
            args, run = tree_case(op, r, root, made) if op in TREE_OPS else case(op, r)
            try:
                res = run()
            except Exception as e:  # the exception class is the result
                res = type(e).__name__
            out.write("\t".join([op] + args) + "\n")
            out.write((res if op in ("lines", "readline") else escape(res)) + "\n")
    out.flush()


if __name__ == "__main__":
    main()
