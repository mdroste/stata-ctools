"""A small, dependency-free Parquet writer for ctools' import tests.

No Parquet library is assumed to be installed, so this module writes the
format directly from the specification (parquet-format, Thrift compact
protocol footer). It covers what the tests exercise: every physical type,
logical and converted type annotations, nulls, row groups, several pages per
chunk, data pages v1 and v2, PLAIN / dictionary / RLE / BIT_PACKED / DELTA_* /
BYTE_STREAM_SPLIT encodings and the codecs UNCOMPRESSED, SNAPPY, GZIP, ZSTD
(Python 3.14's compression.zstd), LZ4 (Hadoop framing), LZ4_RAW and BROTLI
(stored meta-blocks only). Compressors for Snappy and LZ4 are greedy hash
matchers written here, so the files contain back-references as well as
literals.

Schema nesting is available through the low-level `write_file`, which takes
schema elements and per-leaf levels as given, for tests of nested schemas.
"""

import struct
import zlib

# Physical types
BOOLEAN, INT32, INT64, INT96, FLOAT, DOUBLE, BYTE_ARRAY, FIXED_LEN_BYTE_ARRAY = range(8)
# Repetition
REQUIRED, OPTIONAL, REPEATED = range(3)
# Encodings
PLAIN, PLAIN_DICTIONARY, RLE, BIT_PACKED = 0, 2, 3, 4
DELTA_BINARY_PACKED, DELTA_LENGTH_BYTE_ARRAY, DELTA_BYTE_ARRAY = 5, 6, 7
RLE_DICTIONARY, BYTE_STREAM_SPLIT = 8, 9
# Codecs
UNCOMPRESSED, SNAPPY, GZIP, LZO, BROTLI, LZ4, ZSTD, LZ4_RAW = range(8)
CODECS = {
    "uncompressed": UNCOMPRESSED,
    "snappy": SNAPPY,
    "gzip": GZIP,
    "brotli": BROTLI,
    "lz4": LZ4,
    "zstd": ZSTD,
    "lz4_raw": LZ4_RAW,
}
# Converted types
CT = dict(
    UTF8=0, MAP=1, MAP_KEY_VALUE=2, LIST=3, ENUM=4, DECIMAL=5, DATE=6, TIME_MILLIS=7,
    TIME_MICROS=8, TIMESTAMP_MILLIS=9, TIMESTAMP_MICROS=10, UINT_8=11, UINT_16=12,
    UINT_32=13, UINT_64=14, INT_8=15, INT_16=16, INT_32=17, INT_64=18, JSON=19, BSON=20,
    INTERVAL=21,
)

# ---------------------------------------------------------------------------
# Thrift compact protocol (writer side)
# ---------------------------------------------------------------------------
T_TRUE, T_FALSE, T_BYTE, T_I16, T_I32, T_I64, T_DOUBLE, T_BINARY = 1, 2, 3, 4, 5, 6, 7, 8
T_LIST, T_SET, T_MAP, T_STRUCT = 9, 10, 11, 12


def uleb(n):
    out = bytearray()
    while True:
        b = n & 0x7F
        n >>= 7
        if n:
            out.append(b | 0x80)
        else:
            out.append(b)
            return bytes(out)


def zigzag(n):
    return (n << 1) ^ (n >> 63)


class Struct:
    """Fields as (id, kind, value); kind in bool,i8,i16,i32,i64,bin,struct,list."""

    def __init__(self, *fields):
        self.fields = [f for f in fields if f is not None and f[2] is not None]


def _ctype(kind):
    return {"bool": T_TRUE, "i8": T_BYTE, "i16": T_I16, "i32": T_I32, "i64": T_I64,
            "bin": T_BINARY, "struct": T_STRUCT, "double": T_DOUBLE}[kind]


def _value(kind, v):
    if kind in ("i16", "i32", "i64"):
        return uleb(zigzag(v))
    if kind == "i8":
        return struct.pack("b", v)
    if kind == "double":
        return struct.pack("<d", v)
    if kind == "bin":
        if isinstance(v, str):
            v = v.encode("utf-8")
        return uleb(len(v)) + v
    if kind == "struct":
        return encode_struct(v)
    if kind == "bool":  # inside a list
        return b"\x01" if v else b"\x02"
    raise ValueError(kind)


def encode_struct(s):
    out = bytearray()
    last = 0
    for fid, kind, v in sorted(s.fields, key=lambda f: f[0]):
        if kind.startswith("list:"):
            ctype = T_LIST
        elif kind == "bool":
            ctype = T_TRUE if v else T_FALSE
        else:
            ctype = _ctype(kind)
        delta = fid - last
        if 0 < delta <= 15:
            out.append((delta << 4) | ctype)
        else:
            out.append(ctype)
            out += uleb(zigzag(fid))
        last = fid
        if kind.startswith("list:"):
            ek = kind[5:]
            et = _ctype(ek)
            if et == T_TRUE:
                et = T_TRUE
            n = len(v)
            out.append((n << 4) | et if n < 15 else 0xF0 | et)
            if n >= 15:
                out += uleb(n)
            for item in v:
                out += _value(ek, item)
        elif kind != "bool":
            out += _value(kind, v)
    out.append(0)
    return bytes(out)


# ---------------------------------------------------------------------------
# Logical types
# ---------------------------------------------------------------------------
def _unit(unit):
    return Struct((1 if unit == "ms" else 2 if unit == "us" else 3, "struct", Struct()))


def logical(kind, **kw):
    """LogicalType union as a Struct."""
    ids = dict(STRING=1, MAP=2, LIST=3, ENUM=4, DECIMAL=5, DATE=6, TIME=7, TIMESTAMP=8,
               INTEGER=10, UNKNOWN=11, JSON=12, BSON=13, UUID=14, FLOAT16=15)
    if kind == "DECIMAL":
        body = Struct((1, "i32", kw["scale"]), (2, "i32", kw["precision"]))
    elif kind in ("TIME", "TIMESTAMP"):
        body = Struct((1, "bool", kw["utc"]), (2, "struct", _unit(kw["unit"])))
    elif kind == "INTEGER":
        body = Struct((1, "i8", kw["bits"]), (2, "bool", kw["signed"]))
    else:
        body = Struct()
    return Struct((ids[kind], "struct", body))


def schema_element(name, ptype=None, repetition=None, children=None, converted=None,
                   logical_type=None, type_length=None, scale=None, precision=None,
                   field_id=None):
    return Struct(
        (1, "i32", ptype),
        (2, "i32", type_length),
        (3, "i32", repetition),
        (4, "bin", name),
        (5, "i32", children),
        (6, "i32", CT[converted] if isinstance(converted, str) else converted),
        (7, "i32", scale),
        (8, "i32", precision),
        (9, "i32", field_id),
        (10, "struct", logical_type),
    )


# ---------------------------------------------------------------------------
# Bit packing and RLE / bit-packed hybrid
# ---------------------------------------------------------------------------
def bitpack(values, width):
    """LSB-first bit packing (RLE hybrid, DELTA miniblocks)."""
    acc = 0
    nbits = 0
    out = bytearray()
    for v in values:
        acc |= (v & ((1 << width) - 1)) << nbits
        nbits += width
        while nbits >= 8:
            out.append(acc & 0xFF)
            acc >>= 8
            nbits -= 8
    if nbits:
        out.append(acc & 0xFF)
    return bytes(out)


def bitpack_msb(values, width):
    """The deprecated BIT_PACKED level encoding: MSB first."""
    acc = 0
    nbits = 0
    out = bytearray()
    for v in values:
        for b in range(width - 1, -1, -1):
            acc = (acc << 1) | ((v >> b) & 1)
            nbits += 1
            if nbits == 8:
                out.append(acc)
                acc = nbits = 0
    if nbits:
        out.append(acc << (8 - nbits))
    return bytes(out)


def rle_hybrid(values, width, mode="mixed"):
    """RLE/bit-packed hybrid without length prefix.

    mode: "mixed" uses RLE runs of 8 or more equal values and bit-packed groups
    otherwise; "rle" uses only RLE runs; "packed" only bit-packed groups.
    """
    nbytes = (width + 7) // 8
    out = bytearray()
    values = list(values)
    n = len(values)

    def rle_run(v, count):
        return uleb(count << 1) + int(v).to_bytes(nbytes, "little") if nbytes else uleb(count << 1)

    def packed(chunk):
        groups = (len(chunk) + 7) // 8
        padded = chunk + [0] * (groups * 8 - len(chunk))
        return uleb((groups << 1) | 1) + bitpack(padded, width)

    if mode == "rle":
        i = 0
        while i < n:
            j = i
            while j < n and values[j] == values[i]:
                j += 1
            out += rle_run(values[i], j - i)
            i = j
        return bytes(out)
    if mode == "packed":
        return packed(values) if n else b""
    i = 0
    pending = []
    while i < n:
        j = i
        while j < n and values[j] == values[i]:
            j += 1
        if j - i >= 8 and len(pending) % 8 == 0:
            if pending:
                out += packed(pending)
                pending = []
            out += rle_run(values[i], j - i)
            i = j
        else:
            pending.append(values[i])
            i += 1
    if pending:
        out += packed(pending)
    return bytes(out)


# ---------------------------------------------------------------------------
# DELTA encodings
# ---------------------------------------------------------------------------
def delta_binary_packed(values, bits=64, block=128, miniblocks=4):
    mask = (1 << bits) - 1
    out = bytearray(uleb(block) + uleb(miniblocks) + uleb(len(values)))
    first = values[0] if values else 0
    out += uleb(zigzag(first) & ((1 << 64) - 1))
    per_mini = block // miniblocks
    deltas = []
    for a, b in zip(values, values[1:]):
        d = (b - a) & mask
        if d >> (bits - 1):
            d -= 1 << bits
        deltas.append(d)
    for start in range(0, len(deltas), block):
        chunk = deltas[start:start + block]
        mn = min(chunk)
        out += uleb(zigzag(mn) & ((1 << 64) - 1))
        rel = [(d - mn) & mask for d in chunk]
        widths = []
        packs = []
        for m in range(miniblocks):
            part = rel[m * per_mini:(m + 1) * per_mini]
            if not part:
                widths.append(0)
                continue
            w = max(part).bit_length()
            widths.append(w)
            part = part + [0] * (per_mini - len(part))
            packs.append(bitpack(part, w))
        out += bytes(widths)
        for p in packs:
            out += p
    return bytes(out)


def delta_length_byte_array(values):
    return delta_binary_packed([len(v) for v in values], 32) + b"".join(values)


def delta_byte_array(values):
    prefixes, suffixes = [], []
    prev = b""
    for v in values:
        k = 0
        while k < len(v) and k < len(prev) and v[k] == prev[k]:
            k += 1
        prefixes.append(k)
        suffixes.append(v[k:])
        prev = v
    return delta_binary_packed(prefixes, 32) + delta_length_byte_array(suffixes)


# ---------------------------------------------------------------------------
# PLAIN values
# ---------------------------------------------------------------------------
def plain(ptype, values, type_length=None):
    if ptype == BOOLEAN:
        return bitpack([1 if v else 0 for v in values], 1)
    if ptype == INT32:
        return b"".join(struct.pack("<i", v) for v in values)
    if ptype == INT64:
        return b"".join(struct.pack("<q", v) for v in values)
    if ptype == INT96:
        return b"".join(v if isinstance(v, bytes) else struct.pack("<qi", *v) for v in values)
    if ptype == FLOAT:
        return b"".join(struct.pack("<f", v) for v in values)
    if ptype == DOUBLE:
        return b"".join(struct.pack("<d", v) for v in values)
    if ptype == BYTE_ARRAY:
        return b"".join(struct.pack("<I", len(v)) + v for v in values)
    if ptype == FIXED_LEN_BYTE_ARRAY:
        for v in values:
            assert len(v) == type_length
        return b"".join(values)
    raise ValueError(ptype)


def width_of(ptype, type_length):
    return {INT32: 4, INT64: 8, FLOAT: 4, DOUBLE: 8, INT96: 12}.get(ptype, type_length)


def encode_values(ptype, values, encoding, type_length=None):
    if encoding == PLAIN:
        return plain(ptype, values, type_length)
    if encoding == RLE:  # booleans
        assert ptype == BOOLEAN
        body = rle_hybrid([1 if v else 0 for v in values], 1)
        return struct.pack("<I", len(body)) + body
    if encoding == DELTA_BINARY_PACKED:
        return delta_binary_packed(values, 32 if ptype == INT32 else 64)
    if encoding == DELTA_LENGTH_BYTE_ARRAY:
        return delta_length_byte_array(values)
    if encoding == DELTA_BYTE_ARRAY:
        return delta_byte_array(values)
    if encoding == BYTE_STREAM_SPLIT:
        raw = plain(ptype, values, type_length)
        w = width_of(ptype, type_length)
        n = len(values)
        return bytes(raw[i * w + k] for k in range(w) for i in range(n))
    raise ValueError(encoding)


# ---------------------------------------------------------------------------
# Compression
# ---------------------------------------------------------------------------
def snappy_compress(data):
    out = bytearray(uleb(len(data)))
    n = len(data)
    table = {}
    lit_start = 0
    i = 0

    def literal(a, b):
        while a < b:
            m = min(b - a, 65536)
            L = m - 1
            if L < 60:
                out.append(L << 2)
            elif L < 256:
                out.extend(bytes([60 << 2, L]))
            else:
                out.extend(bytes([61 << 2]) + struct.pack("<H", L))
            out.extend(data[a:a + m])
            a += m

    while i + 4 <= n:
        key = data[i:i + 4]
        cand = table.get(key)
        table[key] = i
        if cand is not None and i - cand <= 65535:
            length = 4
            while i + length < n and data[cand + length] == data[i + length]:
                length += 1
            literal(lit_start, i)
            off = i - cand
            rest = length
            while rest > 0:
                m = min(rest, 64)
                if rest - m and rest - m < 4:
                    m = rest - 4
                if 4 <= m <= 11 and off < 2048:
                    out.append(((off >> 8) << 5) | ((m - 4) << 2) | 1)
                    out.append(off & 0xFF)
                else:
                    out.append(((m - 1) << 2) | 2)
                    out += struct.pack("<H", off)
                rest -= m
            i += length
            lit_start = i
        else:
            i += 1
    literal(lit_start, n)
    return bytes(out)


def lz4_block_compress(data):
    n = len(data)
    out = bytearray()
    table = {}
    anchor = 0
    i = 0
    limit = n - 12  # last match must start 12 bytes before the end

    def lenext(x):
        b = bytearray()
        while x >= 255:
            b.append(255)
            x -= 255
        b.append(x)
        return b

    while i < limit:
        key = data[i:i + 4]
        cand = table.get(key)
        table[key] = i
        if cand is not None and i - cand <= 65535:
            length = 4
            while i + length < n - 5 and data[cand + length] == data[i + length]:
                length += 1
            lit = i - anchor
            ml = length - 4
            token = (min(lit, 15) << 4) | min(ml, 15)
            out.append(token)
            if lit >= 15:
                out += lenext(lit - 15)
            out += data[anchor:i]
            out += struct.pack("<H", i - cand)
            if ml >= 15:
                out += lenext(ml - 15)
            i += length
            anchor = i
        else:
            i += 1
    lit = n - anchor
    out.append(min(lit, 15) << 4)
    if lit >= 15:
        out += lenext(lit - 15)
    out += data[anchor:]
    return bytes(out)


def brotli_stored(data):
    """Brotli stream of uncompressed meta-blocks (RFC 7932 section 9.2)."""
    bits = []

    def put(v, n):
        for k in range(n):
            bits.append((v >> k) & 1)

    def flush_bytes():
        while len(bits) % 8:
            bits.append(0)
        out = bytearray()
        for k in range(0, len(bits), 8):
            out.append(sum(bits[k + j] << j for j in range(8)))
        bits.clear()
        return bytes(out)

    out = bytearray()
    put(0, 1)  # WBITS = 16
    pos = 0
    while pos < len(data):
        m = min(65536, len(data) - pos)
        put(0, 1)  # ISLAST
        put(0, 2)  # MNIBBLES = 4
        put(m - 1, 16)
        put(1, 1)  # ISUNCOMPRESSED
        out += flush_bytes()
        out += data[pos:pos + m]
        pos += m
    put(1, 1)  # ISLAST
    put(1, 1)  # ISLASTEMPTY
    out += flush_bytes()
    return bytes(out)


def compress(codec, data):
    if codec == UNCOMPRESSED:
        return data
    if codec == SNAPPY:
        return snappy_compress(data)
    if codec == GZIP:
        c = zlib.compressobj(6, zlib.DEFLATED, 31)
        return c.compress(data) + c.flush()
    if codec == ZSTD:
        from compression import zstd

        return zstd.compress(data, level=3)
    if codec == LZ4_RAW:
        return lz4_block_compress(data)
    if codec == LZ4:  # Hadoop framing: BE uncompressed size, BE compressed size
        block = lz4_block_compress(data)
        return struct.pack(">II", len(data), len(block)) + block
    if codec == BROTLI:
        return brotli_stored(data)
    raise ValueError(codec)


# ---------------------------------------------------------------------------
# Pages and column chunks
# ---------------------------------------------------------------------------
def page_header(ptype, usize, csize, body):
    return encode_struct(Struct((1, "i32", ptype), (2, "i32", usize), (3, "i32", csize), body))


def levels_v1(levels, max_level, encoding):
    if max_level == 0:
        return b""
    width = max_level.bit_length()
    if encoding == BIT_PACKED:
        return bitpack_msb(levels, width)
    body = rle_hybrid(levels, width)
    return struct.pack("<I", len(body)) + body


class Leaf:
    """One leaf column's data for one row group.

    levels: definition levels (None when max_def == 0); rep: repetition levels.
    values: the non-null values in order.
    """

    def __init__(self, ptype, values, def_levels=None, rep_levels=None, max_def=0, max_rep=0,
                 path=("x",), type_length=None, encoding=PLAIN, dictionary=False,
                 page_version=1, codec=UNCOMPRESSED, page_rows=None, level_encoding=RLE,
                 dict_encoding=RLE_DICTIONARY, dict_page_encoding=PLAIN, v2_compressed=True,
                 index_mode="mixed", stats=None, dict_values=None, plain_from_page=None):
        self.__dict__.update(locals())
        del self.__dict__["self"]


def build_chunk(leaf, offset):
    """Return (bytes, ColumnMetaData Struct) for one column chunk at offset."""
    out = bytearray()
    nlev = len(leaf.def_levels) if leaf.def_levels is not None else len(leaf.values)
    defs = leaf.def_levels if leaf.def_levels is not None else [0] * nlev
    reps = leaf.rep_levels if leaf.rep_levels is not None else [0] * nlev
    encodings = set()
    dict_offset = None
    values = leaf.values
    index = None
    total_u = 0
    if leaf.dictionary:
        uniq = list(leaf.dict_values or [])
        pos = {(v if not isinstance(v, float) else struct.pack("<d", v)): i for i, v in enumerate(uniq)}
        index = []
        for v in values:
            key = v if not isinstance(v, float) else struct.pack("<d", v)
            if key not in pos:
                pos[key] = len(uniq)
                uniq.append(v)
            index.append(pos[key])
        body = plain(leaf.ptype, uniq, leaf.type_length)
        comp = compress(leaf.codec, body)
        hdr = page_header(2, len(body), len(comp),
                          (7, "struct", Struct((1, "i32", len(uniq)), (2, "i32", leaf.dict_page_encoding))))
        dict_offset = offset
        out += hdr + comp
        total_u += len(hdr) + len(body)
        encodings.add(leaf.dict_page_encoding)
    data_offset = offset + len(out)
    # Split levels into pages by row (rep level 0 starts a row)
    starts = [k for k in range(nlev) if reps[k] == 0]
    page_rows = leaf.page_rows or max(len(starts), 1)
    bounds = []
    for p in range(0, max(len(starts), 1), page_rows):
        a = starts[p] if starts else 0
        b = starts[p + page_rows] if p + page_rows < len(starts) else nlev
        bounds.append((a, b))
    vpos = 0
    for pno, (a, b) in enumerate(bounds):
        pdefs = defs[a:b]
        nvals = sum(1 for d in pdefs if d == leaf.max_def)
        # a writer whose dictionary grew too large falls back to PLAIN pages
        fallback = leaf.plain_from_page is not None and pno >= leaf.plain_from_page
        if leaf.dictionary and not fallback:
            idx = index[vpos:vpos + nvals]
            width = max(max(idx, default=0).bit_length(), 0)
            vbytes = bytes([width]) + rle_hybrid(idx, width, leaf.index_mode)
            enc = leaf.dict_encoding
        else:
            enc = PLAIN if fallback else leaf.encoding
            vbytes = encode_values(leaf.ptype, values[vpos:vpos + nvals], enc, leaf.type_length)
        encodings.add(enc)
        vpos += nvals
        nrows = sum(1 for r in reps[a:b] if r == 0)
        if leaf.page_version == 1:
            lev = levels_v1(reps[a:b], leaf.max_rep, leaf.level_encoding) + \
                levels_v1(pdefs, leaf.max_def, leaf.level_encoding)
            body = lev + vbytes
            comp = compress(leaf.codec, body)
            dph = Struct((1, "i32", b - a), (2, "i32", enc), (3, "i32", leaf.level_encoding),
                         (4, "i32", leaf.level_encoding), (5, "struct", leaf.stats))
            hdr = page_header(0, len(body), len(comp), (5, "struct", dph))
            out += hdr + comp
            total_u += len(hdr) + len(body)
            if leaf.max_def or leaf.max_rep:
                encodings.add(leaf.level_encoding)
        else:
            rl = rle_hybrid(reps[a:b], leaf.max_rep.bit_length()) if leaf.max_rep else b""
            dl = rle_hybrid(pdefs, leaf.max_def.bit_length()) if leaf.max_def else b""
            comp = compress(leaf.codec, vbytes) if leaf.v2_compressed else vbytes
            nnull = sum(1 for d in pdefs if d < leaf.max_def)
            dph = Struct((1, "i32", b - a), (2, "i32", nnull), (3, "i32", nrows), (4, "i32", enc),
                         (5, "i32", len(dl)), (6, "i32", len(rl)),
                         (7, "bool", None if leaf.v2_compressed else False),
                         (8, "struct", leaf.stats))
            usize = len(rl) + len(dl) + len(vbytes)
            hdr = page_header(3, usize, len(rl) + len(dl) + len(comp), (8, "struct", dph))
            out += hdr + rl + dl + comp
            total_u += len(hdr) + usize
            if leaf.max_def or leaf.max_rep:
                encodings.add(RLE)
    meta = Struct(
        (1, "i32", leaf.ptype),
        (2, "list:i32", sorted(encodings)),
        (3, "list:bin", list(leaf.path)),
        (4, "i32", leaf.codec),
        (5, "i64", nlev),
        (6, "i64", total_u),
        (7, "i64", len(out)),
        (9, "i64", data_offset),
        (11, "i64", dict_offset),
        (12, "struct", leaf.stats),
    )
    return bytes(out), meta


def write_file(path, schema, row_groups, key_value=None, created_by="ctools parquet_writer"):
    """schema: list of schema_element Structs (root first). row_groups: list of
    (num_rows, [Leaf, ...]) in leaf order."""
    out = bytearray(b"PAR1")
    groups = []
    total_rows = 0
    for nrows, leaves in row_groups:
        chunks = []
        start = len(out)
        for leaf in leaves:
            off = len(out)
            data, meta = build_chunk(leaf, off)
            out += data
            chunks.append(Struct((2, "i64", off), (3, "struct", meta)))
        size = len(out) - start
        groups.append(Struct((1, "list:struct", chunks), (2, "i64", size), (3, "i64", nrows),
                             (5, "i64", start), (6, "i64", size), (7, "i16", len(groups))))
        total_rows += nrows
    kv = None
    if key_value:
        kv = [Struct((1, "bin", k), (2, "bin", v)) for k, v in key_value]
    meta = Struct((1, "i32", 1), (2, "list:struct", schema), (3, "i64", total_rows),
                  (4, "list:struct", groups), (5, "list:struct", kv), (6, "bin", created_by))
    footer = encode_struct(meta)
    out += footer + struct.pack("<I", len(footer)) + b"PAR1"
    with open(path, "wb") as f:
        f.write(out)
    return bytes(out)


# ---------------------------------------------------------------------------
# Flat tables
# ---------------------------------------------------------------------------
class Column:
    """A flat column: values may contain None (requires optional)."""

    def __init__(self, name, ptype, values, optional=True, converted=None, logical_type=None,
                 type_length=None, scale=None, precision=None, **leaf_options):
        self.__dict__.update(locals())
        del self.__dict__["self"]


def write_table(path, columns, row_group_rows=None, key_value=None, created_by="ctools parquet_writer"):
    n = len(columns[0].values) if columns else 0
    schema = [schema_element("schema", children=len(columns))]
    for c in columns:
        schema.append(schema_element(c.name, c.ptype, OPTIONAL if c.optional else REQUIRED,
                                     converted=c.converted, logical_type=c.logical_type,
                                     type_length=c.type_length, scale=c.scale,
                                     precision=c.precision))
    groups = []
    step = row_group_rows or max(n, 1)
    for start in range(0, max(n, 1), step):
        stop = min(start + step, n)
        leaves = []
        for c in columns:
            vals = c.values[start:stop]
            if c.optional:
                defs = [0 if v is None else 1 for v in vals]
                present = [v for v in vals if v is not None]
                leaves.append(Leaf(c.ptype, present, defs, None, 1, 0, (c.name,), c.type_length,
                                   **c.leaf_options))
            else:
                assert all(v is not None for v in vals)
                leaves.append(Leaf(c.ptype, list(vals), None, None, 0, 0, (c.name,),
                                   c.type_length, **c.leaf_options))
        groups.append((stop - start, leaves))
        if n == 0:
            break
    return write_file(path, schema, groups, key_value, created_by)


# ---------------------------------------------------------------------------
# ARROW:schema metadata (an Arrow IPC Schema message, base64 encoded)
# ---------------------------------------------------------------------------
class FBTable:
    """A flatbuffer table: fields as (slot, kind, value); kind is one of
    u8, bool, i16, i32, i64, str, table, vec (a list of FBTable)."""

    def __init__(self, *fields):
        self.fields = [f for f in fields if f[2] is not None]


_FB_SIZE = {"u8": 1, "bool": 1, "i16": 2, "i32": 4, "i64": 8, "str": 4, "table": 4, "vec": 4}
_FB_FMT = {"u8": "<B", "bool": "<B", "i16": "<h", "i32": "<i", "i64": "<q"}


def _fb_write(buf, obj):
    """Appends obj (and everything it refers to, after it); returns its offset."""
    while len(buf) % 8:
        buf.append(0)
    if isinstance(obj, (str, bytes)):
        data = obj.encode() if isinstance(obj, str) else obj
        pos = len(buf)
        buf += struct.pack("<I", len(data)) + data + b"\0"
        return pos
    if isinstance(obj, list):
        pos = len(buf)
        buf += struct.pack("<I", len(obj)) + b"\0\0\0\0" * len(obj)
        for i, item in enumerate(obj):
            at = pos + 4 + 4 * i
            child = _fb_write(buf, item)
            struct.pack_into("<I", buf, at, child - at)
        return pos
    nslots = max([f[0] for f in obj.fields], default=-1) + 1
    # inline layout: soffset, then fields largest first
    order = sorted(obj.fields, key=lambda f: -_FB_SIZE[f[1]])
    layout, size = {}, 4
    for slot, kind, _ in order:
        s = _FB_SIZE[kind]
        while size % s:
            size += 1
        layout[slot] = size
        size += s
    while size % 4:
        size += 1
    vt = struct.pack("<HH", 4 + 2 * nslots, size) + b"".join(
        struct.pack("<H", layout.get(i, 0)) for i in range(nslots))
    # the table must start 8-aligned (i64 fields) right after its vtable
    pad = (-(len(buf) + len(vt))) % 8
    buf += b"\0" * pad
    vpos = len(buf)
    buf += vt
    tpos = len(buf)
    buf += b"\0" * size
    struct.pack_into("<i", buf, tpos, tpos - vpos)
    for slot, kind, value in obj.fields:
        at = tpos + layout[slot]
        if kind in _FB_FMT:
            struct.pack_into(_FB_FMT[kind], buf, at, int(value))
    for slot, kind, value in obj.fields:
        if kind in ("str", "table", "vec"):
            at = tpos + layout[slot]
            child = _fb_write(buf, value)
            struct.pack_into("<I", buf, at, child - at)
    return tpos


def fb_buffer(root):
    buf = bytearray(8)
    pos = _fb_write(buf, root)
    struct.pack_into("<I", buf, 0, pos)
    return bytes(buf)


ARROW_TYPES = dict(Null=1, Int=2, FloatingPoint=3, Binary=4, Utf8=5, Bool=6, Decimal=7, Date=8,
                   Time=9, Timestamp=10, Interval=11, List=12, Struct_=13, FixedSizeBinary=15,
                   Map=17, Duration=18, LargeBinary=19, LargeUtf8=20, LargeList=21,
                   BinaryView=23, Utf8View=24)


def arrow_int(bits, signed=True):
    return FBTable((0, "i32", bits), (1, "bool", signed))


def arrow_field(name, type_name, type_table=None, nullable=True, dictionary=None, children=(),
                dict_id=0):
    """dictionary: (index bits, signed, ordered) for a dictionary-encoded field;
    dict_id: its dictionary id (Arrow numbers them 0, 1, ... across fields)."""
    dic = None
    if dictionary:
        bits, signed, ordered = dictionary
        dic = FBTable((0, "i64", dict_id), (1, "table", arrow_int(bits, signed)), (2, "bool", ordered))
    return FBTable((0, "str", name), (1, "bool", nullable), (2, "u8", ARROW_TYPES[type_name]),
                   (3, "table", type_table or FBTable()), (4, "table", dic),
                   (5, "vec", list(children)))


def arrow_schema(fields):
    """The ARROW:schema key-value pair for a list of arrow_field tables."""
    import base64

    schema = FBTable((0, "i16", 0), (1, "vec", list(fields)))
    message = FBTable((0, "i16", 4), (1, "u8", 1), (2, "table", schema), (3, "i64", 0))
    body = fb_buffer(message)
    body += b"\0" * ((-len(body)) % 8)
    framed = struct.pack("<Ii", 0xFFFFFFFF, len(body)) + body
    return ("ARROW:schema", base64.b64encode(framed).decode())
