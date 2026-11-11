"""Writes the Parquet fixtures of the s2_parquet validation case.

    python validation/support/parquet_fixtures.py [directory]

The default directory is validation/fixtures/parquet. Every file is made by
parquet_writer.py (no Parquet library is needed), deterministically: the
same Python, zlib and libzstd versions give the same bytes. Python 3.14 or
later is required for its compression.zstd module. The brotli_*.parquet and
corrupt_brotli.parquet files also need pyarrow, for libbrotli (and to write
two of them with Arrow's own writer); without it they are skipped.

The files cover each Parquet physical and logical type, nulls, the values
that change the Stata storage type import parquet chooses, row groups and
row ranges, every encoding and codec, nested schemas, column names Stata
must rewrite, Arrow's ARROW:schema metadata (dictionary columns, large and
duration types), and damaged files.
"""

import random
import struct
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import parquet_writer as pw  # noqa: E402
from parquet_writer import (  # noqa: E402
    BIT_PACKED, BOOLEAN, BROTLI, BYTE_ARRAY, BYTE_STREAM_SPLIT, CODECS, DELTA_BINARY_PACKED,
    DELTA_BYTE_ARRAY, DELTA_LENGTH_BYTE_ARRAY, DOUBLE, FIXED_LEN_BYTE_ARRAY, FLOAT, GZIP,
    INT32, INT64, INT96, LZO, OPTIONAL, PLAIN_DICTIONARY, REPEATED, REQUIRED, RLE, SNAPPY,
    ZSTD, Column, FBTable, Leaf, arrow_field, arrow_int, arrow_schema, logical,
    schema_element, write_file, write_table,
)

STR = logical("STRING")


def I(bits, signed=True):
    return logical("INTEGER", bits=bits, signed=signed)


def DEC(scale, precision):
    return dict(logical_type=logical("DECIMAL", scale=scale, precision=precision),
                converted="DECIMAL", scale=scale, precision=precision)


def TS(unit, utc=True):
    return logical("TIMESTAMP", utc=utc, unit=unit)


def TM(unit, utc=True):
    return logical("TIME", utc=utc, unit=unit)


def types(d):
    """Every flat type, with nulls"""
    cols = [
        Column("t_bool", BOOLEAN, [True, False, None, True]),
        Column("t_int8", INT32, [1, -5, None, 100], logical_type=I(8), converted="INT_8"),
        Column("t_int16", INT32, [1, -300, None, 32740], logical_type=I(16)),
        Column("t_int32", INT32, [1, -70000, None, 2147483620]),
        Column("t_int32req", INT32, [7, 8, 9, 10], optional=False),
        Column("t_int64", INT64, [1, -2147483647, None, 2147483620]),
        Column("t_uint8", INT32, [0, 7, None, 100], logical_type=I(8, False), converted="UINT_8"),
        Column("t_uint16", INT32, [0, 7, None, 32740], logical_type=I(16, False)),
        Column("t_uint32", INT32, [0, 7, None, 5], logical_type=I(32, False), converted="UINT_32"),
        Column("t_uint64", INT64, [0, 7, None, 5], logical_type=I(64, False), converted="UINT_64"),
        Column("t_float", FLOAT, [1.5, 0.1, None, -2.0]),
        Column("t_double", DOUBLE, [1.25, 1e300, None, -0.1]),
        Column("t_half", FIXED_LEN_BYTE_ARRAY, [struct.pack("<e", 1.5), struct.pack("<e", -2.25),
                                              None, struct.pack("<e", 65504)],
               type_length=2, logical_type=logical("FLOAT16")),
        Column("t_utf8", BYTE_ARRAY, [b"a", "héllo".encode(), None, b""], logical_type=STR,
               converted="UTF8"),
        Column("t_utf8ct", BYTE_ARRAY, [b"abc", b"x", None, b"yz"], converted="UTF8"),
        Column("t_json", BYTE_ARRAY, [b'{"a":1}', b"[]", None, b"0"], logical_type=logical("JSON")),
        Column("t_enum", BYTE_ARRAY, [b"red", b"blue", None, b"x"], logical_type=logical("ENUM")),
        Column("t_bson", BYTE_ARRAY, [b"\x05\x00\x00\x00\x00", b"x", None, b"y"],
               logical_type=logical("BSON")),
        Column("t_binary", BYTE_ARRAY, [b"\x00\x01z", b"ab", None, b"q"]),
        Column("t_flba", FIXED_LEN_BYTE_ARRAY, [b"abcd", b"\x00\x01\x02\x03", None, b"wxyz"],
               type_length=4),
        Column("t_uuid", FIXED_LEN_BYTE_ARRAY, [bytes(range(16)), bytes(16), None, b"u" * 16],
               type_length=16, logical_type=logical("UUID")),
        Column("t_interval", FIXED_LEN_BYTE_ARRAY, [struct.pack("<III", 1, 2, 3)] * 2 + [None] +
               [bytes(12)], type_length=12, converted="INTERVAL"),
        Column("t_date", INT32, [0, 18000, None, -3653], logical_type=logical("DATE"),
               converted="DATE"),
        Column("t_time_ms", INT32, [0, 3723004, None, 86399999], logical_type=TM("ms"),
               converted="TIME_MILLIS"),
        Column("t_time_us", INT64, [0, 3723004005, None, -1001], logical_type=TM("us", False)),
        Column("t_time_ns", INT64, [0, 3723004005006, None, 1], logical_type=TM("ns")),
        Column("t_ts_ms", INT64, [0, 1700000000123, None, -1], logical_type=TS("ms"),
               converted="TIMESTAMP_MILLIS"),
        Column("t_ts_ms_local", INT64, [0, 1700000000123, None, 5], logical_type=TS("ms", False)),
        Column("t_ts_us", INT64, [0, 1700000000123456, None, -1001], logical_type=TS("us"),
               converted="TIMESTAMP_MICROS"),
        Column("t_ts_ns", INT64, [0, 1700000000123456789, None, -1000001], logical_type=TS("ns", False)),
        Column("t_ts_converted", INT64, [0, 1700000000123, None, 1], converted="TIMESTAMP_MICROS"),
        Column("t_int96", INT96, [(0, 2440588), (3723004005006, 2460000), None, (-1, 2440588)]),
        Column("t_dec_int32", INT32, [12345, -5, None, 0], **DEC(2, 9)),
        Column("t_dec_int64", INT64, [123456789012, -5, None, 1], **DEC(3, 18)),
        Column("t_dec_flba", FIXED_LEN_BYTE_ARRAY, [(12345).to_bytes(5, "big", signed=True),
                                                  (-7).to_bytes(5, "big", signed=True), None,
                                                  bytes(5)], type_length=5, **DEC(4, 10)),
        Column("t_dec_bytes", BYTE_ARRAY, [(12345).to_bytes(3, "big", signed=True), b"\xf9", None,
                                         b"\x00"], **DEC(1, 20)),
        Column("t_dec_wide", FIXED_LEN_BYTE_ARRAY, [(10 ** 40).to_bytes(20, "big", signed=True),
                                                  (-1).to_bytes(20, "big", signed=True), None,
                                                  bytes(20)], type_length=20, **DEC(3, 45)),
        Column("t_null", INT32, [None] * 4, logical_type=logical("UNKNOWN")),
        Column("t_long_text", BYTE_ARRAY, [b"x" * 2046, b"y", None, b"z" * 10], logical_type=STR),
        Column("t_text2045", BYTE_ARRAY, [b"x" * 2045, b"y", None, b""], logical_type=STR),
        Column("t_text_nul", BYTE_ARRAY, [b"a\x00b", b"c", None, b"\xff\xfe"], logical_type=STR),
        Column("t_all_null_text", BYTE_ARRAY, [None] * 4, logical_type=STR),
        Column("t_all_null_num", INT64, [None] * 4),
    ]
    write_table(d / "types.parquet", cols)


def widen(d):
    """Values that change the storage type, or do not fit it"""
    cols = [
        Column("i8_101", INT32, [-127, 101], logical_type=I(8)),
        Column("i8_min", INT32, [-128, 0], logical_type=I(8)),
        Column("i8_wrap", INT32, [1, 1000], logical_type=I(8)),
        Column("u8_255", INT32, [0, 255], logical_type=I(8, False)),
        Column("i16_big", INT32, [-32767, 32741], logical_type=I(16)),
        Column("i16_min", INT32, [-32768, 0], logical_type=I(16)),
        Column("u16_big", INT32, [0, 65535], logical_type=I(16, False)),
        Column("i32_big", INT32, [2147483621, 0]),
        Column("i32_min", INT32, [-2147483648, 0]),
        Column("u32_big", INT32, [0, -1], logical_type=I(32, False)),
        Column("i64_big", INT64, [2147483621, -5]),
        Column("i64_min", INT64, [-2147483648, 7]),
        Column("i64_extremes", INT64, [2 ** 63 - 1, -2 ** 63]),
        Column("u64_big", INT64, [2 ** 31, -1], logical_type=I(64, False)),
        Column("float_big", FLOAT, [1.5, 3.4e38]),
        Column("float_neg", FLOAT, [-3.4e38, 0.0]),
        Column("float_inf", FLOAT, [float("inf"), 1.0]),
        Column("float_nan", FLOAT, [float("nan"), 1.0]),
        Column("double_odd", DOUBLE, [9e307, -0.0]),
        Column("double_inf", DOUBLE, [float("-inf"), float("nan")]),
        Column("double_tiny", DOUBLE, [1e-320, 8.98846567431158e307]),
        Column("date_far", INT32, [2 ** 31 - 1, -2 ** 31], logical_type=logical("DATE")),
        Column("ts_far", INT64, [2 ** 63 - 1, -2 ** 63], logical_type=TS("ms", False)),
        Column("int96_odd", INT96, [(0, 0), (5000000, 0)]),
    ]
    write_table(d / "widen.parquet", cols)


def rowgroups(d):
    """Types follow the values in the row range; text widths the whole column"""
    lens = [1, 2, 3, 4, 5, 6, 7, 20, 9, 10]
    cols = [
        Column("id", INT32, list(range(1, 11))),
        Column("text", BYTE_ARRAY, [b"x" * k for k in lens], logical_type=STR),
        Column("u16", INT32, [1] * 7 + [40000, 1, 1], logical_type=I(16, False)),
        Column("i64", INT64, [1] * 7 + [2 ** 40, 1, 1]),
        Column("f", FLOAT, [1.0] * 7 + [3e38, 1.0, 1.0]),
        Column("i8", INT32, [1, 2, 3, 127, 5, 6, 7, 8, 9, -128], logical_type=I(8)),
    ]
    write_table(d / "rowgroups.parquet", cols, row_group_rows=3)


def encodings(d):
    random.seed(20261104)
    n = 300
    ints = [random.randint(-10 ** 6, 10 ** 6) if i % 7 else None for i in range(n)]
    big = [random.randint(-2 ** 62, 2 ** 62) for i in range(n)]
    words = [(b"w%d" % (i * 7 % 50)) if i % 11 else None for i in range(n)]
    flt = [i / 3 for i in range(n)]
    cols = [
        Column("plain", INT32, ints),
        Column("dict_rle", INT32, [v % 9 if v is not None else None for v in ints], dictionary=True,
               index_mode="rle"),
        Column("dict_packed", BYTE_ARRAY, words, logical_type=STR, dictionary=True,
               index_mode="packed", page_rows=70),
        Column("dict_plain_dictionary", BYTE_ARRAY, words, logical_type=STR, dictionary=True,
               dict_encoding=PLAIN_DICTIONARY, dict_page_encoding=PLAIN_DICTIONARY),
        Column("dict_fallback", BYTE_ARRAY, words, logical_type=STR, dictionary=True,
               page_rows=60, plain_from_page=2),
        Column("delta32", INT32, ints, encoding=DELTA_BINARY_PACKED, page_rows=130),
        Column("delta64", INT64, big, encoding=DELTA_BINARY_PACKED),
        Column("delta_length", BYTE_ARRAY, words, logical_type=STR, encoding=DELTA_LENGTH_BYTE_ARRAY),
        Column("delta_bytes", BYTE_ARRAY, sorted(w for w in words if w) + [None] * words.count(None),
               logical_type=STR, encoding=DELTA_BYTE_ARRAY),
        Column("delta_flba", FIXED_LEN_BYTE_ARRAY, [b"%03d" % (i // 3) for i in range(n)],
               type_length=3, encoding=DELTA_BYTE_ARRAY),
        Column("split_float", FLOAT, flt, encoding=BYTE_STREAM_SPLIT),
        Column("split_double", DOUBLE, [v / 7 for v in flt], encoding=BYTE_STREAM_SPLIT),
        Column("split_int32", INT32, ints, encoding=BYTE_STREAM_SPLIT),
        Column("split_int64", INT64, big, encoding=BYTE_STREAM_SPLIT),
        Column("bool_plain", BOOLEAN, [i % 3 == 0 if i % 5 else None for i in range(n)]),
        Column("bool_rle", BOOLEAN, [i % 3 == 0 if i % 5 else None for i in range(n)], encoding=RLE,
               page_version=2),
        Column("levels_bit_packed", INT32, ints, level_encoding=BIT_PACKED),
        Column("v2_snappy", INT32, ints, page_version=2, codec=SNAPPY, page_rows=100),
        Column("v2_uncompressed_page", BYTE_ARRAY, words, logical_type=STR, page_version=2,
               codec=GZIP, v2_compressed=False),
        Column("v2_dict", BYTE_ARRAY, words, logical_type=STR, page_version=2, dictionary=True,
               codec=ZSTD, page_rows=90),
    ]
    write_table(d / "encodings.parquet", cols, row_group_rows=160)


def codecs(d):
    random.seed(7)
    n = 400
    words = [("w%d" % i) * random.randint(1, 4) for i in range(60)]
    text = [(" ".join(random.choice(words) for _ in range(random.randint(0, 5)))).encode()
            for _ in range(n)]
    x = [round(random.gauss(0, 1), 3) for _ in range(n)]
    for name, codec in CODECS.items():
        write_table(d / f"codec_{name}.parquet", [
            Column("text", BYTE_ARRAY, text, logical_type=STR, codec=codec, page_rows=150),
            Column("x", DOUBLE, x, codec=codec),
            Column("k", INT32, [i % 17 for i in range(n)], codec=codec, dictionary=True),
        ])
    # LZO: no compressor is needed, as no reader decodes it
    orig = pw.compress
    pw.compress = lambda codec, data: data if codec == LZO else orig(codec, data)
    write_table(d / "codec_lzo.parquet", [Column("a", INT32, [1, 2], codec=LZO)])
    pw.compress = orig
    # Zstd at a high level: Huffman literals, FSE tables, several blocks
    from compression import zstd

    random.seed(9)
    vocab = [("v%d" % i) * random.randint(1, 3) for i in range(3000)]
    long_text = [(" ".join(random.choice(vocab) for _ in range(random.randint(0, 6)))).encode()
                 for _ in range(6000)]
    pw.compress = lambda codec, data: (zstd.compress(data, level=19) if codec == ZSTD
                                       else orig(codec, data))
    write_table(d / "zstd_level19.parquet", [
        Column("text", BYTE_ARRAY, long_text, logical_type=STR, codec=ZSTD),
        Column("k", INT64, [i * 1000 for i in range(6000)], codec=ZSTD,
               encoding=DELTA_BINARY_PACKED)])
    pw.compress = orig


ENGLISH = ("the people of the world have been working with data for a long time and "
           "they know that the information in a table is more useful when it is organized "
           "by column; this is one of the reasons that the format was designed, and the "
           "first version of the specification describes how pages are compressed and "
           "encoded. Please contact the members of our community for more details about "
           "the history of this project, which has been available since 2013.").split()


def brotli(d):
    """Brotli pages compressed by libbrotli (through pyarrow) at qualities 0-11:
    every physical type, nulls, dictionary pages, data pages v1 and v2, several
    pages and row groups; files written by pyarrow itself; and truncated
    Brotli pages. pyarrow is used only for these files."""
    try:
        import pyarrow as pa
        import pyarrow.parquet as pq
    except ImportError:
        print("pyarrow is not installed: the brotli_*.parquet fixtures were not written")
        return
    random.seed(7932)
    n = 2500
    calls = iter(range(10 ** 9))
    qualities = (11, 0, 5, 1, 9, 4, 10, 2, 6, 3, 7, 8)
    orig = pw.compress

    def libbrotli(codec, data):
        if codec != BROTLI:
            return orig(codec, data)
        q = qualities[next(calls) % len(qualities)]
        return pa.Codec("brotli", compression_level=q).compress(data, asbytes=True)

    text = [" ".join(random.choice(ENGLISH) for _ in range(random.randint(0, 40))).encode()
            if i % 9 else None for i in range(n)]
    cats = [random.choice([b"north", b"south", b"east", b"west", None]) for _ in range(n)]
    ints = [random.randint(-10 ** 6, 10 ** 6) if i % 11 else None for i in range(n)]
    pw.compress = libbrotli
    try:
        cols = [
            Column("b_bool", BOOLEAN, [i % 3 == 0 if i % 7 else None for i in range(n)],
                   codec=BROTLI),
            Column("b_int32", INT32, ints, codec=BROTLI, page_rows=700),
            Column("b_int32_v2", INT32, ints, codec=BROTLI, page_version=2, page_rows=500),
            Column("b_int64", INT64, [random.randint(-2 ** 60, 2 ** 60) for _ in range(n)],
                   codec=BROTLI, optional=False),
            Column("b_int96", INT96, [(random.randint(0, 86399 * 10 ** 9), 2460000 + i % 400)
                                      if i % 13 else None for i in range(n)], codec=BROTLI),
            Column("b_float", FLOAT, [round(random.gauss(0, 10), 2) for _ in range(n)],
                   codec=BROTLI),
            Column("b_double", DOUBLE, [round(random.gauss(0, 1), 4) if i % 5 else None
                                        for i in range(n)], codec=BROTLI, page_rows=600),
            Column("b_text", BYTE_ARRAY, text, logical_type=STR, codec=BROTLI, page_rows=400),
            Column("b_text_v2", BYTE_ARRAY, text, logical_type=STR, codec=BROTLI,
                   page_version=2, page_rows=900),
            Column("b_dict_text", BYTE_ARRAY, cats, logical_type=STR, codec=BROTLI,
                   dictionary=True, page_rows=800),
            Column("b_dict_v2", INT32, [v % 13 if v is not None else None for v in ints],
                   codec=BROTLI, dictionary=True, page_version=2),
            Column("b_flba", FIXED_LEN_BYTE_ARRAY, [b"%06d" % (i * 37 % 1000000)
                                                    for i in range(n)], type_length=6,
                   codec=BROTLI),
            Column("b_delta_ts", INT64, [1700000000000000 + 1000003 * i for i in range(n)],
                   logical_type=TS("us"), encoding=DELTA_BINARY_PACKED, codec=BROTLI),
        ]
        write_table(d / "brotli_types.parquet", cols, row_group_rows=900)
        # Damaged: every Brotli page loses its last bytes
        pw.compress = lambda codec, data: libbrotli(codec, data)[:-3]
        write_table(d / "corrupt_brotli.parquet",
                    [Column("a", INT32, list(range(1000)), codec=BROTLI),
                     Column("s", BYTE_ARRAY, [b"x%d" % i for i in range(1000)],
                            logical_type=STR, codec=BROTLI)])
    finally:
        pw.compress = orig
    # Written by pyarrow itself (Arrow's own page layout and metadata)
    table = pa.table({
        "id": pa.array(range(n), pa.int32()),
        "big": pa.array([i * 1000003 - 5 for i in range(n)], pa.int64()),
        "x": pa.array([None if i % 17 == 0 else round(random.random(), 6) for i in range(n)]),
        "flag": pa.array([None if i % 19 == 0 else i % 2 == 0 for i in range(n)]),
        "words": pa.array([None if t is None else t.decode() for t in text], pa.string()),
        "cat": pa.array([None if c is None else c.decode() for c in cats], pa.string()),
    })
    for version, level in (("1.0", 11), ("2.0", 4)):
        pq.write_table(table, d / f"brotli_pyarrow_v{version[0]}.parquet", compression="brotli",
                       compression_level=level, data_page_version=version,
                       row_group_size=1000, use_dictionary=["cat"])


def nested(d):
    S = schema_element
    sch = [S("schema", children=5), S("id", INT32, REQUIRED),
           S("st", repetition=OPTIONAL, children=2), S("x", INT32, OPTIONAL), S("y", DOUBLE, OPTIONAL),
           S("lst", repetition=OPTIONAL, children=1, converted="LIST", logical_type=logical("LIST")),
           S("list", repetition=REPEATED, children=1), S("element", INT32, OPTIONAL),
           S("m", repetition=OPTIONAL, children=1, converted="MAP", logical_type=logical("MAP")),
           S("key_value", repetition=REPEATED, children=2, converted="MAP_KEY_VALUE"),
           S("key", BYTE_ARRAY, REQUIRED, logical_type=STR), S("value", INT32, OPTIONAL),
           S("rep", INT32, REPEATED)]
    leaves = [Leaf(INT32, [1, 2], path=("id",)),
              Leaf(INT32, [5], [2, 0], max_def=2, path=("st", "x")),
              Leaf(DOUBLE, [1.5], [2, 0], max_def=2, path=("st", "y")),
              Leaf(INT32, [7, 8, 9], [3, 3, 3], [0, 1, 0], max_def=3, max_rep=1,
                   path=("lst", "list", "element")),
              Leaf(BYTE_ARRAY, [b"k1", b"k2"], [2, 2], [0, 0], max_def=2, max_rep=1,
                   path=("m", "key_value", "key")),
              Leaf(INT32, [3, 4], [3, 3], [0, 0], max_def=3, max_rep=1,
                   path=("m", "key_value", "value")),
              Leaf(INT32, [7, 8, 9], [1, 1, 1], [0, 1, 0], max_def=1, max_rep=1, path=("rep",))]
    write_file(d / "nested.parquet", sch, [(2, leaves)])
    sch = [S("schema", children=1), S("st", repetition=OPTIONAL, children=1), S("x", INT32, OPTIONAL)]
    write_file(d / "nested_only.parquet", sch,
               [(2, [Leaf(INT32, [5], [2, 0], max_def=2, path=("st", "x"))])])


def names(d):
    good = ["ok", "has space", "1digit", "x" * 40, "_", "café", "a-b", "", "日本",
            "1", "ó", "Ⅷx", "a²", "١a", "_1a", "A b", "e", "x" * 31 + "éé"]
    write_table(d / "names.parquet", [Column(n, INT32, [i]) for i, n in enumerate(good)])
    write_table(d / "names_dup.parquet", [Column(n, INT32, [i]) for i, n in
                                          enumerate(["a", "a b", "ab", "B"])])
    write_table(d / "names_trunc.parquet", [Column("y" * 33, INT32, [1]), Column("y" * 34, INT32, [2])])
    write_table(d / "names_reserved.parquet", [Column("ok", INT32, [1]), Column("in", INT32, [2])])
    write_table(d / "names_select.parquet", [Column("u", INT32, [None], logical_type=logical("UNKNOWN")),
                                             Column("1", INT32, [5]), Column("ok", INT32, [6]),
                                             Column("b", INT32, [7])])


def arrow(d):
    ids = iter(range(100))
    dic = lambda name, kind="Utf8", bits=32: arrow_field(name, kind, dictionary=(bits, True, False),
                                                         dict_id=next(ids))
    cats = [b"lo", b"mid", b"hi", None, b"mid", b"lo"]
    cols = [
        Column("cat", BYTE_ARRAY, cats, logical_type=STR, dictionary=True),
        Column("cat8", BYTE_ARRAY, cats, logical_type=STR, dictionary=True),
        Column("cat_plain", BYTE_ARRAY, [b"q", b"p", b"q", b"r", None, b"p"], logical_type=STR),
        Column("cat_empty", BYTE_ARRAY, [b"", b"a", b"", None, b"a", b"b"], logical_type=STR,
               dictionary=True),
        Column("cat_unused", BYTE_ARRAY, [b"b", b"a", b"b", b"b", None, b"a"], logical_type=STR,
               dictionary=True, dict_values=[b"zz", b"a", b"b", b"an unused, wider entry"]),
        Column("cat_null", BYTE_ARRAY, [None] * 6, logical_type=STR, dictionary=True),
        Column("cat_binary", BYTE_ARRAY, [b"ab", b"cd", b"ab", None, b"\x00x", b"cd"], dictionary=True),
        Column("large_text", BYTE_ARRAY, [b"abc", None, b"d", b"", b"e", b"f"], logical_type=STR),
        Column("large_bytes", BYTE_ARRAY, [b"q", b"r", None, b"s", b"t", b"u"]),
        Column("duration", INT64, [1500000, -1, None, 0, 1, 2]),
        Column("plain_int", INT32, [1, 2, 3, 4, 5, 6]),
    ]
    fields = [dic("cat"), dic("cat8", bits=8), dic("cat_plain"), dic("cat_empty"), dic("cat_unused"),
              dic("cat_null"), dic("cat_binary", "Binary"), arrow_field("large_text", "LargeUtf8"),
              arrow_field("large_bytes", "LargeBinary"),
              arrow_field("duration", "Duration", FBTable((0, "i16", 2))),
              arrow_field("plain_int", "Int", arrow_int(32))]
    write_table(d / "arrow.parquet", cols, key_value=[arrow_schema(fields)])
    # One dictionary per row group, unified in order
    write_table(d / "arrow_rowgroups.parquet",
                [Column("cat", BYTE_ARRAY, [b"a", b"b", b"c", b"z", b"a", b"e", b"q", b"r", b"s"],
                        logical_type=STR, dictionary=True)],
                row_group_rows=3, key_value=[arrow_schema([dic("cat")])])
    # Past 100 entries a dictionary column is an int
    write_table(d / "arrow_int.parquet",
                [Column("cat", BYTE_ARRAY, [b"c%03d" % i for i in range(101)] + [b"x" * 15],
                        logical_type=STR, dictionary=True)],
                key_value=[arrow_schema([dic("cat")])])
    # Field counts that disagree: Arrow ignores the stored schema
    write_table(d / "arrow_mismatch.parquet",
                [Column("a", BYTE_ARRAY, [b"x", b"y"], logical_type=STR), Column("b", INT32, [1, 2])],
                key_value=[arrow_schema([dic("a")])])
    # Arrow rejects one dictionary id for two value types
    write_table(d / "arrow_conflict.parquet",
                [Column("a", BYTE_ARRAY, [b"x"], logical_type=STR, dictionary=True),
                 Column("b", BYTE_ARRAY, [b"y"], dictionary=True)],
                key_value=[arrow_schema([arrow_field("a", "Utf8", dictionary=(32, True, False)),
                                         arrow_field("b", "Binary", dictionary=(32, True, False))])])
    write_table(d / "arrow_garbage.parquet", [Column("a", INT32, [1, 2])],
                key_value=[("ARROW:schema", "!!!not base64")])


def damaged(d):
    write_table(d / "empty.parquet", [Column("a", INT32, []), Column("s", BYTE_ARRAY, [],
                                                                      logical_type=STR)])
    data = bytearray(write_table(d / "corrupt_snappy.parquet",
                                 [Column("a", INT32, list(range(100)), codec=SNAPPY)]))
    data[40] ^= 0xFF
    data[41] ^= 0xFF
    (d / "corrupt_snappy.parquet").write_bytes(bytes(data))
    data = write_table(d / "truncated.parquet", [Column("a", INT32, list(range(100)))])
    (d / "truncated.parquet").write_bytes(data[:-20])
    (d / "not_parquet.parquet").write_bytes(b"id,value\n1,2\n")


def multi(d):
    random.seed(3)
    n = 3000
    cols = [Column("id", INT32, list(range(n)), optional=False, page_rows=700, codec=ZSTD),
            Column("x", DOUBLE, [round(random.random(), 4) if i % 13 else None for i in range(n)],
                   page_rows=900, codec=SNAPPY),
            Column("s", BYTE_ARRAY, [b"%d" % (i * 7919 % 1000) for i in range(n)], logical_type=STR,
                   dictionary=True, page_rows=600, plain_from_page=1, codec=GZIP),
            Column("t", INT64, [1700000000000 + 997 * i for i in range(n)], logical_type=TS("ms"),
                   encoding=DELTA_BINARY_PACKED, codec=ZSTD)]
    write_table(d / "multi.parquet", cols, row_group_rows=1100)


def main():
    root = Path(__file__).resolve().parents[2]
    out = Path(sys.argv[1]) if len(sys.argv) > 1 else root / "validation" / "fixtures" / "parquet"
    out.mkdir(parents=True, exist_ok=True)
    for make in (types, widen, rowgroups, encodings, codecs, brotli, nested, names, arrow, damaged,
                 multi):
        make(out)
    for p in sorted(out.glob("*.parquet")):
        print(f"{p.stat().st_size:8d} {p.name}")


if __name__ == "__main__":
    main()
