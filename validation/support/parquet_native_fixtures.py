"""Full-file native corpus: independent expected cells and hand-framed damage.

Valid files include the repository writer and checked-in PyArrow-produced files.
Corrupt bodies/metadata below are constructed explicitly, not by the C decoder.
"""
import random
import shutil
import struct
from pathlib import Path
from validation import ROOT
from validation.support import parquet_writer as pw


def reference(path, columns):
    with path.open('wb') as f:
        f.write(struct.pack('<II',len(columns[0][1]),len(columns)))
        for kind, values in columns:
            f.write(bytes([kind]))
            for value in values:
                f.write(bytes([value is not None]))
                if value is None:
                    continue
                if kind==2:
                    if isinstance(value,str): value=value.encode()
                    f.write(struct.pack('<I',len(value))+value)
                else:
                    f.write(struct.pack('<q' if kind==1 else '<d',value))


def corpus(folder):
    folder=Path(folder);folder.mkdir(parents=True,exist_ok=True)
    cases=[]
    def add(name, columns=None, external=None, sweep=True):
        path=folder/(name+'.parquet')
        if external: shutil.copyfile(ROOT/'validation/fixtures/parquet'/external,path)
        if columns: reference(path.with_suffix('.ref'),columns)
        cases.append((name,columns is not None,sweep))
        return path
    x=[None if i%4==0 else i-5 for i in range(12)]
    s=[(b'a',b'ab','東京'.encode(),b'',None,b'a\0z')[i%6] for i in range(12)]
    for codec in (pw.UNCOMPRESSED,pw.SNAPPY,pw.GZIP,pw.BROTLI,pw.LZ4,pw.LZ4_RAW):
        for version in (1,2):
            for encoding in ('plain','dict','delta','split'):
                name=f'good_{codec}_{version}_{encoding}'
                opts=dict(codec=codec,page_version=version,page_rows=3)
                xo=dict(dictionary=encoding=='dict',encoding=pw.DELTA_BINARY_PACKED if encoding=='delta' else
                        pw.BYTE_STREAM_SPLIT if encoding=='split' else pw.PLAIN)
                so=dict(dictionary=encoding=='dict',encoding=pw.DELTA_BYTE_ARRAY if encoding=='delta' else
                        pw.DELTA_LENGTH_BYTE_ARRAY if encoding=='split' else pw.PLAIN)
                pw.write_table(add(name,[(0,x),(2,s)]),[
                    pw.Column('x',pw.INT32,x,**opts,**xo),
                    pw.Column('s',pw.BYTE_ARRAY,s,logical_type=pw.logical('STRING'),**opts,**so)],row_group_rows=7)
    # Exact source values of the independently compressed checked-in Zstd file.
    rng=random.Random(7)
    words=[('w%d'%i)*rng.randint(1,4) for i in range(60)]
    text=[' '.join(rng.choice(words) for _ in range(rng.randint(0,5))).encode() for _ in range(400)]
    nums=[round(rng.gauss(0,1),3) for _ in range(400)]
    add('good_zstd',[(2,text),(0,nums),(0,[i%17 for i in range(400)])],'codec_zstd.parquet',False)
    # PyArrow is the producer, including its own metadata/page framing. The
    # first two columns are analytic; all remaining columns are decoded too.
    for version in (1,2):
        add(f'good_arrow{version}',[(0,list(range(2500))),(1,[i*1000003-5 for i in range(2500)])],
            f'brotli_pyarrow_v{version}.parquet',False)

    def raw(name, body, physical=pw.INT32, encoding=pw.PLAIN, count=2, codec=0,
            usize=None,csize=None,dictionary=None,optional=False,groups=None,total=None):
        dph=pw.Struct((1,'i32',count),(2,'i32',encoding),(3,'i32',pw.RLE),(4,'i32',pw.RLE))
        page=pw.page_header(0,len(body) if usize is None else usize,len(body) if csize is None else csize,
                            (5,'struct',dph))+body
        dict_off=None;data_off=4
        if dictionary is not None:
            dp=pw.page_header(2,len(dictionary),len(dictionary),
                (7,'struct',pw.Struct((1,'i32',1),(2,'i32',pw.PLAIN))))+dictionary
            dict_off=4;data_off+=len(dp);page=dp+page
        schema=[pw.schema_element('schema',children=1),
                pw.schema_element('x',physical,pw.OPTIONAL if optional else pw.REQUIRED)]
        metadata=pw.Struct((1,'i32',physical),(2,'list:i32',[encoding]),(3,'list:bin',['x']),
            (4,'i32',codec),(5,'i64',count),(6,'i64',len(page)),(7,'i64',len(page)),
            (9,'i64',data_off),(11,'i64',dict_off))
        rows=groups or [count]
        rgs=[pw.Struct((1,'list:struct',[pw.Struct((2,'i64',4),(3,'struct',metadata))]),
                       (2,'i64',len(page)),(3,'i64',n)) for n in rows]
        footer=pw.encode_struct(pw.Struct((1,'i32',1),(2,'list:struct',schema),
            (3,'i64',sum(rows) if total is None else total),(4,'list:struct',rgs)))
        add(name).write_bytes(b'PAR1'+page+footer+struct.pack('<I',len(footer))+b'PAR1')
    raw('bad_plain',b'\0')
    raw('bad_page_size',bytes(8),csize=2147483647)
    raw('bad_page_allocation',bytes(8),usize=2147483647)
    raw('bad_dictionary',b'\x01\x04\x01',encoding=pw.RLE_DICTIONARY,dictionary=struct.pack('<i',42))
    raw('bad_levels',struct.pack('<I',2)+b'\x04\x02'+bytes(8),optional=True)
    raw('bad_delta_prefix',b'\x80\x01\x04\x00\x00',physical=pw.BYTE_ARRAY,encoding=pw.DELTA_BYTE_ARRAY)
    raw('bad_row_total',bytes(8),groups=[2],total=3)
    raw('bad_row_overflow',bytes(8),groups=[2**63-1,1],total=0)
    raw('bad_negative_rows',bytes(8),groups=[-1],total=0)
    raw('bad_huge_rows',bytes(8),groups=[2**31],total=2**31)
    for codec in (pw.SNAPPY,pw.GZIP,pw.BROTLI,pw.LZ4,pw.LZ4_RAW,pw.ZSTD):
        raw('bad_codec_'+str(codec),b'\0',codec=codec,usize=8)
    good=(folder/'good_0_1_plain.parquet').read_bytes()
    add('bad_magic').write_bytes(b'NOPE'+good[4:])
    add('bad_footer_size').write_bytes(good[:-8]+b'\xff'*4+b'PAR1')
    add('bad_truncated').write_bytes(good[:-9])
    header='static const struct { const char *name; int valid,sweep; } cases[]={\n'
    header+=''.join('{"%s",%d,%d},\n'%(name,valid,sweep) for name,valid,sweep in cases)+'};\n'
    (folder/'cases.h').write_text(header,encoding='utf-8')
    return folder
