"""Run native delimited option comparisons through the stata shell alias."""
import re
import shlex
import subprocess
from pathlib import Path
ROOT = Path(__file__).resolve().parents[1]

def fixtures(directory):
    text = 'ID,Text\n1,café Ω 日本語\n2,é\n'
    for name, encoding in [('utf32le','utf-32le'),('utf32be','utf-32be'),('utf16be','utf-16be')]:
        (directory/(name+'.csv')).write_bytes(text.encode(encoding))
    (directory/'utf32bom.csv').write_bytes(text.encode('utf-32'))
    for name,encoding,value in [('cp1251','cp1251','Привет'),('sjis','shift_jis','日本語'),('latin2','iso8859_2','Łódź')]:
        (directory/(name+'.csv')).write_bytes(('ID,Text\n1,'+value+'\n2,abc\n').encode(encoding))
    (directory/'multidelim.txt').write_text('1||2;3\n4;5||6\n7;"8|9";10\n')
    (directory/'multiline.csv').write_text('id,text\n1,"a\nb\nc"\n2,d\n')
    (directory/'range.csv').write_text('id,value\n1,textual\n2,1.5\n3,2.25\n4,99\n')
    (directory/'empty.csv').write_bytes(b'')
    (directory/'words.txt').write_text('idtabvalue\n1tab2\n3tab4\n')
    for name, separator in [('tab','\t'),('semi',';'),('pipe','|'),('colon',':'),('space',' ')]:
        (directory/('autodelim_'+name+'.csv')).write_text(separator.join(('id','value'))+'\n'+separator.join(('1','2'))+'\n'+separator.join(('3','4'))+'\n')
    for width in (20,40,50,100,200,1000):
        for pct in (1,28,29,34,35,50,63,64,65,66,67,70,90,93,94,95):
            (directory/f'storagew{width}p{pct}.csv').write_text('ID,Text\n'+''.join(
                f'{i},'+('x'*(width-3)+f'{i:03d}' if i<=pct else f's{i}')+'\n' for i in range(1,101)))
    numbers = {
        'us':'1,234.56', 'german':'1.234,56', 'french':'1\u202f234,56',
        'french_nbsp':'1\u00a0234,56', 'russian':'1\u00a0234,56', 'swiss':'1’234.56',
        'arabic':'١٬٢٣٤٫٥٦', 'arabic_minus':'\u061c-١٢٣٤٫٥٦',
        'persian':'۱۲۳۴٫۵۶', 'persian_minus':'\u200e−۱۲۳۴٫۵۶', 'hindi':'१२३४.५६',
        'bad_sign':'+1,234.56', 'scientific':'1.23e3', 'infinity':'∞',
        'bad_groups':'1.2,3', 'groups':'1,,23', 'trailing_group':'123,',
    }
    for name,value in numbers.items():
        (directory/(name+'.csv')).write_text('id;value\n1;'+value+'\n2;1\n')
    (directory/'quote_numbers.csv').write_text('id,value\n1,"123"\n2,1""2\n3,12"\n4,"34\n')
    (directory/'quotes.csv').write_text('id,text\n1,a"b"c\n2,"a""b"\n3,c"\n4,"d\n')
    (directory/'long.csv').write_text('id,text\n1,"'+('éΩ""' * 3000)+'"\n2,"'+('x' * 5000)+'"\n')
    regional = [
        ('western','cp1252','La France et les États-Unis ont une économie européenne. '),
        ('central','iso8859_2','Polska gospodarka i społeczeństwo rozwijają się. '),
        ('cyrillic','cp1251','Российская экономика и международная торговля. '),
        ('koi8','koi8_r','Российская экономика и международная торговля. '),
        ('arabic','cp1256','الاقتصاد الدولي والتجارة العربية في العالم. '),
        ('greek','iso8859_7','Η οικονομία και το διεθνές εμπόριο στην Ελλάδα. '),
        ('hebrew','iso8859_8','כלכלה בינלאומית ומסחר בישראל. '),
        ('turkish','iso8859_9','Türkiye ekonomisi ve uluslararası ticaret. '),
        ('japanese','shift_jis','日本の経済と世界の貿易についての研究。'),
        ('eucjp','euc_jp','日本の経済と世界の貿易についての研究。'),
        ('korean','euc_kr','한국의 경제와 국제 무역에 대한 연구. '),
        ('chinese','gb18030','中国经济与国际贸易的研究。'),
        ('traditional','big5','國際經濟與貿易的研究。'),
        ('iso2022','iso2022_jp','日本の経済と世界の貿易についての研究。'),
    ]
    for name, encoding, value in regional:
        (directory/('detect_'+name+'.csv')).write_bytes(('id,text\n1,"'+value*100+'"\n2,abc\n').encode(encoding))

def main():
    directory = ROOT / 'temp/io_delimited_options'
    directory.mkdir(parents=True, exist_ok=True)
    fixtures(directory)
    log = directory / 'regressions.log'
    log.unlink(missing_ok=True)
    driver = directory / 'driver.do'
    driver.write_text(f'''clear all
capture log close _all
log using "{log}", text replace
capture noisily do "{ROOT}/validation/validate_io_delimited_options.do"
local rc = _rc
di "CSV_OPTIONS_DRIVER_RC=`rc'"
log close
exit, clear
''')
    subprocess.run(['/bin/zsh', '-lic', 'stata -q -b do '+shlex.quote(str(driver))],
                   cwd=ROOT, check=True)
    text = log.read_text(errors='replace')
    result = re.findall(r'CSV_OPTIONS_DRIVER_RC=([0-9]+)', text)
    summaries = re.findall(r'CSV_OPTIONS_PASSED=([0-9]+) FAILED=([0-9]+)', text)
    if not result or result[-1] != '0' or not summaries or summaries[-1][1] != '0':
        raise SystemExit('FAIL delimited option regressions; see '+str(log))
    print('PASS '+summaries[-1][0]+' exact option comparisons; '+str(log))

if __name__ == '__main__':
    main()
