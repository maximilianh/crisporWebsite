import sys
from StringIO import StringIO
import re
import pandas

raw = sys.stdin.read()
array_name_pattern = r"[a-zA-Z_0-9]+\s+ARRAY"
tables = filter(lambda x: x, re.split(array_name_pattern, raw))
table_names = [line for line in raw.split("\n")
               if re.match(array_name_pattern, line)]
assert len(table_names) == len(tables)
tables = map(StringIO, tables)
tables = map(pandas.read_table, tables)
print tables[0]
tables = map(lambda x:x.dropna(axis=1,how='all'), tables)
writer = pandas.ExcelWriter(sys.argv[1])
map(lambda (table,name): table.to_excel(writer, name), zip(tables, table_names))
writer.save()
