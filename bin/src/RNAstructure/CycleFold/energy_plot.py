# coding: utf-8
import pandas as p
import numpy as np
import math
import energy as e
from matplotlib import pyplot as plt

#read the data
data = [x.split() for x in open("xia_duplexes.txt")]
d = p.DataFrame(data)
d.columns = ["sequence", "dG", "dH", 'dS', 'tm', 'dG_calc', 'dH_calc', 'dS_calc', 'tm_calc', 'reference']
d[0] = d.apply(lambda x:x[0].replace('p','').replace('/',''), axis='columns')

#do some claculations
d['length'] = d.apply(lambda x: len(x[0]), axis='columns')
d['ncm_probability'] = d.apply(lambda x:e.efncm_duplex(x[0]), axis=1)

#relationship between ncm probability and length
#with a regression line out to 0
plt.style.use('ggplot')
m, b = np.polyfit(d['length'], map(lambda x:math.log(x), d['ncm_probability']), 1)
plt.plot(np.arange(11), m*np.arange(11) + b, '-')
plt.scatter(d['length'], map(lambda x:math.log(x), d['ncm_probability']))
plt.ylabel("log probability from NCM model")
plt.xlabel("length, nt")
plt.annotate("regression line: y = %.3f * x + %0.3f" % (m,b), xy=(2,5))
plt.savefig("log_probabilities.png")
plt.show()
