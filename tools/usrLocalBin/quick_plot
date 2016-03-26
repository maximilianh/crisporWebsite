#!/usr/bin/env python
"""
quick_plot
9 Oct 2013
Dent Earl (dent.earl a gmail.com)

A quick plotting program for creating fast sketches of data.

"""
##############################
# Copyright (C) 2013 by
# Dent Earl (dearl@soe.ucsc.edu, dent.earl@gmail.com)
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
# THE SOFTWARE.
##############################
# plotting boilerplate / cargo cult
import matplotlib
matplotlib.use('Agg')
#####
# the param pdf.fonttype allows for text to be editable in Illustrator.
# Use either Output Type 3 (Type3) or Type 42 (TrueType)
matplotlib.rcParams['pdf.fonttype'] = 42
import matplotlib.backends.backend_pdf as pltBack
import matplotlib.lines as lines
import matplotlib.patches as patches
import matplotlib.pyplot as plt
import matplotlib.pylab as pylab
from matplotlib.ticker import ScalarFormatter, FormatStrFormatter
import numpy
##############################
from argparse import ArgumentParser
import os
from scipy.stats import scoreatpercentile, linregress
import sys


class BadInput(Exception):
  pass


class Data(object):
  """Class data holds data for plotting
  """
  def __init__(self):
    self.data = None  # this will be an n by 2 numpy array.
    self.label = ''


def InitArguments(parser):
  """Initialize arguments for the program.

  Args:
    parser: an argparse parser object
  """
  parser.add_argument('files', nargs='+', help='files to plot')
  parser.add_argument('--out', dest='out', default='my_plot',
                      type=str,
                      help=('path/filename where figure will be created. No '
                            'extension needed. default=%(default)s'))
  parser.add_argument('--mode', dest='mode', default='line', type=str,
                      help=('plotting mode. may be in (line, scatter, '
                            'column, bar, hist, tick, point) '
                            'default=%(default)s'))
  parser.add_argument('--alpha', default=1.0, type=float,
                      help='alpha value for markers in --mode scatter')
  parser.add_argument('--dot_size', '--markersize', dest='markersize',
                      default=2.0, type=float,
                      help='value for markers in --mode scatter')
  parser.add_argument('--lineWidth', dest='linewidth', default=2.0,
                      type=float,
                      help='Line width for the plot. default=%(default)s')
  parser.add_argument('--logy', dest='is_log_y', default=False,
                      action='store_true',
                      help='Put the y-axis into log. default=%(default)s')
  parser.add_argument('--logx', dest='is_log_x', default=False,
                      action='store_true',
                      help='Put the x-axis into log. default=%(default)s')
  parser.add_argument('--title', dest='title', type=str,
                      default='sentinel_value',
                      help='Plot title.')
  parser.add_argument('--xlabel', dest='xlabel', type=str,
                      default='sentinel_value',
                      help='X-axis label.')
  parser.add_argument('--ylabel', dest='ylabel', type=str,
                      default='sentinel_value',
                      help='Y-axis label.')
  parser.add_argument('--height', dest='height', default=4.0, type=float,
                      help='height of image, in inches. default=%(default)s')
  parser.add_argument('--width', dest='width', default=9.0, type=float,
                      help='width of image, in inches. default=%(default)s')
  parser.add_argument('--dpi', dest='dpi', default=300,
                      type=int,
                      help=('dots per inch of raster outputs, i.e. '
                            'if --outFormat is all or png. '
                            'default=%(default)s'))
  parser.add_argument('--out_format', dest='out_format', default='pdf',
                      type=str,
                      help=('output format [pdf|png|eps|all]. '
                            'default=%(default)s'))
  parser.add_argument('--no_legend', dest='is_legend', default=True,
                      action='store_false',
                      help=('Turns off the filename / color legend. '
                            'Helpful for large numbers of files.'))
  parser.add_argument('--regression', dest='regression', default=False,
                      action='store_true',
                      help='turn on a simple linear regression line')
  parser.add_argument('--jitter', dest='jitter', default=False,
                      action='store_true',
                      help='turn on jitter for certain plotting modes')


def CheckArguments(args, parser):
  """Verify that input arguments are correct and sufficient.

  Args:
    args: an argparse arguments object
    parser: an argparse parser object
  """
  if len(args.files) > 0:
    for f in args.files:
      if not os.path.exists(f):
        parser.error('File %s does not exist.\n' % f)
  else:
    parser.error('File paths must be passed in on command line!')
  if args.dpi < 72:
    parser.error('--dpi %d less than screen res, 72. Must be >= 72.'
                 % args.dpi)
  if args.out_format not in ('pdf', 'png', 'eps', 'all'):
    parser.error('Unrecognized --out_format %s. Choose one from: '
                 'pdf png eps all.' % args.out_format)
  if args.mode not in ('line', 'scatter', 'bar', 'column', 'hist',
                       'tick', 'point'):
    parser.error('Unrecognized --mode %s. Choose one from: '
                 'line scatter bar column hist tick point.' % args.mode)
  if (args.out.endswith('.png') or args.out.endswith('.pdf') or
      args.out.endswith('.eps')):
    args.out = args.out[:-4]
  args.ymax = -sys.maxint
  args.ymin = sys.maxint
  args.xmax = args.ymax
  args.xmin = args.ymin
  # TODO: allow for a way to override the color list
  args.color_list = ['#1f77b4', # d blue
                     '#aec7e8', # l blue
                     '#ff7f0e', # d orange
                     '#ffbb78', # l orange
                     '#2ca02c', # d green
                     '#98df8a', # l green
                     '#d62728', # d red
                     '#ff9896', # l red
                     '#9467bd', # d purple
                     '#c5b0d5', # l purple
                     '#8c564b', # d brown
                     '#c49c94', # l brown
                     '#e377c2', # d lavender
                     '#f7b6d2', # l lavender
                     '#7f7f7f', # d gray
                     '#c7c7c7', # l gray
                     '#bcbd22', # d olive
                     '#dbdb8d', # l olive
                     '#17becf', # d aqua
                     '#9edae5'  # l aqua
                     ]


def InitImage(args):
  """Initialize a new image.

  Args:
    args: an argparse arguments object

  Returns:
    fig: a matplotlib figure object
    pdf: a matplotlib pdf drawing (backend) object
  """
  pdf = None
  if args.out_format == 'pdf' or args.out_format == 'all':
    pdf = pltBack.PdfPages(args.out + '.pdf')
  fig = plt.figure(figsize=(args.width, args.height),
                   dpi=args.dpi, facecolor='w')
  return (fig, pdf)


def EstablishAxes(fig, args):
  """Create a single axis on the figure object

  Args:
    fig: a matplotlib figure object
    args: an argparse arguments object

  Returns:
    ax: a matplotlib axis object
  Raises:
    ValueError: If an unknown spine location is passed.
  """
  args.axLeft = 0.095
  args.axWidth = 0.84
  args.axBottom = 0.17
  args.axHeight = 0.76
  ax = fig.add_axes([args.axLeft, args.axBottom,
                     args.axWidth, args.axHeight])
  ax.yaxis.set_major_locator(pylab.NullLocator())
  ax.xaxis.set_major_locator(pylab.NullLocator())
  for loc, spine in ax.spines.iteritems():
    if loc in ['left', 'bottom']:
      spine.set_position(('outward', 10))
    elif loc in ['right', 'top']:
      spine.set_color('none')
    else:
      raise ValueError('unknown spine location: %s' % loc)
  ax.xaxis.set_ticks_position('bottom')
  ax.yaxis.set_ticks_position('left')
  return ax


def WriteImage(fig, pdf, args):
  """Write the image to disk.

  Args:
    fig: a matplotlib figure object
    pdf: a matplotlib pdf drawing (backend) object
    args: an argparse arguments object
  """
  if args.out_format == 'pdf':
    fig.savefig(pdf, format = 'pdf')
    pdf.close()
  elif args.out_format == 'png':
    fig.savefig(args.out + '.png', format='png', dpi=args.dpi)
  elif args.out_format == 'all':
    fig.savefig(pdf, format='pdf')
    pdf.close()
    fig.savefig(args.out + '.png', format='png', dpi=args.dpi)
    fig.savefig(args.out + '.eps', format='eps')
  elif args.out_format == 'eps':
    fig.savefig(args.out + '.eps', format='eps')


def PlotTwoDimension(data_list, ax, args):
  """Plot two dimensional data.

  Args:
    data_list: a list of Data objects
    ax: a matplotlib axis object
    args: an argparse arguments object
  """
  if args.mode == 'scatter':
    marker = 'o'
    args.linewidth = 0.0
    alpha = args.alpha
  else:
    marker = None
    alpha = 1.0
  for i, data in enumerate(data_list, 0):
    if len(data.data[0]) != len(data.data[1]):
      data.data[0] = range(1, len(data.data[1]) + 1)
      args.xmin = 1
      args.xmax=len(data.data[1])
    ax.add_line(
      lines.Line2D(xdata=data.data[0],
                   ydata=data.data[1],
                   color=args.color_list[i],
                   marker=marker,
                   markersize=args.markersize,
                   markerfacecolor=args.color_list[i],
                   markeredgecolor='None',
                   alpha=alpha,
                   linewidth=args.linewidth))
    if args.regression:
      data.data[0].sort()
      data.data[1].sort()
      rxlist = numpy.array(data.data[0])
      rylist = numpy.array(data.data[1])
      A = numpy.array([xlist, numpy.ones(len(rxlist))])
      try:
        w = numpy.linalg.lstsq(A.T, rylist)[0]
      except ValueError:
        sys.stderr.write('Warning, unable to perform regression!')
        return
      fitline = (w[0] * rxlist) + w[1]
      ax.add_line(
        lines.Line2D(xdata=[rxlist[0], rxlist[-1]],
                     ydata=[fitline[0], fitline[-1]],
                     color='red',
                     linestyle='--'))
      slope, intercept, r, p, stderr = linregress(rxlist, rylist)
      ax.text(x=(rxlist[0] + rxlist[-1]) / 2.0,
              y=(fitline[0] + fitline[-1]) / 2.0,
              s=r'%f * x + %f, $r^2$=%f' % (w[0], w[1], r * r))


def ReadFiles(args):
  """Read and parse all input files.

  Args:
    args: an argparse arguments object

  Returns:
    data_list: a list of numpy arrays of dimension (n by 1) or (n by 2)
      depending on the --mode flag.
  """
  data_list = []
  for a_file in args.files:
    f = open(a_file, 'r')
    xlist = []
    ylist = []
    i = 0
    for line in f:
      line = line.strip()
      if line.startswith('#'):
        continue
      d = line.split()
      if d[0].lower() == 'nan':
        continue
      if len(d) > 1:
        # two values
        if d[0].lower() == 'nan' or d[1].lower() == 'nan':
          continue
        xlist.append(float(d[0]))
        ylist.append(float(d[1]))
        if args.xmin > xlist[-1]:
          args.xmin = xlist[-1]
        if args.xmax < xlist[-1]:
          args.xmax = xlist[-1]
      else:
        ylist.append(float(d[0]))
        if args.mode in ('scatter', 'line'):
          xlist.append(i)
          i += 1
      if args.ymin > ylist[-1]:
        args.ymin = ylist[-1]
      if args.ymax < ylist[-1]:
        args.ymax = ylist[-1]
    f.close()
    d = Data()
    d.label = os.path.basename(a_file)
    if ylist == []:
      d.data = numpy.array(xlist)
    else:
      d.data = numpy.array([xlist, ylist])
    data_list.append(d)
  return data_list


def PlotOneDimension(data_list, ax, args):
  """Plot one dimensional data.

  Args:
    data_list: a list of Data objects.
    ax: a matplotlib axis object.
    args: an argparse arguments object.
  """
  if args.mode == 'bar' or args.mode == 'column':
    PlotColumns(data_list, ax, args)
  elif args.mode == 'hist':
    PlotHistogram(data_list, ax, args)
  elif args.mode == 'tick':
    PlotTicks(data_list, ax, args)
  elif args.mode == 'point':
    PlotPoints(data_list, ax, args)


def PlotHistogram(data_list, ax, args):
  """Plot one dimensional data as histogram.

  Args:
    data_list: a list of Data objects.
    ax: a matplotlib axis object.
    args: an argparse arguments object.
  """
  width = 2.0 / 3.0 / len(data_list)
  datas = []
  for data in data_list:
    datas.append(data.data[1])
  n, bins, patch_groups = ax.hist(
    datas, color=args.color_list[0:len(datas)], histtype='bar')
  for pg in patch_groups:
    if isinstance(pg, matplotlib.container.BarContainer):
      # if there are multiple files, pg will be a BarContainer
      for patch in pg.patches:
        patch.set_edgecolor('none')
    else:
      # if there is only one file to plot, pg is a Rectangle
      pg.set_edgecolor('white')


def PlotColumns(data_list, ax, args):
  """Plot one dimensional data as column / bar plot.

  Args:
    data_list: a list of Data objects.
    ax: a matplotlib axis object.
    args: an argparse arguments object.
  """
  width = 2.0 / 3.0 / len(data_list)
  for i, data in enumerate(data_list, 0):
    data.data[0] = range(0, len(data.data[1]))
    data.data[0] = numpy.add(data.data[0], width * i)  # offset
    args.xmin = 0
    args.xmax = len(data.data[1])
    rects = ax.bar(data.data[0],
                   data.data[1],
                   width,
                   color=args.color_list[i],
                   linewidth=0.0,
                   alpha=1.0)
    ax.xaxis.set_ticklabels([])


def GetTickYValues(i, args):
  """Produce the lower and upper y values for a Tick plot.

  Args:
    i: Integer offset of this set of values.
    args: an argparse arguments object.

  Returns:
    y0, y1: the lower and upper y values for a Tick Plot
  """
  if args.jitter:
    lo = numpy.random.uniform(low=0.0, high=0.3)
    return i + lo, i + lo + 0.1
  else:
    return i, i + 0.8


def PlotTicks(data_list, ax, args):
  """Plot one dimensional data as tick marks on a line.

  Args:
    data_list: a list of Data objects
    ax: a matplotlib axis object
    args: an argparse arguments object
  """
  data_min = min(map(numpy.min, map(lambda x: x.data[1], data_list)))
  data_max = max(map(numpy.max, map(lambda x: x.data[1], data_list)))
  data_range = data_max - data_min
  if data_range == 0.0:
    data_min, data_max, data_range = -0.5, 0.5, 1.0
  data_min -= data_range * 0.1
  data_max += data_range * 0.1
  for i, data in enumerate(data_list, 0):
    for d in data.data[1]:
      y0, y1 = GetTickYValues(i, args)
      ax.add_line(
        lines.Line2D(xdata=[d, d],
                     ydata=[y0, y1],
                     color=args.color_list[i],
                     marker=None,
                     markersize=args.markersize,
                     markerfacecolor=args.color_list[i],
                     markeredgecolor='None',
                     alpha=args.alpha,
                     linewidth=args.linewidth))
  ax.set_ylim([0.0, len(data_list)])
  ax.set_xlim([data_min, data_max])
  ax.yaxis.set_ticks_position('none')
  ax.yaxis.set_ticks([])


def GetPointYValues(n, i, args):
  """Produce the y values for a Point plot.

  Args:
    n: number of values to produce.
    i: Integer offset of this set of values.
    args: an argparse arguments object.

  Returns:
    y: a list of y values for a Point plot.
  """
  if args.jitter:
    return numpy.random.uniform(low=i, high=i + 0.5, size=n)
  else:
    return [i] * n


def PlotPoints(data_list, ax, args):
  """Plot one dimensional data as points on a line.

  Args:
    data_list: a list of Data objects
    ax: a matplotlib axis object
    args: an argparse arguments object
  """
  data_min = min(map(numpy.min, map(lambda x: x.data[1], data_list)))
  data_max = max(map(numpy.max, map(lambda x: x.data[1], data_list)))
  data_range = data_max - data_min
  if data_range == 0.0:
    data_min, data_max, data_range = -0.5, 0.5, 1.0
  data_min -= data_range * 0.1
  data_max += data_range * 0.1
  for i, data in enumerate(data_list, 0):
    data.data[0] = GetPointYValues(len(data.data[1]), i, args)
    ax.add_line(
      lines.Line2D(xdata=data.data[1],
                   ydata=data.data[0],
                   color=args.color_list[i],
                   marker='o',
                   markersize=args.markersize,
                   markerfacecolor=args.color_list[i],
                   markeredgecolor='None',
                   alpha=args.alpha,
                   linewidth=0.0))
  ax.set_ylim([-0.5, len(data_list)])
  ax.set_xlim([data_min, data_max])
  ax.yaxis.set_ticks_position('none')
  ax.yaxis.set_ticks([])


def PlotData(data_list, ax, args):
  """Plot all of the data according to input arguments.

  Args:
    data_list: a list of Data objects.
    ax: a matplotlib axis object.
    args: an argparse argument object.
  """
  if args.mode in ('scatter', 'line'):
    PlotTwoDimension(data_list, ax, args)
  elif args.mode in ('bar', 'column', 'hist', 'tick', 'point'):
    PlotOneDimension(data_list, ax, args)


def MakeProxyPlots(args):
  """Make some proxy plots for use with legends.

  Proxy plots are plots that are not actually drawn but whose
  colors are used for correctly populating a legend.

  Args:
    args: an argparse argument object.

  Returns:
    proxy_plots: A list of matplotlib plot objects.
  """
  proxy_plots = []
  for i, afile in enumerate(args.files, 0):
    proxy_plots.append(
      plt.Rectangle(
        (0, 0), 1, 1, fc=args.color_list[i % len(args.color_list)],
        ec=args.color_list[i % len(args.color_list)]))
  return proxy_plots


def MakeLegendLabels(args):
  """Make labels for use with legends.

  Args:
    args: an argparse argument object

  Returns:
    legend_labels: A list of strings.
  """
  legend_labels = []
  for afile in args.files:
    legend_labels.append(os.path.basename(afile))
  return legend_labels


def CleanAxis(ax, args):
  """Clean the axis up, apply scales, add legend.

  Args:
    ax: a matplotlib axis object
    args: an argparse argument object
  """
  # logarithmic axes
  if args.is_log_y:
    ax.set_yscale('log')
  else:
    arange = args.ymax - args.ymin
    if args.mode not in ('hist', 'tick', 'point'):
      ax.set_ylim([args.ymin - arange * 0.05, args.ymax + arange * 0.05])
  if args.is_log_x:
    ax.set_xscale('log')
  else:
    arange = args.xmax - args.xmin
    if args.mode not in ('hist', 'tick', 'point'):
      ax.set_xlim([args.xmin - arange * 0.05, args.xmax + arange * 0.05])
  # labels
  if args.xlabel != 'sentinel_value':
    ax.set_xlabel(args.xlabel)
  if args.ylabel != 'sentinel_value':
    if args.mode in ('tick', 'point'):
      sys.stderr.write('Warning, --ylabel specified while '
                       '--mode=tick, ylabel not displayed\n')
    else:
      ax.set_ylabel(args.ylabel)
  if args.title != 'sentinel_value':
    ax.set_title(args.title)
  # legend
  if args.is_legend:
    proxy_plots = MakeProxyPlots(args)
    legend_labels = MakeLegendLabels(args)
    leg = plt.legend(proxy_plots, legend_labels, 'upper right', numpoints=1)
    leg._drawFrame = False


def main():
  usage = '%(prog)s file1 file2 file3... [options]\n\n'
  description = ('%(prog)s is a tool to produce quick plots. col1 '
                 'of input file is x value col2 is y value. If '
                 'the --mode is column/bar/hist then only col1 is '
                 'used.')
  parser = ArgumentParser(usage=usage, description=description)
  InitArguments(parser)
  args = parser.parse_args()
  CheckArguments(args, parser)
  fig, pdf = InitImage(args)
  ax = EstablishAxes(fig, args)

  data_list = ReadFiles(args)
  PlotData(data_list, ax, args)

  CleanAxis(ax, args)
  WriteImage(fig, pdf, args)


if __name__ == '__main__':
    main()
