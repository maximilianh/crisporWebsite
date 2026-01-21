These examples illustrate how a PS file can be made much smaller by replacing multiple drawing commands with postscript procedures (user-defined functions).

File                     Size   %Size
------------------------ -----  -----
ps-0-default.ps          25627   100%
ps-1-optimize-colors.ps  17799  69.5%
ps-2-optimize-procs.ps    5185  20.2%
ps-3-optimize-nuclist.ps  3323  13.0%

ps-0-default.ps - The full-size file with drawing instructions.

ps-1-optimize-colors.ps
  Colors replaced with defined procedures 
  (e.g. "0.00 0.00 0.00 setrgbcolor" replaced with "blk")

ps-2-optimize-procs.ps
  Larger drawing routines replaced with procedures (e.g. drawNuc drawBkbn etc)

ps-3-optimize-nuclist.ps
  Replaced drawing routines with explicit listing of nucleotides and "helices" 
  (which are expanded to pairs) and number labels.
  Procedures loop through lists, drawing nucs and the backbone.
  This allows nucleotide positions to be listed
  once instead of once-per-draw-command.
  