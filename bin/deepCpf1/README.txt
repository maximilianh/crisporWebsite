This is code from the paper  https://www.nature.com/articles/nbt.4061

Deep learning improves prediction of CRISPR–Cpf1 guide RNA activity

Hui Kwon Kim, Seonwoo Min, Myungjae Song, Soobin Jung, Jae Woo Choi, Younggwang Kim, Sangeun Lee, Sungroh Yoon & Hyongbum (Henry) Kim

Source code:
https://github.com/MyungjaeSong/Paired-Library

The code on github does not seem to have a license.

Website: http://deepcrispr.info/

BUT:
The code here is actually not from GitHub but was sent by authors by email,
as I was unable to reproduce the scores from output_example.txt.
The author of the code sent these very quickly by email, mswzeus@gmail.com

The updated code here contains different models, the *100 factor for CA array
and I added a line to force the theano backend (the scores do not work with
the default tensorflow backend).

I also modified DeepCpf1.py to make it possible to call it from my code.
