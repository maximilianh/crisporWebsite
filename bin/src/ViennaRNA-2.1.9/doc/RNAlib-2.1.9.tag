<?xml version='1.0' encoding='UTF-8' standalone='yes' ?>
<tagfile>
  <compound kind="file">
    <name>2Dfold.h</name>
    <path>/home/mescalin/ronny/WORK/ViennaRNA/H/</path>
    <filename>2Dfold_8h</filename>
    <member kind="function">
      <type>TwoDfold_vars *</type>
      <name>get_TwoDfold_variables</name>
      <anchorfile>group__kl__neighborhood__mfe.html</anchorfile>
      <anchor>gac9284f132cf0eaa0a2f43590eda05488</anchor>
      <arglist>(const char *seq, const char *structure1, const char *structure2, int circ)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>destroy_TwoDfold_variables</name>
      <anchorfile>group__kl__neighborhood__mfe.html</anchorfile>
      <anchor>ga05bf4f31d216b1b160fd2d3d68e9b487</anchor>
      <arglist>(TwoDfold_vars *our_variables)</arglist>
    </member>
    <member kind="function">
      <type>TwoDfold_solution *</type>
      <name>TwoDfoldList</name>
      <anchorfile>group__kl__neighborhood__mfe.html</anchorfile>
      <anchor>ga47da790166020558d27323aef489703e</anchor>
      <arglist>(TwoDfold_vars *vars, int distance1, int distance2)</arglist>
    </member>
    <member kind="function">
      <type>char *</type>
      <name>TwoDfold_backtrack_f5</name>
      <anchorfile>group__kl__neighborhood__mfe.html</anchorfile>
      <anchor>gaf4dc05bf8fc1ea53acd7aeb798ba80c2</anchor>
      <arglist>(unsigned int j, int k, int l, TwoDfold_vars *vars)</arglist>
    </member>
  </compound>
  <compound kind="file">
    <name>2Dpfold.h</name>
    <path>/home/mescalin/ronny/WORK/ViennaRNA/H/</path>
    <filename>2Dpfold_8h</filename>
    <member kind="function">
      <type>TwoDpfold_vars *</type>
      <name>get_TwoDpfold_variables</name>
      <anchorfile>group__kl__neighborhood__pf.html</anchorfile>
      <anchor>ga1aca740e2a75ab2b2951538266e53d64</anchor>
      <arglist>(const char *seq, const char *structure1, char *structure2, int circ)</arglist>
    </member>
    <member kind="function">
      <type>TwoDpfold_vars *</type>
      <name>get_TwoDpfold_variables_from_MFE</name>
      <anchorfile>group__kl__neighborhood__pf.html</anchorfile>
      <anchor>gacc2f66da7ee62096cab629fce7112216</anchor>
      <arglist>(TwoDfold_vars *mfe_vars)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>destroy_TwoDpfold_variables</name>
      <anchorfile>group__kl__neighborhood__pf.html</anchorfile>
      <anchor>gafe994291458ee2ac34d3eb825ef62a15</anchor>
      <arglist>(TwoDpfold_vars *vars)</arglist>
    </member>
    <member kind="function">
      <type>TwoDpfold_solution *</type>
      <name>TwoDpfoldList</name>
      <anchorfile>group__kl__neighborhood__pf.html</anchorfile>
      <anchor>ga3e1cd3b24eb635c65181182cbb4ae3eb</anchor>
      <arglist>(TwoDpfold_vars *vars, int maxDistance1, int maxDistance2)</arglist>
    </member>
    <member kind="function">
      <type>char *</type>
      <name>TwoDpfold_pbacktrack</name>
      <anchorfile>group__kl__neighborhood__stochbt.html</anchorfile>
      <anchor>gae251288f50dd4ae7d315af0085775f71</anchor>
      <arglist>(TwoDpfold_vars *vars, int d1, int d2)</arglist>
    </member>
    <member kind="function">
      <type>char *</type>
      <name>TwoDpfold_pbacktrack5</name>
      <anchorfile>group__kl__neighborhood__stochbt.html</anchorfile>
      <anchor>ga13430ac6a7f90df426774f131647d2c7</anchor>
      <arglist>(TwoDpfold_vars *vars, int d1, int d2, unsigned int length)</arglist>
    </member>
  </compound>
  <compound kind="file">
    <name>alifold.h</name>
    <path>/home/mescalin/ronny/WORK/ViennaRNA/H/</path>
    <filename>alifold_8h</filename>
    <member kind="function">
      <type>void</type>
      <name>update_alifold_params</name>
      <anchorfile>alifold_8h.html</anchorfile>
      <anchor>ac484c6bd429bafbd353b91044508d8e9</anchor>
      <arglist>(void)</arglist>
    </member>
    <member kind="function">
      <type>float</type>
      <name>alifold</name>
      <anchorfile>group__consensus__mfe__fold.html</anchorfile>
      <anchor>ga4cf00f0659e5f0480335d69e797f05b1</anchor>
      <arglist>(const char **strings, char *structure)</arglist>
    </member>
    <member kind="function">
      <type>float</type>
      <name>circalifold</name>
      <anchorfile>group__consensus__mfe__fold.html</anchorfile>
      <anchor>gadbd3b0b1c144cbfb4efe704b2b260f96</anchor>
      <arglist>(const char **strings, char *structure)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>free_alifold_arrays</name>
      <anchorfile>group__consensus__mfe__fold.html</anchorfile>
      <anchor>ga72095e4554b5d577250ea14c42acc49e</anchor>
      <arglist>(void)</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>get_mpi</name>
      <anchorfile>group__consensus__fold.html</anchorfile>
      <anchor>gaa2d600be90844094ec145ea14a314d2f</anchor>
      <arglist>(char *Alseq[], int n_seq, int length, int *mini)</arglist>
    </member>
    <member kind="function">
      <type>float **</type>
      <name>readribosum</name>
      <anchorfile>group__consensus__fold.html</anchorfile>
      <anchor>ga5e125c9586fcd4e2e1559fe76f7289cc</anchor>
      <arglist>(char *name)</arglist>
    </member>
    <member kind="function">
      <type>float</type>
      <name>energy_of_alistruct</name>
      <anchorfile>group__consensus__fold.html</anchorfile>
      <anchor>ga1c48869c03b49a342bf4cbdd61900081</anchor>
      <arglist>(const char **sequences, const char *structure, int n_seq, float *energy)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>encode_ali_sequence</name>
      <anchorfile>group__consensus__fold.html</anchorfile>
      <anchor>gaa3e40277c837d6f7603afe319884c786</anchor>
      <arglist>(const char *sequence, short *S, short *s5, short *s3, char *ss, unsigned short *as, int circ)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>alloc_sequence_arrays</name>
      <anchorfile>group__consensus__fold.html</anchorfile>
      <anchor>ga8a560930f7f2582cc3967723a86cfdfa</anchor>
      <arglist>(const char **sequences, short ***S, short ***S5, short ***S3, unsigned short ***a2s, char ***Ss, int circ)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>free_sequence_arrays</name>
      <anchorfile>group__consensus__fold.html</anchorfile>
      <anchor>ga298a420a8c879202e2617b3f724fde38</anchor>
      <arglist>(unsigned int n_seq, short ***S, short ***S5, short ***S3, unsigned short ***a2s, char ***Ss)</arglist>
    </member>
    <member kind="function">
      <type>float</type>
      <name>alipf_fold_par</name>
      <anchorfile>group__consensus__pf__fold.html</anchorfile>
      <anchor>ga4d2ff54d8210fc7cceeeff389d4dbd1d</anchor>
      <arglist>(const char **sequences, char *structure, plist **pl, pf_paramT *parameters, int calculate_bppm, int is_constrained, int is_circular)</arglist>
    </member>
    <member kind="function">
      <type>float</type>
      <name>alipf_fold</name>
      <anchorfile>group__consensus__pf__fold.html</anchorfile>
      <anchor>gad32ded7d753ccaf211ab35782d1f42a9</anchor>
      <arglist>(const char **sequences, char *structure, plist **pl)</arglist>
    </member>
    <member kind="function">
      <type>float</type>
      <name>alipf_circ_fold</name>
      <anchorfile>group__consensus__pf__fold.html</anchorfile>
      <anchor>ga6b4dde1d43b79ab3753508c46cf50363</anchor>
      <arglist>(const char **sequences, char *structure, plist **pl)</arglist>
    </member>
    <member kind="function">
      <type>double *</type>
      <name>export_ali_bppm</name>
      <anchorfile>group__consensus__pf__fold.html</anchorfile>
      <anchor>gadaaf83394216413505e48d913dbc1b4e</anchor>
      <arglist>(void)</arglist>
    </member>
    <member kind="function">
      <type>char *</type>
      <name>alipbacktrack</name>
      <anchorfile>group__consensus__stochbt.html</anchorfile>
      <anchor>ga0df40248788f0fb17ebdc59d74116d1c</anchor>
      <arglist>(double *prob)</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>get_alipf_arrays</name>
      <anchorfile>group__consensus__fold.html</anchorfile>
      <anchor>ga0cc49457fd79eeb04d4a7f97c868b09b</anchor>
      <arglist>(short ***S_p, short ***S5_p, short ***S3_p, unsigned short ***a2s_p, char ***Ss_p, double **qb_p, double **qm_p, double **q1k_p, double **qln_p, short **pscore)</arglist>
    </member>
    <member kind="variable">
      <type>double</type>
      <name>cv_fact</name>
      <anchorfile>group__consensus__fold.html</anchorfile>
      <anchor>gaf3cbac6ff5d706d6e414677841ddf94c</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>double</type>
      <name>nc_fact</name>
      <anchorfile>group__consensus__fold.html</anchorfile>
      <anchor>ga502948a122a2af5b914355b1f3ea2f61</anchor>
      <arglist></arglist>
    </member>
  </compound>
  <compound kind="file">
    <name>cofold.h</name>
    <path>/home/mescalin/ronny/WORK/ViennaRNA/H/</path>
    <filename>cofold_8h</filename>
    <member kind="function">
      <type>float</type>
      <name>cofold</name>
      <anchorfile>group__mfe__cofold.html</anchorfile>
      <anchor>gabc8517f22cfe70595ee81fc837910d52</anchor>
      <arglist>(const char *sequence, char *structure)</arglist>
    </member>
    <member kind="function">
      <type>float</type>
      <name>cofold_par</name>
      <anchorfile>group__mfe__cofold.html</anchorfile>
      <anchor>gafe430060533f14b11fc611f60b3f1f6f</anchor>
      <arglist>(const char *string, char *structure, paramT *parameters, int is_constrained)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>free_co_arrays</name>
      <anchorfile>group__mfe__cofold.html</anchorfile>
      <anchor>gaafb33d7473eb9af9d1b168ca8761c41a</anchor>
      <arglist>(void)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>update_cofold_params</name>
      <anchorfile>group__mfe__cofold.html</anchorfile>
      <anchor>ga4fcbf34e77b99bfbb2333d2ab0c41a57</anchor>
      <arglist>(void)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>export_cofold_arrays_gq</name>
      <anchorfile>group__mfe__cofold.html</anchorfile>
      <anchor>ga5f5bf4df35d0554f6ace9579f8744c48</anchor>
      <arglist>(int **f5_p, int **c_p, int **fML_p, int **fM1_p, int **fc_p, int **ggg_p, int **indx_p, char **ptype_p)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>export_cofold_arrays</name>
      <anchorfile>group__mfe__cofold.html</anchorfile>
      <anchor>ga5cb6b59983f1f74ccc00b9b9c4e84482</anchor>
      <arglist>(int **f5_p, int **c_p, int **fML_p, int **fM1_p, int **fc_p, int **indx_p, char **ptype_p)</arglist>
    </member>
    <member kind="function">
      <type>SOLUTION *</type>
      <name>zukersubopt</name>
      <anchorfile>group__subopt__zuker.html</anchorfile>
      <anchor>ga0d5104e3ecf119d8eabd40aa5fe47f90</anchor>
      <arglist>(const char *string)</arglist>
    </member>
    <member kind="function">
      <type>SOLUTION *</type>
      <name>zukersubopt_par</name>
      <anchorfile>group__subopt__zuker.html</anchorfile>
      <anchor>ga6d98a9450d1affadf144ac79f543da8c</anchor>
      <arglist>(const char *string, paramT *parameters)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>get_monomere_mfes</name>
      <anchorfile>cofold_8h.html</anchorfile>
      <anchor>a4958b517c613e4d2afd5bce6c1060a79</anchor>
      <arglist>(float *e1, float *e2)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>initialize_cofold</name>
      <anchorfile>cofold_8h.html</anchorfile>
      <anchor>afee0c32208aa2ac97338b6e3fbad7fa5</anchor>
      <arglist>(int length)</arglist>
    </member>
  </compound>
  <compound kind="file">
    <name>convert_epars.h</name>
    <path>/home/mescalin/ronny/WORK/ViennaRNA/H/</path>
    <filename>convert__epars_8h</filename>
    <member kind="define">
      <type>#define</type>
      <name>VRNA_CONVERT_OUTPUT_ALL</name>
      <anchorfile>group__energy__parameters__convert.html</anchorfile>
      <anchor>ga8dc6aee5a806c49b71557152f9616bc4</anchor>
      <arglist></arglist>
    </member>
    <member kind="define">
      <type>#define</type>
      <name>VRNA_CONVERT_OUTPUT_HP</name>
      <anchorfile>group__energy__parameters__convert.html</anchorfile>
      <anchor>gaf66fe2cb11dfcfd32d791049c254a8a4</anchor>
      <arglist></arglist>
    </member>
    <member kind="define">
      <type>#define</type>
      <name>VRNA_CONVERT_OUTPUT_STACK</name>
      <anchorfile>group__energy__parameters__convert.html</anchorfile>
      <anchor>gad23522d63f8d4c50d5a5deee9bee3ef2</anchor>
      <arglist></arglist>
    </member>
    <member kind="define">
      <type>#define</type>
      <name>VRNA_CONVERT_OUTPUT_MM_HP</name>
      <anchorfile>group__energy__parameters__convert.html</anchorfile>
      <anchor>gaa892c7b4957459090f3e08da298cc347</anchor>
      <arglist></arglist>
    </member>
    <member kind="define">
      <type>#define</type>
      <name>VRNA_CONVERT_OUTPUT_MM_INT</name>
      <anchorfile>group__energy__parameters__convert.html</anchorfile>
      <anchor>ga4ff223fb1f9c62cd92d9ab811ad03d55</anchor>
      <arglist></arglist>
    </member>
    <member kind="define">
      <type>#define</type>
      <name>VRNA_CONVERT_OUTPUT_MM_INT_1N</name>
      <anchorfile>group__energy__parameters__convert.html</anchorfile>
      <anchor>gaf5d3743219f83c6348155cd81e755bbb</anchor>
      <arglist></arglist>
    </member>
    <member kind="define">
      <type>#define</type>
      <name>VRNA_CONVERT_OUTPUT_MM_INT_23</name>
      <anchorfile>group__energy__parameters__convert.html</anchorfile>
      <anchor>ga78382ec622ba99e0ac2262317bdd7316</anchor>
      <arglist></arglist>
    </member>
    <member kind="define">
      <type>#define</type>
      <name>VRNA_CONVERT_OUTPUT_MM_MULTI</name>
      <anchorfile>group__energy__parameters__convert.html</anchorfile>
      <anchor>gae67af9f1cdf7baf2865481282a5d1034</anchor>
      <arglist></arglist>
    </member>
    <member kind="define">
      <type>#define</type>
      <name>VRNA_CONVERT_OUTPUT_MM_EXT</name>
      <anchorfile>group__energy__parameters__convert.html</anchorfile>
      <anchor>gaf14ead7ef1fdbe725ade653750fc51e3</anchor>
      <arglist></arglist>
    </member>
    <member kind="define">
      <type>#define</type>
      <name>VRNA_CONVERT_OUTPUT_DANGLE5</name>
      <anchorfile>group__energy__parameters__convert.html</anchorfile>
      <anchor>ga036ffd996d8c8a9acf631760dd1da24b</anchor>
      <arglist></arglist>
    </member>
    <member kind="define">
      <type>#define</type>
      <name>VRNA_CONVERT_OUTPUT_DANGLE3</name>
      <anchorfile>group__energy__parameters__convert.html</anchorfile>
      <anchor>ga34a8a5479ef885834ef32f3fb43d79bc</anchor>
      <arglist></arglist>
    </member>
    <member kind="define">
      <type>#define</type>
      <name>VRNA_CONVERT_OUTPUT_INT_11</name>
      <anchorfile>group__energy__parameters__convert.html</anchorfile>
      <anchor>ga079aafefd5f8ab57ee5120099a34bd25</anchor>
      <arglist></arglist>
    </member>
    <member kind="define">
      <type>#define</type>
      <name>VRNA_CONVERT_OUTPUT_INT_21</name>
      <anchorfile>group__energy__parameters__convert.html</anchorfile>
      <anchor>gacf770881d9034431ebe741642342a1f9</anchor>
      <arglist></arglist>
    </member>
    <member kind="define">
      <type>#define</type>
      <name>VRNA_CONVERT_OUTPUT_INT_22</name>
      <anchorfile>group__energy__parameters__convert.html</anchorfile>
      <anchor>gaa307671e2631cdacad9cbe4c6583b05f</anchor>
      <arglist></arglist>
    </member>
    <member kind="define">
      <type>#define</type>
      <name>VRNA_CONVERT_OUTPUT_BULGE</name>
      <anchorfile>group__energy__parameters__convert.html</anchorfile>
      <anchor>ga7092fe0be4de6f02cc0bf08e81af726a</anchor>
      <arglist></arglist>
    </member>
    <member kind="define">
      <type>#define</type>
      <name>VRNA_CONVERT_OUTPUT_INT</name>
      <anchorfile>group__energy__parameters__convert.html</anchorfile>
      <anchor>gac5c2289fdf8ff1b980976d1613ff943a</anchor>
      <arglist></arglist>
    </member>
    <member kind="define">
      <type>#define</type>
      <name>VRNA_CONVERT_OUTPUT_ML</name>
      <anchorfile>group__energy__parameters__convert.html</anchorfile>
      <anchor>gaf2c8755d64eff3852aa45df9ac80a4fe</anchor>
      <arglist></arglist>
    </member>
    <member kind="define">
      <type>#define</type>
      <name>VRNA_CONVERT_OUTPUT_MISC</name>
      <anchorfile>group__energy__parameters__convert.html</anchorfile>
      <anchor>ga46d5b1535ae86060b6317565b7c6b40b</anchor>
      <arglist></arglist>
    </member>
    <member kind="define">
      <type>#define</type>
      <name>VRNA_CONVERT_OUTPUT_SPECIAL_HP</name>
      <anchorfile>group__energy__parameters__convert.html</anchorfile>
      <anchor>gaa1ff48a79642d69579d1766561ec6db6</anchor>
      <arglist></arglist>
    </member>
    <member kind="define">
      <type>#define</type>
      <name>VRNA_CONVERT_OUTPUT_VANILLA</name>
      <anchorfile>group__energy__parameters__convert.html</anchorfile>
      <anchor>ga0d4e8a836bb4864ab5129c085dbf592d</anchor>
      <arglist></arglist>
    </member>
    <member kind="define">
      <type>#define</type>
      <name>VRNA_CONVERT_OUTPUT_NINIO</name>
      <anchorfile>group__energy__parameters__convert.html</anchorfile>
      <anchor>ga2eb0462f16939ddacdaf751a88d675ce</anchor>
      <arglist></arglist>
    </member>
    <member kind="define">
      <type>#define</type>
      <name>VRNA_CONVERT_OUTPUT_DUMP</name>
      <anchorfile>group__energy__parameters__convert.html</anchorfile>
      <anchor>gac86976e9c2a55b3a6481ea60044f6098</anchor>
      <arglist></arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>convert_parameter_file</name>
      <anchorfile>group__energy__parameters__convert.html</anchorfile>
      <anchor>gafbe538bc4eb2cf2a33326e1010005f8a</anchor>
      <arglist>(const char *iname, const char *oname, unsigned int options)</arglist>
    </member>
  </compound>
  <compound kind="file">
    <name>data_structures.h</name>
    <path>/home/mescalin/ronny/WORK/ViennaRNA/H/</path>
    <filename>data__structures_8h</filename>
    <class kind="struct">plist</class>
    <class kind="struct">cpair</class>
    <class kind="struct">COORDINATE</class>
    <class kind="struct">sect</class>
    <class kind="struct">bondT</class>
    <class kind="struct">bondTEn</class>
    <class kind="struct">model_detailsT</class>
    <class kind="struct">paramT</class>
    <class kind="struct">pf_paramT</class>
    <class kind="struct">PAIR</class>
    <class kind="struct">INTERVAL</class>
    <class kind="struct">SOLUTION</class>
    <class kind="struct">cofoldF</class>
    <class kind="struct">ConcEnt</class>
    <class kind="struct">pairpro</class>
    <class kind="struct">pair_info</class>
    <class kind="struct">move_t</class>
    <class kind="struct">intermediate_t</class>
    <class kind="struct">path_t</class>
    <class kind="struct">pu_contrib</class>
    <class kind="struct">interact</class>
    <class kind="struct">pu_out</class>
    <class kind="struct">constrain</class>
    <class kind="struct">duplexT</class>
    <class kind="struct">folden</class>
    <class kind="struct">snoopT</class>
    <class kind="struct">dupVar</class>
    <class kind="struct">TwoDfold_solution</class>
    <class kind="struct">TwoDfold_vars</class>
    <class kind="struct">TwoDpfold_solution</class>
    <class kind="struct">TwoDpfold_vars</class>
    <member kind="define">
      <type>#define</type>
      <name>MAXALPHA</name>
      <anchorfile>data__structures_8h.html</anchorfile>
      <anchor>a05a5ffe718aa431d97419a12fb082379</anchor>
      <arglist></arglist>
    </member>
    <member kind="define">
      <type>#define</type>
      <name>MAXDOS</name>
      <anchorfile>data__structures_8h.html</anchorfile>
      <anchor>a5ec740b80afb4906ba4311dbd8ddbd89</anchor>
      <arglist></arglist>
    </member>
  </compound>
  <compound kind="file">
    <name>dist_vars.h</name>
    <path>/home/mescalin/ronny/WORK/ViennaRNA/H/</path>
    <filename>dist__vars_8h</filename>
    <class kind="struct">Postorder_list</class>
    <class kind="struct">Tree</class>
    <class kind="struct">swString</class>
    <member kind="variable">
      <type>int</type>
      <name>edit_backtrack</name>
      <anchorfile>dist__vars_8h.html</anchorfile>
      <anchor>aa03194c513af6b860e7b33e370b82bdb</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>char *</type>
      <name>aligned_line</name>
      <anchorfile>dist__vars_8h.html</anchorfile>
      <anchor>ac1605fe3448ad0a0b809c4fb8f6a854a</anchor>
      <arglist>[4]</arglist>
    </member>
    <member kind="variable">
      <type>int</type>
      <name>cost_matrix</name>
      <anchorfile>dist__vars_8h.html</anchorfile>
      <anchor>ab65d8ff14c6937612212526a60f59b3c</anchor>
      <arglist></arglist>
    </member>
  </compound>
  <compound kind="file">
    <name>duplex.h</name>
    <path>/home/mescalin/ronny/WORK/ViennaRNA/H/</path>
    <filename>duplex_8h</filename>
  </compound>
  <compound kind="file">
    <name>edit_cost.h</name>
    <path>/home/mescalin/ronny/WORK/ViennaRNA/H/</path>
    <filename>edit__cost_8h</filename>
  </compound>
  <compound kind="file">
    <name>energy_const.h</name>
    <path>/home/mescalin/ronny/WORK/ViennaRNA/H/</path>
    <filename>energy__const_8h</filename>
    <member kind="define">
      <type>#define</type>
      <name>GASCONST</name>
      <anchorfile>energy__const_8h.html</anchorfile>
      <anchor>ab1e4a8d82f24ed5db01dde5f25269cf1</anchor>
      <arglist></arglist>
    </member>
    <member kind="define">
      <type>#define</type>
      <name>K0</name>
      <anchorfile>energy__const_8h.html</anchorfile>
      <anchor>a307c72605e3713972b4f4fb2d53ea20e</anchor>
      <arglist></arglist>
    </member>
    <member kind="define">
      <type>#define</type>
      <name>INF</name>
      <anchorfile>energy__const_8h.html</anchorfile>
      <anchor>a12c2040f25d8e3a7b9e1c2024c618cb6</anchor>
      <arglist></arglist>
    </member>
    <member kind="define">
      <type>#define</type>
      <name>FORBIDDEN</name>
      <anchorfile>energy__const_8h.html</anchorfile>
      <anchor>a5064c29ab2d1e20c2304b3c67562774d</anchor>
      <arglist></arglist>
    </member>
    <member kind="define">
      <type>#define</type>
      <name>BONUS</name>
      <anchorfile>energy__const_8h.html</anchorfile>
      <anchor>a96a9822fa134450197dd454b1478a193</anchor>
      <arglist></arglist>
    </member>
    <member kind="define">
      <type>#define</type>
      <name>NBPAIRS</name>
      <anchorfile>energy__const_8h.html</anchorfile>
      <anchor>a5e75221c779d618eab81e096f37e32ce</anchor>
      <arglist></arglist>
    </member>
    <member kind="define">
      <type>#define</type>
      <name>TURN</name>
      <anchorfile>energy__const_8h.html</anchorfile>
      <anchor>ae646250fd59311356c7e5722a81c3a96</anchor>
      <arglist></arglist>
    </member>
    <member kind="define">
      <type>#define</type>
      <name>MAXLOOP</name>
      <anchorfile>energy__const_8h.html</anchorfile>
      <anchor>ad1bd6eabac419670ddd3c9ed82145988</anchor>
      <arglist></arglist>
    </member>
  </compound>
  <compound kind="file">
    <name>findpath.h</name>
    <path>/home/mescalin/ronny/WORK/ViennaRNA/H/</path>
    <filename>findpath_8h</filename>
    <member kind="function">
      <type>int</type>
      <name>find_saddle</name>
      <anchorfile>findpath_8h.html</anchorfile>
      <anchor>ad0e14268e309af773ecd1fce6244ee50</anchor>
      <arglist>(const char *seq, const char *struc1, const char *struc2, int max)</arglist>
    </member>
    <member kind="function">
      <type>path_t *</type>
      <name>get_path</name>
      <anchorfile>findpath_8h.html</anchorfile>
      <anchor>a0ff35d65c892a3403af937c00a867ef9</anchor>
      <arglist>(const char *seq, const char *s1, const char *s2, int maxkeep)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>free_path</name>
      <anchorfile>findpath_8h.html</anchorfile>
      <anchor>a326e6d1640bbfd035e3869f5f4c188f7</anchor>
      <arglist>(path_t *path)</arglist>
    </member>
  </compound>
  <compound kind="file">
    <name>fold.h</name>
    <path>/home/mescalin/ronny/WORK/ViennaRNA/H/</path>
    <filename>fold_8h</filename>
    <member kind="function">
      <type>float</type>
      <name>fold_par</name>
      <anchorfile>group__mfe__fold.html</anchorfile>
      <anchor>gadb973133c241d57c04b253df35e4d34e</anchor>
      <arglist>(const char *sequence, char *structure, paramT *parameters, int is_constrained, int is_circular)</arglist>
    </member>
    <member kind="function">
      <type>float</type>
      <name>fold</name>
      <anchorfile>group__mfe__fold.html</anchorfile>
      <anchor>gaadafcb0f140795ae62e5ca027e335a9b</anchor>
      <arglist>(const char *sequence, char *structure)</arglist>
    </member>
    <member kind="function">
      <type>float</type>
      <name>circfold</name>
      <anchorfile>group__mfe__fold.html</anchorfile>
      <anchor>ga4ac63ab3e8d9a80ced28b8052d94e423</anchor>
      <arglist>(const char *sequence, char *structure)</arglist>
    </member>
    <member kind="function">
      <type>float</type>
      <name>energy_of_structure</name>
      <anchorfile>group__eval.html</anchorfile>
      <anchor>gaf93986cb3cb29770ec9cca69c9fab8cf</anchor>
      <arglist>(const char *string, const char *structure, int verbosity_level)</arglist>
    </member>
    <member kind="function">
      <type>float</type>
      <name>energy_of_struct_par</name>
      <anchorfile>group__eval.html</anchorfile>
      <anchor>gab5169ea4f72f250e43811463a33f4e40</anchor>
      <arglist>(const char *string, const char *structure, paramT *parameters, int verbosity_level)</arglist>
    </member>
    <member kind="function">
      <type>float</type>
      <name>energy_of_circ_structure</name>
      <anchorfile>group__eval.html</anchorfile>
      <anchor>gaeb14f3664aec67fc03268ac75253f0f8</anchor>
      <arglist>(const char *string, const char *structure, int verbosity_level)</arglist>
    </member>
    <member kind="function">
      <type>float</type>
      <name>energy_of_circ_struct_par</name>
      <anchorfile>group__eval.html</anchorfile>
      <anchor>ga75dc765ee4a1177832bc817c94cf88e5</anchor>
      <arglist>(const char *string, const char *structure, paramT *parameters, int verbosity_level)</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>energy_of_structure_pt</name>
      <anchorfile>group__eval.html</anchorfile>
      <anchor>ga8831445966b761417e713360791299d8</anchor>
      <arglist>(const char *string, short *ptable, short *s, short *s1, int verbosity_level)</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>energy_of_struct_pt_par</name>
      <anchorfile>group__eval.html</anchorfile>
      <anchor>gada4701dd7519b29da75ceac147601f4e</anchor>
      <arglist>(const char *string, short *ptable, short *s, short *s1, paramT *parameters, int verbosity_level)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>free_arrays</name>
      <anchorfile>group__mfe__fold.html</anchorfile>
      <anchor>ga107fdfe5fd641868156bfd849f6866c7</anchor>
      <arglist>(void)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>parenthesis_structure</name>
      <anchorfile>fold_8h.html</anchorfile>
      <anchor>a87b7869bd1d8dc79c60775c74e009e9b</anchor>
      <arglist>(char *structure, bondT *bp, int length)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>parenthesis_zuker</name>
      <anchorfile>fold_8h.html</anchorfile>
      <anchor>a325f3835c68f34fe833b2b7a5828857f</anchor>
      <arglist>(char *structure, bondT *bp, int length)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>update_fold_params</name>
      <anchorfile>group__mfe__fold.html</anchorfile>
      <anchor>ga41bf8f6fa15b94471f7095cad9f0ccf3</anchor>
      <arglist>(void)</arglist>
    </member>
    <member kind="function">
      <type>float</type>
      <name>energy_of_move</name>
      <anchorfile>fold_8h.html</anchorfile>
      <anchor>a539ecaed89730f7644c202f304d7529b</anchor>
      <arglist>(const char *string, const char *structure, int m1, int m2)</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>energy_of_move_pt</name>
      <anchorfile>fold_8h.html</anchorfile>
      <anchor>a49e0ee561be69faf0568213546f6a53f</anchor>
      <arglist>(short *pt, short *s, short *s1, int m1, int m2)</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>loop_energy</name>
      <anchorfile>fold_8h.html</anchorfile>
      <anchor>a507d4fd93f4b398d793ba2402731388d</anchor>
      <arglist>(short *ptable, short *s, short *s1, int i)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>assign_plist_from_db</name>
      <anchorfile>fold_8h.html</anchorfile>
      <anchor>adaa59b81664e2e36cb9932e891558fae</anchor>
      <arglist>(plist **pl, const char *struc, float pr)</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>LoopEnergy</name>
      <anchorfile>fold_8h.html</anchorfile>
      <anchor>a2163034a25c6115d894b199e97e03f6c</anchor>
      <arglist>(int n1, int n2, int type, int type_2, int si1, int sj1, int sp1, int sq1)</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>HairpinE</name>
      <anchorfile>fold_8h.html</anchorfile>
      <anchor>ab327ce11972f5ac069d52c8dedfdb700</anchor>
      <arglist>(int size, int type, int si1, int sj1, const char *string)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>initialize_fold</name>
      <anchorfile>fold_8h.html</anchorfile>
      <anchor>ac3f0a28d9cb609d388b155445073fd20</anchor>
      <arglist>(int length)</arglist>
    </member>
    <member kind="function">
      <type>float</type>
      <name>energy_of_struct</name>
      <anchorfile>fold_8h.html</anchorfile>
      <anchor>ac2b37fea2145c94d925a3f33378ef87b</anchor>
      <arglist>(const char *string, const char *structure)</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>energy_of_struct_pt</name>
      <anchorfile>fold_8h.html</anchorfile>
      <anchor>a27ce6f68512d43bf1fe14a06c9d76d5c</anchor>
      <arglist>(const char *string, short *ptable, short *s, short *s1)</arglist>
    </member>
    <member kind="function">
      <type>float</type>
      <name>energy_of_circ_struct</name>
      <anchorfile>fold_8h.html</anchorfile>
      <anchor>a657222e2758c46bf13b416ef3032e417</anchor>
      <arglist>(const char *string, const char *structure)</arglist>
    </member>
    <member kind="variable">
      <type>int</type>
      <name>logML</name>
      <anchorfile>fold_8h.html</anchorfile>
      <anchor>a80c3c5fd35e7479704cc91d2d0367743</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>int</type>
      <name>uniq_ML</name>
      <anchorfile>fold_8h.html</anchorfile>
      <anchor>a6c5655c8b272e3e6cab74dd0f540294f</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>int</type>
      <name>cut_point</name>
      <anchorfile>fold_8h.html</anchorfile>
      <anchor>ab9b2c3a37a5516614c06d0ab54b97cda</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>int</type>
      <name>eos_debug</name>
      <anchorfile>group__eval.html</anchorfile>
      <anchor>ga567530678f6260a1a649a5beca5da4c5</anchor>
      <arglist></arglist>
    </member>
  </compound>
  <compound kind="file">
    <name>fold_vars.h</name>
    <path>/home/mescalin/ronny/WORK/ViennaRNA/H/</path>
    <filename>fold__vars_8h</filename>
    <member kind="function">
      <type>void</type>
      <name>set_model_details</name>
      <anchorfile>fold__vars_8h.html</anchorfile>
      <anchor>a4c3257186a796182462f18a5480ac8b3</anchor>
      <arglist>(model_detailsT *md)</arglist>
    </member>
    <member kind="variable">
      <type>int</type>
      <name>fold_constrained</name>
      <anchorfile>fold__vars_8h.html</anchorfile>
      <anchor>a0afc287c2464866d94858c39175154af</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>int</type>
      <name>noLonelyPairs</name>
      <anchorfile>fold__vars_8h.html</anchorfile>
      <anchor>a097eccaabd6ae8b4fef83cccff85bb5d</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>int</type>
      <name>dangles</name>
      <anchorfile>fold__vars_8h.html</anchorfile>
      <anchor>a72b511ed1201f7e23ec437e468790d74</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>int</type>
      <name>noGU</name>
      <anchorfile>fold__vars_8h.html</anchorfile>
      <anchor>abf380d09e4f1ab94fc6af57cf0ad5d32</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>int</type>
      <name>no_closingGU</name>
      <anchorfile>fold__vars_8h.html</anchorfile>
      <anchor>aa8d1c7b92489179e1eafa562b7bdd259</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>int</type>
      <name>tetra_loop</name>
      <anchorfile>fold__vars_8h.html</anchorfile>
      <anchor>a4f6265bdf0ead7ff4628a360adbfd77e</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>int</type>
      <name>energy_set</name>
      <anchorfile>fold__vars_8h.html</anchorfile>
      <anchor>afb1ef1166da85092ae8a325e02dcae71</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>int</type>
      <name>circ</name>
      <anchorfile>fold__vars_8h.html</anchorfile>
      <anchor>af9202a1a09f5828dc731e2d9a10fa111</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>int</type>
      <name>csv</name>
      <anchorfile>fold__vars_8h.html</anchorfile>
      <anchor>af2763d55a74663a5e60652b8880baa5b</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>int</type>
      <name>oldAliEn</name>
      <anchorfile>fold__vars_8h.html</anchorfile>
      <anchor>ac408868ba00671cbc7d1d535105af045</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>int</type>
      <name>ribo</name>
      <anchorfile>fold__vars_8h.html</anchorfile>
      <anchor>a0656afca1d2853f9ee6591172f5638de</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>char *</type>
      <name>RibosumFile</name>
      <anchorfile>fold__vars_8h.html</anchorfile>
      <anchor>a5dbaa0cca2c8c82048a0f0e38e164944</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>char *</type>
      <name>nonstandards</name>
      <anchorfile>fold__vars_8h.html</anchorfile>
      <anchor>a2695d91cc535d09c2eae5c3884e2ec64</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>double</type>
      <name>temperature</name>
      <anchorfile>fold__vars_8h.html</anchorfile>
      <anchor>ab4b11c8d9c758430960896bc3fe82ead</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>int</type>
      <name>james_rule</name>
      <anchorfile>fold__vars_8h.html</anchorfile>
      <anchor>af349001ad3b4d008d0051d935b1b6261</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>int</type>
      <name>logML</name>
      <anchorfile>fold__vars_8h.html</anchorfile>
      <anchor>a80c3c5fd35e7479704cc91d2d0367743</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>int</type>
      <name>cut_point</name>
      <anchorfile>fold__vars_8h.html</anchorfile>
      <anchor>ab9b2c3a37a5516614c06d0ab54b97cda</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>bondT *</type>
      <name>base_pair</name>
      <anchorfile>fold__vars_8h.html</anchorfile>
      <anchor>a0244a629b5ab4f58b77590c3dfd130dc</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>double *</type>
      <name>pr</name>
      <anchorfile>fold__vars_8h.html</anchorfile>
      <anchor>a0f5757427fd5f2f79d6fca0081cd5a52</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>int *</type>
      <name>iindx</name>
      <anchorfile>fold__vars_8h.html</anchorfile>
      <anchor>a92089ae3a51b5d75a14ce9cc29cc8317</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>double</type>
      <name>pf_scale</name>
      <anchorfile>fold__vars_8h.html</anchorfile>
      <anchor>ad3b22044065acc6dee0af68931b52cfd</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>int</type>
      <name>do_backtrack</name>
      <anchorfile>fold__vars_8h.html</anchorfile>
      <anchor>ad512b5dd4dbec60faccfe137bb474489</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>char</type>
      <name>backtrack_type</name>
      <anchorfile>fold__vars_8h.html</anchorfile>
      <anchor>a83bdb43472a259c71e69fa9f70f420c3</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>int</type>
      <name>gquad</name>
      <anchorfile>fold__vars_8h.html</anchorfile>
      <anchor>a25f2bdcdf56e813d288845484a13d704</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>int</type>
      <name>canonicalBPonly</name>
      <anchorfile>fold__vars_8h.html</anchorfile>
      <anchor>a22ae821b8918930e20ffa3fa84802b4b</anchor>
      <arglist></arglist>
    </member>
  </compound>
  <compound kind="file">
    <name>gquad.h</name>
    <path>/home/mescalin/ronny/WORK/ViennaRNA/H/</path>
    <filename>gquad_8h</filename>
    <member kind="function">
      <type>int *</type>
      <name>get_gquad_matrix</name>
      <anchorfile>gquad_8h.html</anchorfile>
      <anchor>a8b0784c14fa1208d0aebbebdc1318b7a</anchor>
      <arglist>(short *S, paramT *P)</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>parse_gquad</name>
      <anchorfile>gquad_8h.html</anchorfile>
      <anchor>ae41763215b9c64d2a7b67f0df8a28078</anchor>
      <arglist>(const char *struc, int *L, int l[3])</arglist>
    </member>
    <member kind="function">
      <type>PRIVATE int</type>
      <name>backtrack_GQuad_IntLoop</name>
      <anchorfile>gquad_8h.html</anchorfile>
      <anchor>a54475a8eb898fa1e8af8ab5f5375f3be</anchor>
      <arglist>(int c, int i, int j, int type, short *S, int *ggg, int *index, int *p, int *q, paramT *P)</arglist>
    </member>
    <member kind="function">
      <type>PRIVATE int</type>
      <name>backtrack_GQuad_IntLoop_L</name>
      <anchorfile>gquad_8h.html</anchorfile>
      <anchor>a118ec7289f1936bd810be7fe50b98212</anchor>
      <arglist>(int c, int i, int j, int type, short *S, int **ggg, int maxdist, int *p, int *q, paramT *P)</arglist>
    </member>
  </compound>
  <compound kind="file">
    <name>inverse.h</name>
    <path>/home/mescalin/ronny/WORK/ViennaRNA/H/</path>
    <filename>inverse_8h</filename>
    <member kind="function">
      <type>float</type>
      <name>inverse_fold</name>
      <anchorfile>group__inverse__fold.html</anchorfile>
      <anchor>ga7af026de55d4babad879f2c92559cbbc</anchor>
      <arglist>(char *start, const char *target)</arglist>
    </member>
    <member kind="function">
      <type>float</type>
      <name>inverse_pf_fold</name>
      <anchorfile>group__inverse__fold.html</anchorfile>
      <anchor>gaeef52ecbf2a2450ad585a344f9826806</anchor>
      <arglist>(char *start, const char *target)</arglist>
    </member>
    <member kind="variable">
      <type>char *</type>
      <name>symbolset</name>
      <anchorfile>group__inverse__fold.html</anchorfile>
      <anchor>ga8f791e7740a5a28b9f6fafb4e60301d9</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>float</type>
      <name>final_cost</name>
      <anchorfile>group__inverse__fold.html</anchorfile>
      <anchor>ga7f17d3b169af048d32bb185039a9c09c</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>int</type>
      <name>give_up</name>
      <anchorfile>group__inverse__fold.html</anchorfile>
      <anchor>ga7ec4ba51f86e1717a1e174264e4a75ce</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>int</type>
      <name>inv_verbose</name>
      <anchorfile>group__inverse__fold.html</anchorfile>
      <anchor>gafcfc65fba01b9cca5946726ed9057a63</anchor>
      <arglist></arglist>
    </member>
  </compound>
  <compound kind="file">
    <name>Lfold.h</name>
    <path>/home/mescalin/ronny/WORK/ViennaRNA/H/</path>
    <filename>Lfold_8h</filename>
    <member kind="function">
      <type>float</type>
      <name>Lfold</name>
      <anchorfile>group__local__mfe__fold.html</anchorfile>
      <anchor>ga16e5a70e60835bb969eaecbe6482f1be</anchor>
      <arglist>(const char *string, char *structure, int maxdist)</arglist>
    </member>
    <member kind="function">
      <type>float</type>
      <name>Lfoldz</name>
      <anchorfile>group__local__mfe__fold.html</anchorfile>
      <anchor>gab6d79eecc180f586679f7b85cce5cbe9</anchor>
      <arglist>(const char *string, char *structure, int maxdist, int zsc, double min_z)</arglist>
    </member>
    <member kind="function">
      <type>float</type>
      <name>aliLfold</name>
      <anchorfile>group__local__consensus__fold.html</anchorfile>
      <anchor>ga20a173a3cdb83f5d1778e36c1a6b1f2b</anchor>
      <arglist>(const char **strings, char *structure, int maxdist)</arglist>
    </member>
  </compound>
  <compound kind="file">
    <name>loop_energies.h</name>
    <path>/home/mescalin/ronny/WORK/ViennaRNA/H/</path>
    <filename>loop__energies_8h</filename>
    <member kind="function">
      <type>PRIVATE int</type>
      <name>E_IntLoop</name>
      <anchorfile>loop__energies_8h.html</anchorfile>
      <anchor>a3e5ad89f451254b1fe366d77aa8ff7bd</anchor>
      <arglist>(int n1, int n2, int type, int type_2, int si1, int sj1, int sp1, int sq1, paramT *P)</arglist>
    </member>
    <member kind="function">
      <type>PRIVATE int</type>
      <name>E_Hairpin</name>
      <anchorfile>loop__energies_8h.html</anchorfile>
      <anchor>aa362183cf6db89a10cdb0f5c4bd180c6</anchor>
      <arglist>(int size, int type, int si1, int sj1, const char *string, paramT *P)</arglist>
    </member>
    <member kind="function">
      <type>PRIVATE int</type>
      <name>E_Stem</name>
      <anchorfile>loop__energies_8h.html</anchorfile>
      <anchor>af5a6594eba9b2622cb47076650c69819</anchor>
      <arglist>(int type, int si1, int sj1, int extLoop, paramT *P)</arglist>
    </member>
    <member kind="function">
      <type>PRIVATE double</type>
      <name>exp_E_Stem</name>
      <anchorfile>loop__energies_8h.html</anchorfile>
      <anchor>a76cc24ec96199e04beddad13e7891e21</anchor>
      <arglist>(int type, int si1, int sj1, int extLoop, pf_paramT *P)</arglist>
    </member>
    <member kind="function">
      <type>PRIVATE double</type>
      <name>exp_E_Hairpin</name>
      <anchorfile>loop__energies_8h.html</anchorfile>
      <anchor>a0e128184bb097dc2da33706f33b555a6</anchor>
      <arglist>(int u, int type, short si1, short sj1, const char *string, pf_paramT *P)</arglist>
    </member>
    <member kind="function">
      <type>PRIVATE double</type>
      <name>exp_E_IntLoop</name>
      <anchorfile>loop__energies_8h.html</anchorfile>
      <anchor>aa5e98e524e2a41e290b942b09544bc9e</anchor>
      <arglist>(int u1, int u2, int type, int type2, short si1, short sj1, short sp1, short sq1, pf_paramT *P)</arglist>
    </member>
  </compound>
  <compound kind="file">
    <name>LPfold.h</name>
    <path>/home/mescalin/ronny/WORK/ViennaRNA/H/</path>
    <filename>LPfold_8h</filename>
    <member kind="function">
      <type>void</type>
      <name>update_pf_paramsLP</name>
      <anchorfile>group__local__pf__fold.html</anchorfile>
      <anchor>ga5a019014d37fe6105131dfc2fc447880</anchor>
      <arglist>(int length)</arglist>
    </member>
    <member kind="function">
      <type>plist *</type>
      <name>pfl_fold</name>
      <anchorfile>group__local__pf__fold.html</anchorfile>
      <anchor>gaa1ecd401617ebc748a0220026543c777</anchor>
      <arglist>(char *sequence, int winSize, int pairSize, float cutoffb, double **pU, struct plist **dpp2, FILE *pUfp, FILE *spup)</arglist>
    </member>
    <member kind="function">
      <type>plist *</type>
      <name>pfl_fold_par</name>
      <anchorfile>group__local__pf__fold.html</anchorfile>
      <anchor>gab354507e8028f3e1c52ef96bb1eb9df8</anchor>
      <arglist>(char *sequence, int winSize, int pairSize, float cutoffb, double **pU, struct plist **dpp2, FILE *pUfp, FILE *spup, pf_paramT *parameters)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>putoutpU_prob</name>
      <anchorfile>group__local__pf__fold.html</anchorfile>
      <anchor>ga0bcb751860bbf34e3dfee8c2fbdb3ef3</anchor>
      <arglist>(double **pU, int length, int ulength, FILE *fp, int energies)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>putoutpU_prob_bin</name>
      <anchorfile>group__local__pf__fold.html</anchorfile>
      <anchor>ga9acb00ee10e96b1ca4ea394cd8bcec75</anchor>
      <arglist>(double **pU, int length, int ulength, FILE *fp, int energies)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>init_pf_foldLP</name>
      <anchorfile>LPfold_8h.html</anchorfile>
      <anchor>ae85bf55053e9fb295208be322e0fa07a</anchor>
      <arglist>(int length)</arglist>
    </member>
  </compound>
  <compound kind="file">
    <name>MEA.h</name>
    <path>/home/mescalin/ronny/WORK/ViennaRNA/H/</path>
    <filename>MEA_8h</filename>
    <member kind="function">
      <type>float</type>
      <name>MEA</name>
      <anchorfile>MEA_8h.html</anchorfile>
      <anchor>a396ec6144c6a74fcbab4cea6b42d76c3</anchor>
      <arglist>(plist *p, char *structure, double gamma)</arglist>
    </member>
  </compound>
  <compound kind="file">
    <name>mm.h</name>
    <path>/home/mescalin/ronny/WORK/ViennaRNA/H/</path>
    <filename>mm_8h</filename>
  </compound>
  <compound kind="file">
    <name>naview.h</name>
    <path>/home/mescalin/ronny/WORK/ViennaRNA/H/</path>
    <filename>naview_8h</filename>
  </compound>
  <compound kind="file">
    <name>params.h</name>
    <path>/home/mescalin/ronny/WORK/ViennaRNA/H/</path>
    <filename>params_8h</filename>
    <member kind="function">
      <type>paramT *</type>
      <name>scale_parameters</name>
      <anchorfile>group__energy__parameters.html</anchorfile>
      <anchor>ga527ef619cd8210b84d5d53be1e0e29b6</anchor>
      <arglist>(void)</arglist>
    </member>
    <member kind="function">
      <type>paramT *</type>
      <name>get_scaled_parameters</name>
      <anchorfile>group__energy__parameters.html</anchorfile>
      <anchor>gac2f3ca440b7eaf4d999fb27da949fe72</anchor>
      <arglist>(double temperature, model_detailsT md)</arglist>
    </member>
    <member kind="function">
      <type>pf_paramT *</type>
      <name>get_scaled_pf_parameters</name>
      <anchorfile>group__energy__parameters.html</anchorfile>
      <anchor>gab85f6b6da051f380371deb0d8921bdba</anchor>
      <arglist>(void)</arglist>
    </member>
    <member kind="function">
      <type>pf_paramT *</type>
      <name>get_boltzmann_factors</name>
      <anchorfile>group__energy__parameters.html</anchorfile>
      <anchor>ga6fc2f3eef5a3024d44963ac59a42e39d</anchor>
      <arglist>(double temperature, double betaScale, model_detailsT md, double pf_scale)</arglist>
    </member>
    <member kind="function">
      <type>pf_paramT *</type>
      <name>get_boltzmann_factor_copy</name>
      <anchorfile>group__energy__parameters.html</anchorfile>
      <anchor>gacba212326a051734797e65987260fdd0</anchor>
      <arglist>(pf_paramT *parameters)</arglist>
    </member>
    <member kind="function">
      <type>pf_paramT *</type>
      <name>get_scaled_alipf_parameters</name>
      <anchorfile>group__energy__parameters.html</anchorfile>
      <anchor>gaa6a4297a2b91d6f7ae47dd61ca1862a0</anchor>
      <arglist>(unsigned int n_seq)</arglist>
    </member>
    <member kind="function">
      <type>PUBLIC pf_paramT *</type>
      <name>get_boltzmann_factors_ali</name>
      <anchorfile>group__energy__parameters.html</anchorfile>
      <anchor>gaaa049a8c9f1c2ed4398cb1b5a3d65a66</anchor>
      <arglist>(unsigned int n_seq, double temperature, double betaScale, model_detailsT md, double pf_scale)</arglist>
    </member>
  </compound>
  <compound kind="file">
    <name>part_func.h</name>
    <path>/home/mescalin/ronny/WORK/ViennaRNA/H/</path>
    <filename>part__func_8h</filename>
    <member kind="function">
      <type>float</type>
      <name>pf_fold_par</name>
      <anchorfile>group__pf__fold.html</anchorfile>
      <anchor>ga1839c61275760944b3a007c41d5c0823</anchor>
      <arglist>(const char *sequence, char *structure, pf_paramT *parameters, int calculate_bppm, int is_constrained, int is_circular)</arglist>
    </member>
    <member kind="function">
      <type>float</type>
      <name>pf_fold</name>
      <anchorfile>group__pf__fold.html</anchorfile>
      <anchor>gadc3db3d98742427e7001a7fd36ef28c2</anchor>
      <arglist>(const char *sequence, char *structure)</arglist>
    </member>
    <member kind="function">
      <type>float</type>
      <name>pf_circ_fold</name>
      <anchorfile>group__pf__fold.html</anchorfile>
      <anchor>ga819ce5fca8984004ac81c4a3b04cb735</anchor>
      <arglist>(const char *sequence, char *structure)</arglist>
    </member>
    <member kind="function">
      <type>char *</type>
      <name>pbacktrack</name>
      <anchorfile>group__subopt__stochbt.html</anchorfile>
      <anchor>gac03ca6db186bb3bf0a2a326d7fb3ba03</anchor>
      <arglist>(char *sequence)</arglist>
    </member>
    <member kind="function">
      <type>char *</type>
      <name>pbacktrack_circ</name>
      <anchorfile>group__subopt__stochbt.html</anchorfile>
      <anchor>ga00474051204ac9ad576b3e45174d03ff</anchor>
      <arglist>(char *sequence)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>free_pf_arrays</name>
      <anchorfile>group__pf__fold.html</anchorfile>
      <anchor>gae73db3f49a94f0f72e067ecd12681dbd</anchor>
      <arglist>(void)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>update_pf_params</name>
      <anchorfile>group__pf__fold.html</anchorfile>
      <anchor>ga384e927890f9c034ff09fa66da102d28</anchor>
      <arglist>(int length)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>update_pf_params_par</name>
      <anchorfile>group__pf__fold.html</anchorfile>
      <anchor>ga0733527a94de3b79eee3c3c03c99c1bc</anchor>
      <arglist>(int length, pf_paramT *parameters)</arglist>
    </member>
    <member kind="function">
      <type>double *</type>
      <name>export_bppm</name>
      <anchorfile>group__pf__fold.html</anchorfile>
      <anchor>ga6d463707d5f64bdc4d21515b7dd9b115</anchor>
      <arglist>(void)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>assign_plist_from_pr</name>
      <anchorfile>group__pf__fold.html</anchorfile>
      <anchor>ga03e15e831a31b1154855ab47edbdb019</anchor>
      <arglist>(plist **pl, double *probs, int length, double cutoff)</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>get_pf_arrays</name>
      <anchorfile>group__pf__fold.html</anchorfile>
      <anchor>ga18607e79e106cad827f482eedd2f632e</anchor>
      <arglist>(short **S_p, short **S1_p, char **ptype_p, double **qb_p, double **qm_p, double **q1k_p, double **qln_p)</arglist>
    </member>
    <member kind="function">
      <type>double</type>
      <name>get_subseq_F</name>
      <anchorfile>part__func_8h.html</anchorfile>
      <anchor>a189e2a1ec6cc32c53ea72f7543b0441e</anchor>
      <arglist>(int i, int j)</arglist>
    </member>
    <member kind="function">
      <type>char *</type>
      <name>get_centroid_struct_pl</name>
      <anchorfile>group__centroid__fold.html</anchorfile>
      <anchor>ga9aba0ba1433a6d259331e0fe9fc4a9a6</anchor>
      <arglist>(int length, double *dist, plist *pl)</arglist>
    </member>
    <member kind="function">
      <type>char *</type>
      <name>get_centroid_struct_pr</name>
      <anchorfile>group__centroid__fold.html</anchorfile>
      <anchor>gacdabece4aa1e20c9eaa97acb4c4dcc38</anchor>
      <arglist>(int length, double *dist, double *pr)</arglist>
    </member>
    <member kind="function">
      <type>double</type>
      <name>mean_bp_distance</name>
      <anchorfile>group__pf__fold.html</anchorfile>
      <anchor>ga79cbc375af65f11609feb6b055269e7d</anchor>
      <arglist>(int length)</arglist>
    </member>
    <member kind="function">
      <type>double</type>
      <name>mean_bp_distance_pr</name>
      <anchorfile>group__pf__fold.html</anchorfile>
      <anchor>ga72d84525f0afd3a9d60d830a2f501fa5</anchor>
      <arglist>(int length, double *pr)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>bppm_to_structure</name>
      <anchorfile>part__func_8h.html</anchorfile>
      <anchor>a1f562d463c14d4703d9656056200eb38</anchor>
      <arglist>(char *structure, double *pr, unsigned int length)</arglist>
    </member>
    <member kind="function">
      <type>char</type>
      <name>bppm_symbol</name>
      <anchorfile>part__func_8h.html</anchorfile>
      <anchor>a49962ad6242b8c628de6ca16bb831c1d</anchor>
      <arglist>(const float *x)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>init_pf_fold</name>
      <anchorfile>part__func_8h.html</anchorfile>
      <anchor>a15176e23eceeff8c7d14eabcfec8a2af</anchor>
      <arglist>(int length)</arglist>
    </member>
    <member kind="function">
      <type>char *</type>
      <name>centroid</name>
      <anchorfile>part__func_8h.html</anchorfile>
      <anchor>ae89a63bd83e75a80b2ba36d20b31ce81</anchor>
      <arglist>(int length, double *dist)</arglist>
    </member>
    <member kind="function">
      <type>double</type>
      <name>mean_bp_dist</name>
      <anchorfile>part__func_8h.html</anchorfile>
      <anchor>ae9556ba7ded44fe2321b6f67c3fc02a3</anchor>
      <arglist>(int length)</arglist>
    </member>
    <member kind="function">
      <type>double</type>
      <name>expLoopEnergy</name>
      <anchorfile>part__func_8h.html</anchorfile>
      <anchor>a68ba6f3a48e08ca131ab54621ce3a2d7</anchor>
      <arglist>(int u1, int u2, int type, int type2, short si1, short sj1, short sp1, short sq1)</arglist>
    </member>
    <member kind="function">
      <type>double</type>
      <name>expHairpinEnergy</name>
      <anchorfile>part__func_8h.html</anchorfile>
      <anchor>a7b6ab474cc80accc48010ccfcc59f96b</anchor>
      <arglist>(int u, int type, short si1, short sj1, const char *string)</arglist>
    </member>
    <member kind="variable">
      <type>int</type>
      <name>st_back</name>
      <anchorfile>group__subopt__stochbt.html</anchorfile>
      <anchor>gacd79b1a570e6ad9be24cb11fe8cae30a</anchor>
      <arglist></arglist>
    </member>
  </compound>
  <compound kind="file">
    <name>part_func_co.h</name>
    <path>/home/mescalin/ronny/WORK/ViennaRNA/H/</path>
    <filename>part__func__co_8h</filename>
    <member kind="function">
      <type>cofoldF</type>
      <name>co_pf_fold</name>
      <anchorfile>group__pf__cofold.html</anchorfile>
      <anchor>gaa86a5f998789ed71813d23d7307a791b</anchor>
      <arglist>(char *sequence, char *structure)</arglist>
    </member>
    <member kind="function">
      <type>cofoldF</type>
      <name>co_pf_fold_par</name>
      <anchorfile>group__pf__cofold.html</anchorfile>
      <anchor>gabd873b450832ab5f21101fc5ab354d21</anchor>
      <arglist>(char *sequence, char *structure, pf_paramT *parameters, int calculate_bppm, int is_constrained)</arglist>
    </member>
    <member kind="function">
      <type>double *</type>
      <name>export_co_bppm</name>
      <anchorfile>group__pf__cofold.html</anchorfile>
      <anchor>ga11f0252c1d2c4697253ff4b5bd392d3c</anchor>
      <arglist>(void)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>free_co_pf_arrays</name>
      <anchorfile>group__pf__cofold.html</anchorfile>
      <anchor>gade3ce34ae8214811374b1d28a40dc247</anchor>
      <arglist>(void)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>update_co_pf_params</name>
      <anchorfile>group__pf__cofold.html</anchorfile>
      <anchor>ga6e0f36c1f9b7d9dd4bfbad914c1119e5</anchor>
      <arglist>(int length)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>update_co_pf_params_par</name>
      <anchorfile>group__pf__cofold.html</anchorfile>
      <anchor>ga117d880df45bef444d5e2785ffa40a53</anchor>
      <arglist>(int length, pf_paramT *parameters)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>compute_probabilities</name>
      <anchorfile>group__pf__cofold.html</anchorfile>
      <anchor>ga15ae04ac5ab84e876dcf0093120cb617</anchor>
      <arglist>(double FAB, double FEA, double FEB, struct plist *prAB, struct plist *prA, struct plist *prB, int Alength)</arglist>
    </member>
    <member kind="function">
      <type>ConcEnt *</type>
      <name>get_concentrations</name>
      <anchorfile>group__pf__cofold.html</anchorfile>
      <anchor>ga5545cb936ac4ff93c7d699d46e72e8c7</anchor>
      <arglist>(double FEAB, double FEAA, double FEBB, double FEA, double FEB, double *startconc)</arglist>
    </member>
    <member kind="function">
      <type>plist *</type>
      <name>get_plist</name>
      <anchorfile>part__func__co_8h.html</anchorfile>
      <anchor>a334de3c96e2186abfbdc0eaea6d08b14</anchor>
      <arglist>(struct plist *pl, int length, double cut_off)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>init_co_pf_fold</name>
      <anchorfile>part__func__co_8h.html</anchorfile>
      <anchor>aa12dda9dd6179cdd22bcce87c0682c07</anchor>
      <arglist>(int length)</arglist>
    </member>
    <member kind="variable">
      <type>int</type>
      <name>mirnatog</name>
      <anchorfile>group__pf__cofold.html</anchorfile>
      <anchor>gaff27888c4088cc1f60fd59cbd589474c</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>double</type>
      <name>F_monomer</name>
      <anchorfile>group__pf__cofold.html</anchorfile>
      <anchor>gac2d1851a710a8561390861155ca988fe</anchor>
      <arglist>[2]</arglist>
    </member>
  </compound>
  <compound kind="file">
    <name>part_func_up.h</name>
    <path>/home/mescalin/ronny/WORK/ViennaRNA/H/</path>
    <filename>part__func__up_8h</filename>
    <member kind="function">
      <type>pu_contrib *</type>
      <name>pf_unstru</name>
      <anchorfile>group__up__cofold.html</anchorfile>
      <anchor>ga5b4ee40e190d2f633cd01cf0d2fe93cf</anchor>
      <arglist>(char *sequence, int max_w)</arglist>
    </member>
    <member kind="function">
      <type>interact *</type>
      <name>pf_interact</name>
      <anchorfile>group__up__cofold.html</anchorfile>
      <anchor>ga1aa0aa02bc3a724f87360c03097afd00</anchor>
      <arglist>(const char *s1, const char *s2, pu_contrib *p_c, pu_contrib *p_c2, int max_w, char *cstruc, int incr3, int incr5)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>free_interact</name>
      <anchorfile>group__up__cofold.html</anchorfile>
      <anchor>gadde308fd5f696dc271b1532aa96fd12f</anchor>
      <arglist>(interact *pin)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>free_pu_contrib_struct</name>
      <anchorfile>group__up__cofold.html</anchorfile>
      <anchor>gac20bd61824981d45ce0dc9934aa56df8</anchor>
      <arglist>(pu_contrib *pu)</arglist>
    </member>
  </compound>
  <compound kind="file">
    <name>plot_layouts.h</name>
    <path>/home/mescalin/ronny/WORK/ViennaRNA/H/</path>
    <filename>plot__layouts_8h</filename>
    <member kind="define">
      <type>#define</type>
      <name>VRNA_PLOT_TYPE_SIMPLE</name>
      <anchorfile>plot__layouts_8h.html</anchorfile>
      <anchor>ae6d17b9f0a53cf5205a9181e0f8422e9</anchor>
      <arglist></arglist>
    </member>
    <member kind="define">
      <type>#define</type>
      <name>VRNA_PLOT_TYPE_NAVIEW</name>
      <anchorfile>plot__layouts_8h.html</anchorfile>
      <anchor>a94d4c863ecac2f220f76658afb92f964</anchor>
      <arglist></arglist>
    </member>
    <member kind="define">
      <type>#define</type>
      <name>VRNA_PLOT_TYPE_CIRCULAR</name>
      <anchorfile>plot__layouts_8h.html</anchorfile>
      <anchor>a8c9eac631348da92136c8363ecdd9fb9</anchor>
      <arglist></arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>simple_xy_coordinates</name>
      <anchorfile>plot__layouts_8h.html</anchorfile>
      <anchor>af4b9173e7d3fd361c3c85e6def194123</anchor>
      <arglist>(short *pair_table, float *X, float *Y)</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>simple_circplot_coordinates</name>
      <anchorfile>plot__layouts_8h.html</anchorfile>
      <anchor>ac4ea13d35308f09940178d2b05a248c2</anchor>
      <arglist>(short *pair_table, float *x, float *y)</arglist>
    </member>
    <member kind="variable">
      <type>int</type>
      <name>rna_plot_type</name>
      <anchorfile>plot__layouts_8h.html</anchorfile>
      <anchor>a5964c4581431b098b80027d6e14dcdd4</anchor>
      <arglist></arglist>
    </member>
  </compound>
  <compound kind="file">
    <name>profiledist.h</name>
    <path>/home/mescalin/ronny/WORK/ViennaRNA/H/</path>
    <filename>profiledist_8h</filename>
    <member kind="function">
      <type>float</type>
      <name>profile_edit_distance</name>
      <anchorfile>profiledist_8h.html</anchorfile>
      <anchor>abe75e90e00a1e5dd8862944ed53dad5d</anchor>
      <arglist>(const float *T1, const float *T2)</arglist>
    </member>
    <member kind="function">
      <type>float *</type>
      <name>Make_bp_profile_bppm</name>
      <anchorfile>profiledist_8h.html</anchorfile>
      <anchor>a8822fd5268be115c6e6cdc92009436cc</anchor>
      <arglist>(double *bppm, int length)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>print_bppm</name>
      <anchorfile>profiledist_8h.html</anchorfile>
      <anchor>a8e0b4fe3698b3502945116ecc0ba6160</anchor>
      <arglist>(const float *T)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>free_profile</name>
      <anchorfile>profiledist_8h.html</anchorfile>
      <anchor>a9b0b84a5a45761bf42d7c835dcdb3b85</anchor>
      <arglist>(float *T)</arglist>
    </member>
    <member kind="function">
      <type>float *</type>
      <name>Make_bp_profile</name>
      <anchorfile>profiledist_8h.html</anchorfile>
      <anchor>a904c7eaf4a2413567c00ac4891749d18</anchor>
      <arglist>(int length)</arglist>
    </member>
  </compound>
  <compound kind="file">
    <name>PS_dot.h</name>
    <path>/home/mescalin/ronny/WORK/ViennaRNA/H/</path>
    <filename>PS__dot_8h</filename>
    <member kind="function">
      <type>int</type>
      <name>PS_rna_plot</name>
      <anchorfile>PS__dot_8h.html</anchorfile>
      <anchor>a0873c7cc4cd7a11c9a2cea19dde7e9c9</anchor>
      <arglist>(char *string, char *structure, char *file)</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>PS_rna_plot_a</name>
      <anchorfile>PS__dot_8h.html</anchorfile>
      <anchor>a47856b2504b566588785597b6ebb8271</anchor>
      <arglist>(char *string, char *structure, char *file, char *pre, char *post)</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>gmlRNA</name>
      <anchorfile>PS__dot_8h.html</anchorfile>
      <anchor>a70834bc8c0aad4fe6824ff76ccb8f329</anchor>
      <arglist>(char *string, char *structure, char *ssfile, char option)</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>ssv_rna_plot</name>
      <anchorfile>PS__dot_8h.html</anchorfile>
      <anchor>add368528755f9a830727b680243541df</anchor>
      <arglist>(char *string, char *structure, char *ssfile)</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>svg_rna_plot</name>
      <anchorfile>PS__dot_8h.html</anchorfile>
      <anchor>ae7853539b5df98f294b4af434e979304</anchor>
      <arglist>(char *string, char *structure, char *ssfile)</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>xrna_plot</name>
      <anchorfile>PS__dot_8h.html</anchorfile>
      <anchor>a2f6d5953e6a323df898896b8d6614483</anchor>
      <arglist>(char *string, char *structure, char *ssfile)</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>PS_dot_plot_list</name>
      <anchorfile>PS__dot_8h.html</anchorfile>
      <anchor>a00ea223b5cf02eb2faae5ff29f0d5e12</anchor>
      <arglist>(char *seq, char *filename, plist *pl, plist *mf, char *comment)</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>aliPS_color_aln</name>
      <anchorfile>PS__dot_8h.html</anchorfile>
      <anchor>aab48d4dac655d688abe921389ac2847c</anchor>
      <arglist>(const char *structure, const char *filename, const char *seqs[], const char *names[])</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>PS_dot_plot</name>
      <anchorfile>PS__dot_8h.html</anchorfile>
      <anchor>a689a97a7e3b8a2df14728b8204d9d57b</anchor>
      <arglist>(char *string, char *file)</arglist>
    </member>
  </compound>
  <compound kind="file">
    <name>read_epars.h</name>
    <path>/home/mescalin/ronny/WORK/ViennaRNA/H/</path>
    <filename>read__epars_8h</filename>
    <member kind="function">
      <type>void</type>
      <name>read_parameter_file</name>
      <anchorfile>group__energy__parameters__rw.html</anchorfile>
      <anchor>ga165a142a3c68fb6655c69ef4ab7cd749</anchor>
      <arglist>(const char fname[])</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>write_parameter_file</name>
      <anchorfile>group__energy__parameters__rw.html</anchorfile>
      <anchor>ga8a43459be386a7489feeab68dc2c6c76</anchor>
      <arglist>(const char fname[])</arglist>
    </member>
  </compound>
  <compound kind="file">
    <name>RNAstruct.h</name>
    <path>/home/mescalin/ronny/WORK/ViennaRNA/H/</path>
    <filename>RNAstruct_8h</filename>
    <member kind="function">
      <type>char *</type>
      <name>b2HIT</name>
      <anchorfile>RNAstruct_8h.html</anchorfile>
      <anchor>a07b7e90e712559a1992fba3ac6d21bbd</anchor>
      <arglist>(const char *structure)</arglist>
    </member>
    <member kind="function">
      <type>char *</type>
      <name>b2C</name>
      <anchorfile>RNAstruct_8h.html</anchorfile>
      <anchor>a9c80d92391f2833549a8b6dac92233f0</anchor>
      <arglist>(const char *structure)</arglist>
    </member>
    <member kind="function">
      <type>char *</type>
      <name>b2Shapiro</name>
      <anchorfile>RNAstruct_8h.html</anchorfile>
      <anchor>a5cd2feb367feeacad0c03cb7ddba5f10</anchor>
      <arglist>(const char *structure)</arglist>
    </member>
    <member kind="function">
      <type>char *</type>
      <name>add_root</name>
      <anchorfile>RNAstruct_8h.html</anchorfile>
      <anchor>a880d33066dd95441e5fbb73c57ed1c3e</anchor>
      <arglist>(const char *structure)</arglist>
    </member>
    <member kind="function">
      <type>char *</type>
      <name>expand_Shapiro</name>
      <anchorfile>RNAstruct_8h.html</anchorfile>
      <anchor>abe3d815b420dc4553bfb23511198b4c6</anchor>
      <arglist>(const char *coarse)</arglist>
    </member>
    <member kind="function">
      <type>char *</type>
      <name>expand_Full</name>
      <anchorfile>RNAstruct_8h.html</anchorfile>
      <anchor>a78d73cd54a068ef2812812771cdddc6f</anchor>
      <arglist>(const char *structure)</arglist>
    </member>
    <member kind="function">
      <type>char *</type>
      <name>unexpand_Full</name>
      <anchorfile>RNAstruct_8h.html</anchorfile>
      <anchor>a260c4b622093b76a883bf96628280de1</anchor>
      <arglist>(const char *ffull)</arglist>
    </member>
    <member kind="function">
      <type>char *</type>
      <name>unweight</name>
      <anchorfile>RNAstruct_8h.html</anchorfile>
      <anchor>a09a80253ac7b6bae606871ba7c6e5136</anchor>
      <arglist>(const char *wcoarse)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>unexpand_aligned_F</name>
      <anchorfile>RNAstruct_8h.html</anchorfile>
      <anchor>a1054c4477d53b31d79d4cb132100e87a</anchor>
      <arglist>(char *align[2])</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>parse_structure</name>
      <anchorfile>RNAstruct_8h.html</anchorfile>
      <anchor>a3c79042e6bf6f01706bf30ec9e69e8ac</anchor>
      <arglist>(const char *structure)</arglist>
    </member>
    <member kind="variable">
      <type>int</type>
      <name>loop_size</name>
      <anchorfile>RNAstruct_8h.html</anchorfile>
      <anchor>a3f31e0e48125601bfa57b52f8b038e8e</anchor>
      <arglist>[STRUC]</arglist>
    </member>
    <member kind="variable">
      <type>int</type>
      <name>helix_size</name>
      <anchorfile>RNAstruct_8h.html</anchorfile>
      <anchor>a8218c0d581a3fba2a1a56a196abe19a5</anchor>
      <arglist>[STRUC]</arglist>
    </member>
    <member kind="variable">
      <type>int</type>
      <name>loop_degree</name>
      <anchorfile>RNAstruct_8h.html</anchorfile>
      <anchor>aef14e2f8ab3f61e8e659ba6b9003b08a</anchor>
      <arglist>[STRUC]</arglist>
    </member>
    <member kind="variable">
      <type>int</type>
      <name>loops</name>
      <anchorfile>RNAstruct_8h.html</anchorfile>
      <anchor>a439fcb9f8d4f9f4d2227fde5fbfecb30</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>int</type>
      <name>unpaired</name>
      <anchorfile>RNAstruct_8h.html</anchorfile>
      <anchor>add2f952597e02d66e1116a9d11d252d6</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>int</type>
      <name>pairs</name>
      <anchorfile>RNAstruct_8h.html</anchorfile>
      <anchor>a6341cbb704924824e0236c1dce791032</anchor>
      <arglist></arglist>
    </member>
  </compound>
  <compound kind="file">
    <name>stringdist.h</name>
    <path>/home/mescalin/ronny/WORK/ViennaRNA/H/</path>
    <filename>stringdist_8h</filename>
    <member kind="function">
      <type>swString *</type>
      <name>Make_swString</name>
      <anchorfile>stringdist_8h.html</anchorfile>
      <anchor>a3125991b3a403b3f89230474deb3f22e</anchor>
      <arglist>(char *string)</arglist>
    </member>
    <member kind="function">
      <type>float</type>
      <name>string_edit_distance</name>
      <anchorfile>stringdist_8h.html</anchorfile>
      <anchor>a89e3c335ef17780576d7c0e713830db9</anchor>
      <arglist>(swString *T1, swString *T2)</arglist>
    </member>
  </compound>
  <compound kind="file">
    <name>subopt.h</name>
    <path>/home/mescalin/ronny/WORK/ViennaRNA/H/</path>
    <filename>subopt_8h</filename>
    <member kind="function">
      <type>SOLUTION *</type>
      <name>subopt</name>
      <anchorfile>group__subopt__wuchty.html</anchorfile>
      <anchor>ga700f662506a233e42dd7fda74fafd40e</anchor>
      <arglist>(char *seq, char *structure, int delta, FILE *fp)</arglist>
    </member>
    <member kind="function">
      <type>SOLUTION *</type>
      <name>subopt_par</name>
      <anchorfile>group__subopt__wuchty.html</anchorfile>
      <anchor>ga554dedfcdb249fdf151caade58666e4d</anchor>
      <arglist>(char *seq, char *structure, paramT *parameters, int delta, int is_constrained, int is_circular, FILE *fp)</arglist>
    </member>
    <member kind="function">
      <type>SOLUTION *</type>
      <name>subopt_circ</name>
      <anchorfile>group__subopt__wuchty.html</anchorfile>
      <anchor>ga8634516e4740e0b6c9a46d2bae940340</anchor>
      <arglist>(char *seq, char *sequence, int delta, FILE *fp)</arglist>
    </member>
    <member kind="variable">
      <type>int</type>
      <name>subopt_sorted</name>
      <anchorfile>group__subopt__wuchty.html</anchorfile>
      <anchor>ga873cf8ed69e0437f8efa8b1fec854a0e</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>double</type>
      <name>print_energy</name>
      <anchorfile>group__subopt__wuchty.html</anchorfile>
      <anchor>ga5e57d914bcb5feeecdf520e25313fcfe</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>int</type>
      <name>density_of_states</name>
      <anchorfile>group__dos.html</anchorfile>
      <anchor>ga937634a76b46a22530a74906f1957a9e</anchor>
      <arglist>[MAXDOS+1]</arglist>
    </member>
  </compound>
  <compound kind="file">
    <name>treedist.h</name>
    <path>/home/mescalin/ronny/WORK/ViennaRNA/H/</path>
    <filename>treedist_8h</filename>
    <member kind="function">
      <type>Tree *</type>
      <name>make_tree</name>
      <anchorfile>treedist_8h.html</anchorfile>
      <anchor>a08fe4d5afd385dce593b86eaf010c6e3</anchor>
      <arglist>(char *struc)</arglist>
    </member>
    <member kind="function">
      <type>float</type>
      <name>tree_edit_distance</name>
      <anchorfile>treedist_8h.html</anchorfile>
      <anchor>a3b21f1925f7071f46d93431a835217bb</anchor>
      <arglist>(Tree *T1, Tree *T2)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>print_tree</name>
      <anchorfile>treedist_8h.html</anchorfile>
      <anchor>a21ad4de3ba4055aeef08b28c9ad48894</anchor>
      <arglist>(Tree *t)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>free_tree</name>
      <anchorfile>treedist_8h.html</anchorfile>
      <anchor>acbc1cb9bce582ea945e4a467c76a57aa</anchor>
      <arglist>(Tree *t)</arglist>
    </member>
  </compound>
  <compound kind="file">
    <name>utils.h</name>
    <path>/home/mescalin/ronny/WORK/ViennaRNA/H/</path>
    <filename>utils_8h</filename>
    <member kind="define">
      <type>#define</type>
      <name>VRNA_INPUT_ERROR</name>
      <anchorfile>utils_8h.html</anchorfile>
      <anchor>ad403c9ea58f1836689404c2931419c8c</anchor>
      <arglist></arglist>
    </member>
    <member kind="define">
      <type>#define</type>
      <name>VRNA_INPUT_QUIT</name>
      <anchorfile>utils_8h.html</anchorfile>
      <anchor>a72f3c6ca5c83d2b9baed2922d19c403d</anchor>
      <arglist></arglist>
    </member>
    <member kind="define">
      <type>#define</type>
      <name>VRNA_INPUT_MISC</name>
      <anchorfile>utils_8h.html</anchorfile>
      <anchor>a8e3241b321c9c1a78a69e59e2e019a71</anchor>
      <arglist></arglist>
    </member>
    <member kind="define">
      <type>#define</type>
      <name>VRNA_INPUT_FASTA_HEADER</name>
      <anchorfile>utils_8h.html</anchorfile>
      <anchor>a2f0d8069e93d3ac54d9320d6bdb8e7e7</anchor>
      <arglist></arglist>
    </member>
    <member kind="define">
      <type>#define</type>
      <name>VRNA_INPUT_SEQUENCE</name>
      <anchorfile>utils_8h.html</anchorfile>
      <anchor>a8566d6787972100e68b5a2a159b4cf45</anchor>
      <arglist></arglist>
    </member>
    <member kind="define">
      <type>#define</type>
      <name>VRNA_INPUT_CONSTRAINT</name>
      <anchorfile>utils_8h.html</anchorfile>
      <anchor>ac08a9df45b9721b97a47dbfe7a6e5f85</anchor>
      <arglist></arglist>
    </member>
    <member kind="define">
      <type>#define</type>
      <name>VRNA_INPUT_NO_TRUNCATION</name>
      <anchorfile>utils_8h.html</anchorfile>
      <anchor>a086742158293217a46ae2f71bb296937</anchor>
      <arglist></arglist>
    </member>
    <member kind="define">
      <type>#define</type>
      <name>VRNA_INPUT_NO_REST</name>
      <anchorfile>utils_8h.html</anchorfile>
      <anchor>a7a2e8c50a0c7ce82e60da1016e1367fd</anchor>
      <arglist></arglist>
    </member>
    <member kind="define">
      <type>#define</type>
      <name>VRNA_INPUT_NO_SPAN</name>
      <anchorfile>utils_8h.html</anchorfile>
      <anchor>a0de536599b881c787b0943a2671da476</anchor>
      <arglist></arglist>
    </member>
    <member kind="define">
      <type>#define</type>
      <name>VRNA_INPUT_NOSKIP_BLANK_LINES</name>
      <anchorfile>utils_8h.html</anchorfile>
      <anchor>ab4db885222b3b69608310d7c7e63e286</anchor>
      <arglist></arglist>
    </member>
    <member kind="define">
      <type>#define</type>
      <name>VRNA_INPUT_BLANK_LINE</name>
      <anchorfile>utils_8h.html</anchorfile>
      <anchor>a305474b93ccb79ae4c7754016a8ddd84</anchor>
      <arglist></arglist>
    </member>
    <member kind="define">
      <type>#define</type>
      <name>VRNA_INPUT_NOSKIP_COMMENTS</name>
      <anchorfile>utils_8h.html</anchorfile>
      <anchor>a0f6311f11bed1842e3a527ab27b294c6</anchor>
      <arglist></arglist>
    </member>
    <member kind="define">
      <type>#define</type>
      <name>VRNA_INPUT_COMMENT</name>
      <anchorfile>utils_8h.html</anchorfile>
      <anchor>af2062e0eeefffd3ed639af460b3d4fab</anchor>
      <arglist></arglist>
    </member>
    <member kind="define">
      <type>#define</type>
      <name>VRNA_CONSTRAINT_PIPE</name>
      <anchorfile>utils_8h.html</anchorfile>
      <anchor>a4e8d7120619b21df0309af425acbc9a2</anchor>
      <arglist></arglist>
    </member>
    <member kind="define">
      <type>#define</type>
      <name>VRNA_CONSTRAINT_DOT</name>
      <anchorfile>utils_8h.html</anchorfile>
      <anchor>a55e1d16fd693ae9ec8e987b0750da804</anchor>
      <arglist></arglist>
    </member>
    <member kind="define">
      <type>#define</type>
      <name>VRNA_CONSTRAINT_X</name>
      <anchorfile>utils_8h.html</anchorfile>
      <anchor>a077c56550c915d4516d84a5ed8d051f4</anchor>
      <arglist></arglist>
    </member>
    <member kind="define">
      <type>#define</type>
      <name>VRNA_CONSTRAINT_ANG_BRACK</name>
      <anchorfile>utils_8h.html</anchorfile>
      <anchor>a0512d790f738742cbdcf3f7c87b46f48</anchor>
      <arglist></arglist>
    </member>
    <member kind="define">
      <type>#define</type>
      <name>VRNA_CONSTRAINT_RND_BRACK</name>
      <anchorfile>utils_8h.html</anchorfile>
      <anchor>aa20bfca4bb2903c8548000a33d7bbb53</anchor>
      <arglist></arglist>
    </member>
    <member kind="define">
      <type>#define</type>
      <name>VRNA_CONSTRAINT_MULTILINE</name>
      <anchorfile>utils_8h.html</anchorfile>
      <anchor>a7d725ef525b29891abef3f1ed42599a4</anchor>
      <arglist></arglist>
    </member>
    <member kind="define">
      <type>#define</type>
      <name>VRNA_CONSTRAINT_NO_HEADER</name>
      <anchorfile>utils_8h.html</anchorfile>
      <anchor>a08d12a9a846ea593b7171d277c9f033f</anchor>
      <arglist></arglist>
    </member>
    <member kind="define">
      <type>#define</type>
      <name>VRNA_CONSTRAINT_ALL</name>
      <anchorfile>utils_8h.html</anchorfile>
      <anchor>a0a697f77a6fbb10f34e16fa68ed9e655</anchor>
      <arglist></arglist>
    </member>
    <member kind="define">
      <type>#define</type>
      <name>VRNA_CONSTRAINT_G</name>
      <anchorfile>utils_8h.html</anchorfile>
      <anchor>a99dc6b23dc4080a76e2ed1a81c20e94d</anchor>
      <arglist></arglist>
    </member>
    <member kind="define">
      <type>#define</type>
      <name>VRNA_OPTION_MULTILINE</name>
      <anchorfile>utils_8h.html</anchorfile>
      <anchor>abec89c09874528c6cb73140a4c3d86d7</anchor>
      <arglist></arglist>
    </member>
    <member kind="define">
      <type>#define</type>
      <name>MIN2</name>
      <anchorfile>utils_8h.html</anchorfile>
      <anchor>ae0b9cd0ce090bd69b951aa73e8fa4f7d</anchor>
      <arglist>(A, B)</arglist>
    </member>
    <member kind="define">
      <type>#define</type>
      <name>MAX2</name>
      <anchorfile>utils_8h.html</anchorfile>
      <anchor>a33297b3679c713b0c4d897cd0fe3b122</anchor>
      <arglist>(A, B)</arglist>
    </member>
    <member kind="define">
      <type>#define</type>
      <name>MIN3</name>
      <anchorfile>utils_8h.html</anchorfile>
      <anchor>a721b8d5f3abef17f10293f1f7f8c958e</anchor>
      <arglist>(A, B, C)</arglist>
    </member>
    <member kind="define">
      <type>#define</type>
      <name>MAX3</name>
      <anchorfile>utils_8h.html</anchorfile>
      <anchor>a8d577123d2e66d2b7d0bf9af6e172b93</anchor>
      <arglist>(A, B, C)</arglist>
    </member>
    <member kind="define">
      <type>#define</type>
      <name>XSTR</name>
      <anchorfile>utils_8h.html</anchorfile>
      <anchor>a03943706e48069237cd57f2d35ca987e</anchor>
      <arglist>(s)</arglist>
    </member>
    <member kind="define">
      <type>#define</type>
      <name>STR</name>
      <anchorfile>utils_8h.html</anchorfile>
      <anchor>a6388870e639eee9c0a69446876f1f8cc</anchor>
      <arglist>(s)</arglist>
    </member>
    <member kind="define">
      <type>#define</type>
      <name>FILENAME_MAX_LENGTH</name>
      <anchorfile>utils_8h.html</anchorfile>
      <anchor>afb228174279df9486a5cb56ac0bc79a3</anchor>
      <arglist></arglist>
    </member>
    <member kind="define">
      <type>#define</type>
      <name>FILENAME_ID_LENGTH</name>
      <anchorfile>utils_8h.html</anchorfile>
      <anchor>a33c3b1826b8e2739f09f111ec719ded5</anchor>
      <arglist></arglist>
    </member>
    <member kind="function">
      <type>void *</type>
      <name>space</name>
      <anchorfile>utils_8h.html</anchorfile>
      <anchor>ad7e1e137b3bf1f7108933d302a7f0177</anchor>
      <arglist>(unsigned size)</arglist>
    </member>
    <member kind="function">
      <type>void *</type>
      <name>xrealloc</name>
      <anchorfile>utils_8h.html</anchorfile>
      <anchor>a9037ada838835b1b9db41581a021b0c8</anchor>
      <arglist>(void *p, unsigned size)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>nrerror</name>
      <anchorfile>utils_8h.html</anchorfile>
      <anchor>a127ce946e56b5a5773781cabe68e38c5</anchor>
      <arglist>(const char message[])</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>warn_user</name>
      <anchorfile>utils_8h.html</anchorfile>
      <anchor>af2355fa8746f2f30fbe71db65dea3d51</anchor>
      <arglist>(const char message[])</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>init_rand</name>
      <anchorfile>utils_8h.html</anchorfile>
      <anchor>a8aaa6d9be6f803f496d9b97375c371f3</anchor>
      <arglist>(void)</arglist>
    </member>
    <member kind="function">
      <type>double</type>
      <name>urn</name>
      <anchorfile>utils_8h.html</anchorfile>
      <anchor>aaa328491c84996e445d027fde9800f2e</anchor>
      <arglist>(void)</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>int_urn</name>
      <anchorfile>utils_8h.html</anchorfile>
      <anchor>a68ff0849d44f62fe491800378a5ffcb4</anchor>
      <arglist>(int from, int to)</arglist>
    </member>
    <member kind="function">
      <type>char *</type>
      <name>time_stamp</name>
      <anchorfile>utils_8h.html</anchorfile>
      <anchor>a7afeb906cb36e9d77379eabc6907ac46</anchor>
      <arglist>(void)</arglist>
    </member>
    <member kind="function">
      <type>char *</type>
      <name>random_string</name>
      <anchorfile>utils_8h.html</anchorfile>
      <anchor>a1b95eac365a021572e1c37e5993a89be</anchor>
      <arglist>(int l, const char symbols[])</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>hamming</name>
      <anchorfile>utils_8h.html</anchorfile>
      <anchor>ad9dc7bfc9aa664dc6698f17ce07fc7e7</anchor>
      <arglist>(const char *s1, const char *s2)</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>hamming_bound</name>
      <anchorfile>utils_8h.html</anchorfile>
      <anchor>a96d3c36717d624514055ce201cab1542</anchor>
      <arglist>(const char *s1, const char *s2, int n)</arglist>
    </member>
    <member kind="function">
      <type>char *</type>
      <name>get_line</name>
      <anchorfile>utils_8h.html</anchorfile>
      <anchor>abe51806d14cff0789a8c1df7dbc45b71</anchor>
      <arglist>(FILE *fp)</arglist>
    </member>
    <member kind="function">
      <type>unsigned int</type>
      <name>get_input_line</name>
      <anchorfile>utils_8h.html</anchorfile>
      <anchor>a8ef1835eb83f542396f59f0b205965e5</anchor>
      <arglist>(char **string, unsigned int options)</arglist>
    </member>
    <member kind="function">
      <type>unsigned int</type>
      <name>read_record</name>
      <anchorfile>utils_8h.html</anchorfile>
      <anchor>afd194a69af9d92b5b0412a7627ac1595</anchor>
      <arglist>(char **header, char **sequence, char ***rest, unsigned int options)</arglist>
    </member>
    <member kind="function">
      <type>char *</type>
      <name>pack_structure</name>
      <anchorfile>utils_8h.html</anchorfile>
      <anchor>ac6dfa5e22928c087c6e09ff0054a7ced</anchor>
      <arglist>(const char *struc)</arglist>
    </member>
    <member kind="function">
      <type>char *</type>
      <name>unpack_structure</name>
      <anchorfile>utils_8h.html</anchorfile>
      <anchor>a071c6921efe1eb974f115ee6fefa3c39</anchor>
      <arglist>(const char *packed)</arglist>
    </member>
    <member kind="function">
      <type>short *</type>
      <name>make_pair_table</name>
      <anchorfile>utils_8h.html</anchorfile>
      <anchor>a89c32307ee50a0026f4a3131fac0845a</anchor>
      <arglist>(const char *structure)</arglist>
    </member>
    <member kind="function">
      <type>short *</type>
      <name>copy_pair_table</name>
      <anchorfile>utils_8h.html</anchorfile>
      <anchor>afeaa6d68eef3a99d0a7aa08aa91c6601</anchor>
      <arglist>(const short *pt)</arglist>
    </member>
    <member kind="function">
      <type>short *</type>
      <name>alimake_pair_table</name>
      <anchorfile>utils_8h.html</anchorfile>
      <anchor>a3c81b3967056c3888b8472b65fbb16f5</anchor>
      <arglist>(const char *structure)</arglist>
    </member>
    <member kind="function">
      <type>short *</type>
      <name>make_pair_table_snoop</name>
      <anchorfile>utils_8h.html</anchorfile>
      <anchor>a9aa3bf3b4346bb7fb88efc154dd07a79</anchor>
      <arglist>(const char *structure)</arglist>
    </member>
    <member kind="function">
      <type>int *</type>
      <name>make_loop_index_pt</name>
      <anchorfile>utils_8h.html</anchorfile>
      <anchor>a4358e89f64cc87a563b7ef3855f75bed</anchor>
      <arglist>(short *pt)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>print_tty_input_seq</name>
      <anchorfile>utils_8h.html</anchorfile>
      <anchor>a6bf778117d31b7fd90db435323f4ef74</anchor>
      <arglist>(void)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>print_tty_input_seq_str</name>
      <anchorfile>utils_8h.html</anchorfile>
      <anchor>ae4ef89b662a3e9b5b5f0781d9757aba0</anchor>
      <arglist>(const char *s)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>print_tty_constraint_full</name>
      <anchorfile>utils_8h.html</anchorfile>
      <anchor>ae8ae8a34962b9959be3f6c40f0a80ac1</anchor>
      <arglist>(void)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>print_tty_constraint</name>
      <anchorfile>utils_8h.html</anchorfile>
      <anchor>a4d167deb70bb51723e44374dc981deb2</anchor>
      <arglist>(unsigned int option)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>str_DNA2RNA</name>
      <anchorfile>utils_8h.html</anchorfile>
      <anchor>ad3f18dd83f958f18b2f26ecb99305208</anchor>
      <arglist>(char *sequence)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>str_uppercase</name>
      <anchorfile>utils_8h.html</anchorfile>
      <anchor>a17b796b806f96b70382077fb5bc519bb</anchor>
      <arglist>(char *sequence)</arglist>
    </member>
    <member kind="function">
      <type>int *</type>
      <name>get_iindx</name>
      <anchorfile>utils_8h.html</anchorfile>
      <anchor>a55c0f6b3b07b6adf2ee235ba901fe397</anchor>
      <arglist>(unsigned int length)</arglist>
    </member>
    <member kind="function">
      <type>int *</type>
      <name>get_indx</name>
      <anchorfile>utils_8h.html</anchorfile>
      <anchor>a4d9ee1572c1bfcd02d3d3f2db8a6530f</anchor>
      <arglist>(unsigned int length)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>constrain_ptypes</name>
      <anchorfile>utils_8h.html</anchorfile>
      <anchor>a36c3a6c3218b041f992052767bc74549</anchor>
      <arglist>(const char *constraint, unsigned int length, char *ptype, int *BP, int min_loop_size, unsigned int idx_type)</arglist>
    </member>
    <member kind="variable">
      <type>unsigned short</type>
      <name>xsubi</name>
      <anchorfile>utils_8h.html</anchorfile>
      <anchor>af9a866c8417afda7368bbac939ab3c47</anchor>
      <arglist>[3]</arglist>
    </member>
  </compound>
  <compound kind="file">
    <name>1.8.4_epars.h</name>
    <path>/home/mescalin/ronny/WORK/ViennaRNA/lib/</path>
    <filename>1_88_84__epars_8h</filename>
  </compound>
  <compound kind="file">
    <name>1.8.4_intloops.h</name>
    <path>/home/mescalin/ronny/WORK/ViennaRNA/lib/</path>
    <filename>1_88_84__intloops_8h</filename>
  </compound>
  <compound kind="page">
    <name>mp_parse</name>
    <title>Parsing and Comparing - Functions to Manipulate Structures</title>
    <filename>mp_parse</filename>
  </compound>
  <compound kind="page">
    <name>mp_utils</name>
    <title>Utilities - Odds and Ends</title>
    <filename>mp_utils</filename>
    <docanchor file="mp_utils">toc</docanchor>
    <docanchor file="mp_utils" title="Producing secondary structure graphs">utils_ss</docanchor>
    <docanchor file="mp_utils" title="Producing (colored) dot plots for base pair probabilities">utils_dot</docanchor>
    <docanchor file="mp_utils" title="Producing (colored) alignments">utils_aln</docanchor>
    <docanchor file="mp_utils" title="RNA sequence related utilities">utils_seq</docanchor>
    <docanchor file="mp_utils" title="RNA secondary structure related utilities">utils_struc</docanchor>
    <docanchor file="mp_utils" title="Miscellaneous Utilities">utils_misc</docanchor>
  </compound>
  <compound kind="page">
    <name>mp_example</name>
    <title>Example - A Small Example Program</title>
    <filename>mp_example</filename>
  </compound>
  <compound kind="group">
    <name>folding_routines</name>
    <title>RNA Secondary Structure Folding</title>
    <filename>group__folding__routines.html</filename>
    <subgroup>mfe_fold</subgroup>
    <subgroup>pf_fold</subgroup>
    <subgroup>subopt_fold</subgroup>
    <subgroup>cofold</subgroup>
    <subgroup>consensus_fold</subgroup>
    <subgroup>local_fold</subgroup>
    <subgroup>energy_parameters</subgroup>
    <subgroup>eval</subgroup>
    <subgroup>inverse_fold</subgroup>
    <subgroup>class_fold</subgroup>
  </compound>
  <compound kind="group">
    <name>mfe_fold</name>
    <title>Calculating Minimum Free Energy (MFE) Structures</title>
    <filename>group__mfe__fold.html</filename>
    <subgroup>mfe_cofold</subgroup>
    <subgroup>consensus_mfe_fold</subgroup>
    <subgroup>local_mfe_fold</subgroup>
    <subgroup>kl_neighborhood_mfe</subgroup>
  </compound>
  <compound kind="group">
    <name>pf_fold</name>
    <title>Calculating Partition Functions and Pair Probabilities</title>
    <filename>group__pf__fold.html</filename>
    <subgroup>mea_fold</subgroup>
    <subgroup>centroid_fold</subgroup>
    <subgroup>pf_cofold</subgroup>
    <subgroup>up_cofold</subgroup>
    <subgroup>consensus_pf_fold</subgroup>
    <subgroup>local_pf_fold</subgroup>
    <subgroup>kl_neighborhood_pf</subgroup>
    <file>part_func.h</file>
  </compound>
  <compound kind="group">
    <name>mea_fold</name>
    <title>Compute the structure with maximum expected accuracy (MEA)</title>
    <filename>group__mea__fold.html</filename>
  </compound>
  <compound kind="group">
    <name>centroid_fold</name>
    <title>Compute the centroid structure</title>
    <filename>group__centroid__fold.html</filename>
  </compound>
  <compound kind="group">
    <name>subopt_fold</name>
    <title>Enumerating Suboptimal Structures</title>
    <filename>group__subopt__fold.html</filename>
    <subgroup>subopt_zuker</subgroup>
    <subgroup>subopt_wuchty</subgroup>
    <subgroup>subopt_stochbt</subgroup>
    <file>subopt.h</file>
  </compound>
  <compound kind="group">
    <name>subopt_zuker</name>
    <title>Suboptimal structures according to Zuker et al. 1989</title>
    <filename>group__subopt__zuker.html</filename>
  </compound>
  <compound kind="group">
    <name>subopt_wuchty</name>
    <title>Suboptimal structures within an energy band arround the MFE</title>
    <filename>group__subopt__wuchty.html</filename>
  </compound>
  <compound kind="group">
    <name>subopt_stochbt</name>
    <title>Stochastic backtracking in the Ensemble</title>
    <filename>group__subopt__stochbt.html</filename>
    <subgroup>consensus_stochbt</subgroup>
    <subgroup>kl_neighborhood_stochbt</subgroup>
  </compound>
  <compound kind="group">
    <name>cofold</name>
    <title>Calculate Secondary Structures of two RNAs upon Dimerization</title>
    <filename>group__cofold.html</filename>
    <subgroup>mfe_cofold</subgroup>
    <subgroup>pf_cofold</subgroup>
    <subgroup>up_cofold</subgroup>
  </compound>
  <compound kind="group">
    <name>mfe_cofold</name>
    <title>MFE Structures of two hybridized Sequences</title>
    <filename>group__mfe__cofold.html</filename>
    <file>cofold.h</file>
  </compound>
  <compound kind="group">
    <name>pf_cofold</name>
    <title>Partition Function for two hybridized Sequences</title>
    <filename>group__pf__cofold.html</filename>
    <file>part_func_co.h</file>
  </compound>
  <compound kind="group">
    <name>up_cofold</name>
    <title>Partition Function for two hybridized Sequences as a stepwise Process</title>
    <filename>group__up__cofold.html</filename>
    <file>part_func_up.h</file>
  </compound>
  <compound kind="group">
    <name>consensus_fold</name>
    <title>Predicting Consensus Structures from Alignment(s)</title>
    <filename>group__consensus__fold.html</filename>
    <subgroup>consensus_mfe_fold</subgroup>
    <subgroup>consensus_pf_fold</subgroup>
    <subgroup>consensus_stochbt</subgroup>
    <subgroup>local_consensus_fold</subgroup>
    <file>alifold.h</file>
  </compound>
  <compound kind="group">
    <name>consensus_mfe_fold</name>
    <title>MFE Consensus Structures for Sequence Alignment(s)</title>
    <filename>group__consensus__mfe__fold.html</filename>
  </compound>
  <compound kind="group">
    <name>consensus_pf_fold</name>
    <title>Partition Function and Base Pair Probabilities for Sequence Alignment(s)</title>
    <filename>group__consensus__pf__fold.html</filename>
  </compound>
  <compound kind="group">
    <name>consensus_stochbt</name>
    <title>Stochastic Backtracking of Consensus Structures from Sequence Alignment(s)</title>
    <filename>group__consensus__stochbt.html</filename>
  </compound>
  <compound kind="group">
    <name>local_fold</name>
    <title>Predicting Locally stable structures of large sequences</title>
    <filename>group__local__fold.html</filename>
    <subgroup>local_mfe_fold</subgroup>
    <subgroup>local_pf_fold</subgroup>
    <subgroup>local_consensus_fold</subgroup>
    <file>Lfold.h</file>
  </compound>
  <compound kind="group">
    <name>local_mfe_fold</name>
    <title>Local MFE structure Prediction and Z-scores</title>
    <filename>group__local__mfe__fold.html</filename>
  </compound>
  <compound kind="group">
    <name>local_pf_fold</name>
    <title>Partition functions for locally stable secondary structures</title>
    <filename>group__local__pf__fold.html</filename>
    <file>LPfold.h</file>
  </compound>
  <compound kind="group">
    <name>local_consensus_fold</name>
    <title>Local MFE consensus structures for Sequence Alignments</title>
    <filename>group__local__consensus__fold.html</filename>
  </compound>
  <compound kind="group">
    <name>energy_parameters</name>
    <title>Change and Precalculate Energy Parameter Sets and Boltzmann Factors</title>
    <filename>group__energy__parameters.html</filename>
    <subgroup>energy_parameters_rw</subgroup>
    <file>params.h</file>
  </compound>
  <compound kind="group">
    <name>energy_parameters_rw</name>
    <title>Reading/Writing energy parameter sets from/to File</title>
    <filename>group__energy__parameters__rw.html</filename>
    <subgroup>energy_parameters_convert</subgroup>
    <file>read_epars.h</file>
  </compound>
  <compound kind="group">
    <name>energy_parameters_convert</name>
    <title>Converting energy parameter files</title>
    <filename>group__energy__parameters__convert.html</filename>
    <file>convert_epars.h</file>
  </compound>
  <compound kind="group">
    <name>eval</name>
    <title>Energy evaluation</title>
    <filename>group__eval.html</filename>
  </compound>
  <compound kind="group">
    <name>inverse_fold</name>
    <title>Searching Sequences for Predefined Structures</title>
    <filename>group__inverse__fold.html</filename>
    <file>inverse.h</file>
  </compound>
  <compound kind="group">
    <name>class_fold</name>
    <title>Classified Dynamic Programming</title>
    <filename>group__class__fold.html</filename>
    <subgroup>kl_neighborhood</subgroup>
    <subgroup>dos</subgroup>
  </compound>
  <compound kind="group">
    <name>kl_neighborhood</name>
    <title>Distance based partitioning of the Secondary Structure Space</title>
    <filename>group__kl__neighborhood.html</filename>
    <subgroup>kl_neighborhood_mfe</subgroup>
    <subgroup>kl_neighborhood_pf</subgroup>
    <subgroup>kl_neighborhood_stochbt</subgroup>
  </compound>
  <compound kind="group">
    <name>kl_neighborhood_mfe</name>
    <title>Calculating MFE representatives of a Distance Based Partitioning</title>
    <filename>group__kl__neighborhood__mfe.html</filename>
    <file>2Dfold.h</file>
  </compound>
  <compound kind="group">
    <name>kl_neighborhood_pf</name>
    <title>Calculate Partition Functions of a Distance Based Partitioning</title>
    <filename>group__kl__neighborhood__pf.html</filename>
    <file>2Dpfold.h</file>
  </compound>
  <compound kind="group">
    <name>kl_neighborhood_stochbt</name>
    <title>Stochastic Backtracking of Structures from Distance Based Partitioning</title>
    <filename>group__kl__neighborhood__stochbt.html</filename>
  </compound>
  <compound kind="group">
    <name>dos</name>
    <title>Compute the Density of States</title>
    <filename>group__dos.html</filename>
  </compound>
  <compound kind="group">
    <name>parse</name>
    <title>Parsing and Comparing - Functions to Manipulate Structures</title>
    <filename>group__parse.html</filename>
  </compound>
  <compound kind="struct">
    <name>bondT</name>
    <filename>structbondT.html</filename>
  </compound>
  <compound kind="struct">
    <name>bondTEn</name>
    <filename>structbondTEn.html</filename>
  </compound>
  <compound kind="struct">
    <name>cofoldF</name>
    <filename>structcofoldF.html</filename>
    <member kind="variable">
      <type>double</type>
      <name>F0AB</name>
      <anchorfile>structcofoldF.html</anchorfile>
      <anchor>af6c496438321eb8bb907a21de1915c23</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>double</type>
      <name>FAB</name>
      <anchorfile>structcofoldF.html</anchorfile>
      <anchor>a2ae1245ff4a93cd11f882f490f777cb7</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>double</type>
      <name>FcAB</name>
      <anchorfile>structcofoldF.html</anchorfile>
      <anchor>a4899a4f9b42e416baf46c5fe10751c45</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>double</type>
      <name>FA</name>
      <anchorfile>structcofoldF.html</anchorfile>
      <anchor>a460f3ba205c205e6f5ec27cc2e2eb2b2</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>double</type>
      <name>FB</name>
      <anchorfile>structcofoldF.html</anchorfile>
      <anchor>ad3e5466724f3987be9d6f388b8ee5129</anchor>
      <arglist></arglist>
    </member>
  </compound>
  <compound kind="struct">
    <name>ConcEnt</name>
    <filename>structConcEnt.html</filename>
    <member kind="variable">
      <type>double</type>
      <name>A0</name>
      <anchorfile>structConcEnt.html</anchorfile>
      <anchor>adcf4d93c7efeaa4e6c4154b64d367681</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>double</type>
      <name>B0</name>
      <anchorfile>structConcEnt.html</anchorfile>
      <anchor>add4c33b94b34e847fbf5838b04cce346</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>double</type>
      <name>ABc</name>
      <anchorfile>structConcEnt.html</anchorfile>
      <anchor>ac59c07a31d844e7b05bcdc05c4413b19</anchor>
      <arglist></arglist>
    </member>
  </compound>
  <compound kind="struct">
    <name>constrain</name>
    <filename>structconstrain.html</filename>
  </compound>
  <compound kind="struct">
    <name>COORDINATE</name>
    <filename>structCOORDINATE.html</filename>
  </compound>
  <compound kind="struct">
    <name>cpair</name>
    <filename>structcpair.html</filename>
  </compound>
  <compound kind="struct">
    <name>duplexT</name>
    <filename>structduplexT.html</filename>
  </compound>
  <compound kind="struct">
    <name>dupVar</name>
    <filename>structdupVar.html</filename>
  </compound>
  <compound kind="struct">
    <name>folden</name>
    <filename>structfolden.html</filename>
  </compound>
  <compound kind="struct">
    <name>interact</name>
    <filename>structinteract.html</filename>
    <member kind="variable">
      <type>double *</type>
      <name>Pi</name>
      <anchorfile>structinteract.html</anchorfile>
      <anchor>a1fc8b3860c083f164daa9712690a3a56</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>double *</type>
      <name>Gi</name>
      <anchorfile>structinteract.html</anchorfile>
      <anchor>a54f8183542fff4c32ab7ace49a16c02c</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>double</type>
      <name>Gikjl</name>
      <anchorfile>structinteract.html</anchorfile>
      <anchor>ad58303190f9e085c3ab59890cbf61223</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>double</type>
      <name>Gikjl_wo</name>
      <anchorfile>structinteract.html</anchorfile>
      <anchor>a41793812abae560805414761fec398fe</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>int</type>
      <name>i</name>
      <anchorfile>structinteract.html</anchorfile>
      <anchor>ab6d031a21388be8763b75ea74c937f17</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>int</type>
      <name>k</name>
      <anchorfile>structinteract.html</anchorfile>
      <anchor>a61e457fbf943d57364be6ddf1b4e7b8a</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>int</type>
      <name>j</name>
      <anchorfile>structinteract.html</anchorfile>
      <anchor>a7555cb6363d1479341eb72b9c087aa34</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>int</type>
      <name>l</name>
      <anchorfile>structinteract.html</anchorfile>
      <anchor>a030ab45056342e12cb3955e4defd3904</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>int</type>
      <name>length</name>
      <anchorfile>structinteract.html</anchorfile>
      <anchor>ac9fcb5dca54ec5faa76e02b6488b9524</anchor>
      <arglist></arglist>
    </member>
  </compound>
  <compound kind="struct">
    <name>intermediate_t</name>
    <filename>structintermediate__t.html</filename>
    <member kind="variable">
      <type>short *</type>
      <name>pt</name>
      <anchorfile>structintermediate__t.html</anchorfile>
      <anchor>a9a2b6258aa1af06ea3504631de8dadba</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>int</type>
      <name>Sen</name>
      <anchorfile>structintermediate__t.html</anchorfile>
      <anchor>ac44e091915da58927978d54ef59234c7</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>int</type>
      <name>curr_en</name>
      <anchorfile>structintermediate__t.html</anchorfile>
      <anchor>af84d640df33aea99e959b2e4f61a7367</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>move_t *</type>
      <name>moves</name>
      <anchorfile>structintermediate__t.html</anchorfile>
      <anchor>a94e947f18273bbfe3dd544085b025a7b</anchor>
      <arglist></arglist>
    </member>
  </compound>
  <compound kind="struct">
    <name>INTERVAL</name>
    <filename>structINTERVAL.html</filename>
  </compound>
  <compound kind="struct">
    <name>LIST</name>
    <filename>structLIST.html</filename>
  </compound>
  <compound kind="struct">
    <name>LST_BUCKET</name>
    <filename>structLST__BUCKET.html</filename>
  </compound>
  <compound kind="struct">
    <name>model_detailsT</name>
    <filename>structmodel__detailsT.html</filename>
    <member kind="variable">
      <type>int</type>
      <name>dangles</name>
      <anchorfile>structmodel__detailsT.html</anchorfile>
      <anchor>a89f9df217a4a7f4351a642655976376b</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>int</type>
      <name>special_hp</name>
      <anchorfile>structmodel__detailsT.html</anchorfile>
      <anchor>a9d73fde17b0465311a80f607faa85617</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>int</type>
      <name>noLP</name>
      <anchorfile>structmodel__detailsT.html</anchorfile>
      <anchor>a6fb076173d2cbc4259606ce23eedf17d</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>int</type>
      <name>noGU</name>
      <anchorfile>structmodel__detailsT.html</anchorfile>
      <anchor>a0f982e6904d012e4fe41e99daa797f5d</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>int</type>
      <name>noGUclosure</name>
      <anchorfile>structmodel__detailsT.html</anchorfile>
      <anchor>ae9cedf375cd904e5fb8e56cf3f64bcd9</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>int</type>
      <name>logML</name>
      <anchorfile>structmodel__detailsT.html</anchorfile>
      <anchor>a415556dc150e02d108be81ecc5c48e85</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>int</type>
      <name>circ</name>
      <anchorfile>structmodel__detailsT.html</anchorfile>
      <anchor>a55f083dad18c216505805a8062e63074</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>int</type>
      <name>gquad</name>
      <anchorfile>structmodel__detailsT.html</anchorfile>
      <anchor>aa1ecdce7bc3f375bd8a9a7b738abc0ea</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>int</type>
      <name>canonicalBPonly</name>
      <anchorfile>structmodel__detailsT.html</anchorfile>
      <anchor>a20e38e8e65afedb3be9c16ae03d956c0</anchor>
      <arglist></arglist>
    </member>
  </compound>
  <compound kind="struct">
    <name>move_t</name>
    <filename>structmove__t.html</filename>
  </compound>
  <compound kind="struct">
    <name>PAIR</name>
    <filename>structPAIR.html</filename>
  </compound>
  <compound kind="struct">
    <name>pair_info</name>
    <filename>structpair__info.html</filename>
    <member kind="variable">
      <type>unsigned</type>
      <name>i</name>
      <anchorfile>structpair__info.html</anchorfile>
      <anchor>a7905e17da6a6cc48230ee6205628ed7f</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>unsigned</type>
      <name>j</name>
      <anchorfile>structpair__info.html</anchorfile>
      <anchor>a35665817b5792703ff4325e1bcbe5e21</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>float</type>
      <name>p</name>
      <anchorfile>structpair__info.html</anchorfile>
      <anchor>af0895ea40ec0c23bfe8aa2c3babf0e80</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>float</type>
      <name>ent</name>
      <anchorfile>structpair__info.html</anchorfile>
      <anchor>ab3aa7a54e6976f46e69c6ffcddd0e782</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>short</type>
      <name>bp</name>
      <anchorfile>structpair__info.html</anchorfile>
      <anchor>a23fc316453d179474bed7f6ed2489723</anchor>
      <arglist>[8]</arglist>
    </member>
    <member kind="variable">
      <type>char</type>
      <name>comp</name>
      <anchorfile>structpair__info.html</anchorfile>
      <anchor>a4da3d6c9042500c16c4b06e0bbc48190</anchor>
      <arglist></arglist>
    </member>
  </compound>
  <compound kind="struct">
    <name>pairpro</name>
    <filename>structpairpro.html</filename>
  </compound>
  <compound kind="struct">
    <name>paramT</name>
    <filename>structparamT.html</filename>
    <member kind="variable">
      <type>double</type>
      <name>temperature</name>
      <anchorfile>structparamT.html</anchorfile>
      <anchor>a8ed207b95868d1085bd9c197fbc6924f</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>model_detailsT</type>
      <name>model_details</name>
      <anchorfile>structparamT.html</anchorfile>
      <anchor>aeb912822ef912705bc202b14f9d71ad9</anchor>
      <arglist></arglist>
    </member>
  </compound>
  <compound kind="struct">
    <name>path_t</name>
    <filename>structpath__t.html</filename>
  </compound>
  <compound kind="struct">
    <name>pf_paramT</name>
    <filename>structpf__paramT.html</filename>
    <member kind="variable">
      <type>double</type>
      <name>pf_scale</name>
      <anchorfile>structpf__paramT.html</anchorfile>
      <anchor>aef40322e7ca1adbd9b438aeda0352e8f</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>double</type>
      <name>temperature</name>
      <anchorfile>structpf__paramT.html</anchorfile>
      <anchor>aa0e11e9f1f6e212640baf40d7195a014</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>double</type>
      <name>alpha</name>
      <anchorfile>structpf__paramT.html</anchorfile>
      <anchor>a3d2af9040acfa08295efb50f0219149d</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>model_detailsT</type>
      <name>model_details</name>
      <anchorfile>structpf__paramT.html</anchorfile>
      <anchor>a43ec875779c5e7c8bf5fa7e837ec6d09</anchor>
      <arglist></arglist>
    </member>
  </compound>
  <compound kind="struct">
    <name>plist</name>
    <filename>structplist.html</filename>
  </compound>
  <compound kind="struct">
    <name>Postorder_list</name>
    <filename>structPostorder__list.html</filename>
  </compound>
  <compound kind="struct">
    <name>pu_contrib</name>
    <filename>structpu__contrib.html</filename>
    <member kind="variable">
      <type>double **</type>
      <name>H</name>
      <anchorfile>structpu__contrib.html</anchorfile>
      <anchor>ac9034ac9a84ed0647587659d6e9be1e8</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>double **</type>
      <name>I</name>
      <anchorfile>structpu__contrib.html</anchorfile>
      <anchor>a8ca0da20536780589fb3e3472ca0581f</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>double **</type>
      <name>M</name>
      <anchorfile>structpu__contrib.html</anchorfile>
      <anchor>a1222ebf74f426bbcd843dcc325da207b</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>double **</type>
      <name>E</name>
      <anchorfile>structpu__contrib.html</anchorfile>
      <anchor>accb192ba6b4b91a1cb2f8080934fd428</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>int</type>
      <name>length</name>
      <anchorfile>structpu__contrib.html</anchorfile>
      <anchor>a33d5ada6e861db0c81aa3d5b2989262e</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>int</type>
      <name>w</name>
      <anchorfile>structpu__contrib.html</anchorfile>
      <anchor>a403c1c7f20beeeffba7632fac0cfcbff</anchor>
      <arglist></arglist>
    </member>
  </compound>
  <compound kind="struct">
    <name>pu_out</name>
    <filename>structpu__out.html</filename>
    <member kind="variable">
      <type>int</type>
      <name>len</name>
      <anchorfile>structpu__out.html</anchorfile>
      <anchor>a314b8f43c3ee0bf6060afbeced5dbe6c</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>int</type>
      <name>u_vals</name>
      <anchorfile>structpu__out.html</anchorfile>
      <anchor>a7697bc7a46cd1b8e37e337e708cb6023</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>int</type>
      <name>contribs</name>
      <anchorfile>structpu__out.html</anchorfile>
      <anchor>a638b0de1837cfd441871d005d3ab2938</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>char **</type>
      <name>header</name>
      <anchorfile>structpu__out.html</anchorfile>
      <anchor>ac9e9e30b16e7d04c770460b8487fb09d</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>double **</type>
      <name>u_values</name>
      <anchorfile>structpu__out.html</anchorfile>
      <anchor>a366edbc4170d5c177908e178ff340828</anchor>
      <arglist></arglist>
    </member>
  </compound>
  <compound kind="struct">
    <name>sect</name>
    <filename>structsect.html</filename>
  </compound>
  <compound kind="struct">
    <name>snoopT</name>
    <filename>structsnoopT.html</filename>
  </compound>
  <compound kind="struct">
    <name>SOLUTION</name>
    <filename>structSOLUTION.html</filename>
    <member kind="variable">
      <type>float</type>
      <name>energy</name>
      <anchorfile>structSOLUTION.html</anchorfile>
      <anchor>a4fe8e9027171f2dc4031587d7fab6b87</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>char *</type>
      <name>structure</name>
      <anchorfile>structSOLUTION.html</anchorfile>
      <anchor>a89ae453dfad0509468c39a62c303a63b</anchor>
      <arglist></arglist>
    </member>
  </compound>
  <compound kind="struct">
    <name>struct_en</name>
    <filename>structstruct__en.html</filename>
  </compound>
  <compound kind="struct">
    <name>svm_model</name>
    <filename>structsvm__model.html</filename>
  </compound>
  <compound kind="struct">
    <name>swString</name>
    <filename>structswString.html</filename>
  </compound>
  <compound kind="struct">
    <name>Tree</name>
    <filename>structTree.html</filename>
  </compound>
  <compound kind="struct">
    <name>TwoDfold_solution</name>
    <filename>structTwoDfold__solution.html</filename>
    <member kind="variable">
      <type>int</type>
      <name>k</name>
      <anchorfile>structTwoDfold__solution.html</anchorfile>
      <anchor>a298767110e07490d361bf7da920fd153</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>int</type>
      <name>l</name>
      <anchorfile>structTwoDfold__solution.html</anchorfile>
      <anchor>a64fb28259cf925c3bba7b8d14592363a</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>float</type>
      <name>en</name>
      <anchorfile>structTwoDfold__solution.html</anchorfile>
      <anchor>a3f65891d0c931f88440150bb32bcf753</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>char *</type>
      <name>s</name>
      <anchorfile>structTwoDfold__solution.html</anchorfile>
      <anchor>ac87e00bbdb13e0b6ef45c4f65608b416</anchor>
      <arglist></arglist>
    </member>
  </compound>
  <compound kind="struct">
    <name>TwoDfold_vars</name>
    <filename>structTwoDfold__vars.html</filename>
    <member kind="variable">
      <type>paramT *</type>
      <name>P</name>
      <anchorfile>structTwoDfold__vars.html</anchorfile>
      <anchor>ada74adef5f24b4b35c0b25da8223fe26</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>int</type>
      <name>do_backtrack</name>
      <anchorfile>structTwoDfold__vars.html</anchorfile>
      <anchor>ade5c7e9337a458ae20bac75abdc52d64</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>char *</type>
      <name>ptype</name>
      <anchorfile>structTwoDfold__vars.html</anchorfile>
      <anchor>aedf60b8b26dae05ad266d3e098d18208</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>char *</type>
      <name>sequence</name>
      <anchorfile>structTwoDfold__vars.html</anchorfile>
      <anchor>a3596f3d4d320318c4b8428e2abc7ab56</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>short *</type>
      <name>S1</name>
      <anchorfile>structTwoDfold__vars.html</anchorfile>
      <anchor>ab9ee459ffbfb5d2c138a033516056cdc</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>unsigned int</type>
      <name>maxD1</name>
      <anchorfile>structTwoDfold__vars.html</anchorfile>
      <anchor>a621ed2ab02116f3f8f5e7120dec429eb</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>unsigned int</type>
      <name>maxD2</name>
      <anchorfile>structTwoDfold__vars.html</anchorfile>
      <anchor>a03f198a4abdb3b784486d2ba5c533aa4</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>unsigned int *</type>
      <name>mm1</name>
      <anchorfile>structTwoDfold__vars.html</anchorfile>
      <anchor>aa11f5bcd8c4fe70a91c155c877c855d5</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>unsigned int *</type>
      <name>mm2</name>
      <anchorfile>structTwoDfold__vars.html</anchorfile>
      <anchor>a2eaa93316b6beb17531f0c078806036c</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>int *</type>
      <name>my_iindx</name>
      <anchorfile>structTwoDfold__vars.html</anchorfile>
      <anchor>a1a20cb06b58b75d1a3dbdbc8bc60d0a7</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>unsigned int *</type>
      <name>referenceBPs1</name>
      <anchorfile>structTwoDfold__vars.html</anchorfile>
      <anchor>a536525b98c1b633d4c5f2da4f8d78c18</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>unsigned int *</type>
      <name>referenceBPs2</name>
      <anchorfile>structTwoDfold__vars.html</anchorfile>
      <anchor>aa7abf73c3114cb5f0dc90e702fa9dd0f</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>unsigned int *</type>
      <name>bpdist</name>
      <anchorfile>structTwoDfold__vars.html</anchorfile>
      <anchor>af1106e1a592e2dccc92b3452340549e0</anchor>
      <arglist></arglist>
    </member>
  </compound>
  <compound kind="struct">
    <name>TwoDpfold_solution</name>
    <filename>structTwoDpfold__solution.html</filename>
    <member kind="variable">
      <type>int</type>
      <name>k</name>
      <anchorfile>structTwoDpfold__solution.html</anchorfile>
      <anchor>a40ad24e311b193866111623dd1331567</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>int</type>
      <name>l</name>
      <anchorfile>structTwoDpfold__solution.html</anchorfile>
      <anchor>aeaad6adc35413c76a2e2f18d96a6508c</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>double</type>
      <name>q</name>
      <anchorfile>structTwoDpfold__solution.html</anchorfile>
      <anchor>af0bf3071502b4a4fa81eeb6dfacef94c</anchor>
      <arglist></arglist>
    </member>
  </compound>
  <compound kind="struct">
    <name>TwoDpfold_vars</name>
    <filename>structTwoDpfold__vars.html</filename>
    <member kind="variable">
      <type>char *</type>
      <name>ptype</name>
      <anchorfile>structTwoDpfold__vars.html</anchorfile>
      <anchor>a67f37b8901b8d0a049c216d4c6241b07</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>char *</type>
      <name>sequence</name>
      <anchorfile>structTwoDpfold__vars.html</anchorfile>
      <anchor>a32c15a1e31856588259556c9020f32c6</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>short *</type>
      <name>S1</name>
      <anchorfile>structTwoDpfold__vars.html</anchorfile>
      <anchor>a240311ae1e8e121441651d6101e187ac</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>unsigned int</type>
      <name>maxD1</name>
      <anchorfile>structTwoDpfold__vars.html</anchorfile>
      <anchor>a7292b6cbc1ee5bacf55e842f316c4bef</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>unsigned int</type>
      <name>maxD2</name>
      <anchorfile>structTwoDpfold__vars.html</anchorfile>
      <anchor>a8900622d91454d2d037242e290e42834</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>int *</type>
      <name>my_iindx</name>
      <anchorfile>structTwoDpfold__vars.html</anchorfile>
      <anchor>ac2d3e6abf0cb0e1df363904fc938076e</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>int *</type>
      <name>jindx</name>
      <anchorfile>structTwoDpfold__vars.html</anchorfile>
      <anchor>a0699e194a797532c91b284ab10272384</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>unsigned int *</type>
      <name>referenceBPs1</name>
      <anchorfile>structTwoDpfold__vars.html</anchorfile>
      <anchor>aea15706d27b6b0fc19f5773919f43a8a</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>unsigned int *</type>
      <name>referenceBPs2</name>
      <anchorfile>structTwoDpfold__vars.html</anchorfile>
      <anchor>a1221396d712bf76b7f35297f2ab35a9f</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>unsigned int *</type>
      <name>bpdist</name>
      <anchorfile>structTwoDpfold__vars.html</anchorfile>
      <anchor>accef8eaa05fa57ca33aa22cbc7b7aaff</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>unsigned int *</type>
      <name>mm1</name>
      <anchorfile>structTwoDpfold__vars.html</anchorfile>
      <anchor>a7c9e9af6224d4696118e05835441863d</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>unsigned int *</type>
      <name>mm2</name>
      <anchorfile>structTwoDpfold__vars.html</anchorfile>
      <anchor>affb913470783f9edb12a0bfc22466269</anchor>
      <arglist></arglist>
    </member>
  </compound>
  <compound kind="dir">
    <name>/home/mescalin/ronny/WORK/ViennaRNA/H</name>
    <path>/home/mescalin/ronny/WORK/ViennaRNA/H/</path>
    <filename>dir_d72344b28b4f2089ce25682c4e6eba22.html</filename>
    <file>2Dfold.h</file>
    <file>2Dpfold.h</file>
    <file>ali_plex.h</file>
    <file>alifold.h</file>
    <file>aln_util.h</file>
    <file>cofold.h</file>
    <file>convert_epars.h</file>
    <file>data_structures.h</file>
    <file>dist_vars.h</file>
    <file>duplex.h</file>
    <file>edit_cost.h</file>
    <file>energy_const.h</file>
    <file>energy_par.h</file>
    <file>findpath.h</file>
    <file>fold.h</file>
    <file>fold_vars.h</file>
    <file>gquad.h</file>
    <file>inverse.h</file>
    <file>Lfold.h</file>
    <file>loop_energies.h</file>
    <file>LPfold.h</file>
    <file>MEA.h</file>
    <file>mm.h</file>
    <file>move_set.h</file>
    <file>naview.h</file>
    <file>pair_mat.h</file>
    <file>params.h</file>
    <file>part_func.h</file>
    <file>part_func_co.h</file>
    <file>part_func_up.h</file>
    <file>PKplex.h</file>
    <file>plex.h</file>
    <file>plot_layouts.h</file>
    <file>ProfileAln.h</file>
    <file>profiledist.h</file>
    <file>PS_dot.h</file>
    <file>read_epars.h</file>
    <file>ribo.h</file>
    <file>RNAstruct.h</file>
    <file>snofold.h</file>
    <file>snoop.h</file>
    <file>stringdist.h</file>
    <file>subopt.h</file>
    <file>svm_utils.h</file>
    <file>treedist.h</file>
    <file>utils.h</file>
  </compound>
  <compound kind="dir">
    <name>/home/mescalin/ronny/WORK/ViennaRNA/lib</name>
    <path>/home/mescalin/ronny/WORK/ViennaRNA/lib/</path>
    <filename>dir_97aefd0d527b934f1d99a682da8fe6a9.html</filename>
    <file>1.8.4_epars.h</file>
    <file>1.8.4_intloops.h</file>
    <file>intl11.h</file>
    <file>intl11dH.h</file>
    <file>intl21.h</file>
    <file>intl21dH.h</file>
    <file>intl22.h</file>
    <file>intl22dH.h</file>
    <file>list.h</file>
  </compound>
  <compound kind="page">
    <name>index</name>
    <title>ViennaRNA Package core - RNAlib</title>
    <filename>index</filename>
    <docanchor file="index" title="Introduction">mp_intro</docanchor>
  </compound>
</tagfile>
