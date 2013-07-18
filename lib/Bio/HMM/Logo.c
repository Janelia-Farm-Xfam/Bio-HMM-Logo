#include "hmmer.h"


#define perl_obj(pointer,class) ({                 \
  SV* ref=newSViv(0); SV* obj=newSVrv(ref, class); \
  sv_setiv(obj, (IV) pointer); SvREADONLY_on(obj); \
  ref;                                             \
})
#define c_obj(sv,type) (                           \
  (sv_isobject(sv) && sv_derived_from(sv, #type))  \
    ? ((type*)SvIV(SvRV(sv)))                      \
    : NULL                                         \
  )


AV*
convert_c2d_to_perl2d (float **vals, int row_start, int row_end, int col_start, int col_end) {
   AV *arr2d = newAV();
   AV *row;
   int i,j;

   for ( i = row_start; i <= row_end; ++i ) {
      row = newAV();
      for ( j = col_start; j <= col_end; ++j ) {
        av_push( row, newSVnv(vals[i][j]) );
      }
      av_push( arr2d, newRV_noinc( (SV*) row ) );
   }
   return newRV_noinc( arr2d );
}


SV*
inline_get_MM_array (SV *hmm_p) {
  int           status;
  int           i;
  P7_HMM       *hmm       = c_obj(hmm_p, P7_HMM);
  AV           *mm        = newAV();


  for ( i = 1; i <= hmm->M; ++i ) 
     av_push( mm, newSVnv( (hmm->mm != NULL && hmm->mm[i] == 'm') ? 1 : 0  )   ) ;
  
  
  return newRV_noinc( mm );

ERROR:

  croak("error getting mm line ");
  return NULL;

}


SV*
inline_get_emission_heights (SV *hmm_p) {
  int           status;
  int           i, j;
  float        *rel_ents  = NULL; //passed into hmmlogo_emissionHeightsDivRelent(), but not used
  float        **heights  = NULL;
  P7_HMM       *hmm       = c_obj(hmm_p, P7_HMM);
  ESL_ALPHABET *abc       = hmm->abc;
  P7_BG        *bg        = p7_bg_Create(abc);
  SV           *arr2d_sv;


  ESL_ALLOC(rel_ents, (hmm->M+1) * sizeof(float));
  ESL_ALLOC(heights,  (hmm->M+1) * sizeof(float*));
  for (i = 1; i <= hmm->M; i++)
    ESL_ALLOC(heights[i], abc->K * sizeof(float));

  hmmlogo_emissionHeightsDivRelent(hmm, bg, rel_ents, heights);

  arr2d_sv = convert_c2d_to_perl2d(heights, 1, hmm->M, 0, abc->K -1);

  p7_bg_Destroy(bg);
  free(rel_ents);
  for (i = 1; i <= hmm->M; i++)
    free(heights[i]);
  free(heights);

  return (arr2d_sv);

ERROR:
  if (bg != NULL) p7_bg_Destroy(bg);
  if (rel_ents != NULL) free(rel_ents);
  for (i = 1; i <= hmm->M; i++)
    if (heights[i] != NULL) free(heights[i]);
  if (heights != NULL)  free(heights);

  croak("error getting emission heights");
  return NULL;

}


SV*
inline_get_posscore_heights (SV *hmm_p) {
  int           status;
  int           i, j;
  float        *rel_ents  = NULL; //passed into hmmlogo_emissionHeightsDivRelent(), but not used
  float        **heights  = NULL;
  P7_HMM       *hmm       = c_obj(hmm_p, P7_HMM);
  ESL_ALPHABET *abc       = hmm->abc;
  P7_BG        *bg        = p7_bg_Create(abc);
  SV           *arr2d_sv;


  ESL_ALLOC(rel_ents, (hmm->M+1) * sizeof(float));
  ESL_ALLOC(heights,  (hmm->M+1) * sizeof(float*));
  for (i = 1; i <= hmm->M; i++)
    ESL_ALLOC(heights[i], abc->K * sizeof(float));

  hmmlogo_posScoreHeightsDivRelent(hmm, bg, rel_ents, heights);
  arr2d_sv = convert_c2d_to_perl2d(heights, 1, hmm->M, 0, abc->K -1);

  p7_bg_Destroy(bg);
  free(rel_ents);
  for (i = 1; i <= hmm->M; i++)
    free(heights[i]);
  free(heights);

  return (arr2d_sv);

ERROR:
  if (bg != NULL) p7_bg_Destroy(bg);
  if (rel_ents != NULL) free(rel_ents);
  for (i = 1; i <= hmm->M; i++)
    if (heights[i] != NULL) free(heights[i]);
  if (heights != NULL)  free(heights);

  croak("error getting pos_score heights");
  return NULL;

}

SV*
inline_get_score_heights (SV *hmm_p) {
  int           status;
  int           i, j;
  float        **heights  = NULL;
  P7_HMM       *hmm       = c_obj(hmm_p, P7_HMM);
  ESL_ALPHABET *abc       = hmm->abc;
  P7_BG        *bg        = p7_bg_Create(abc);
  SV           *arr2d_sv;


  ESL_ALLOC(heights,  (hmm->M+1) * sizeof(float*));
  for (i = 1; i <= hmm->M; i++)
    ESL_ALLOC(heights[i], abc->K * sizeof(float));

  hmmlogo_ScoreHeights (hmm, bg, heights );
  arr2d_sv = convert_c2d_to_perl2d(heights, 1, hmm->M, 0, abc->K -1);

  p7_bg_Destroy(bg);
  for (i = 1; i <= hmm->M; i++)
    free(heights[i]);
  free(heights);

  return (arr2d_sv);

ERROR:
  if (bg != NULL) p7_bg_Destroy(bg);
  for (i = 1; i <= hmm->M; i++)
    if (heights[i] != NULL) free(heights[i]);
  if (heights != NULL)  free(heights);

  croak("error getting score heights");
  return NULL;

}


SV*
inline_get_insertP (SV *hmm_p) {
  int           status;
  int           i;
  P7_HMM       *hmm       = c_obj(hmm_p, P7_HMM);
  float        *ins_P     = NULL;
  AV           *p_vals    = newAV();

  ESL_ALLOC(ins_P, (hmm->M+1) * sizeof(float));
  hmmlogo_IndelValues(hmm, ins_P, NULL, NULL);

   for ( i = 1; i <= hmm->M; ++i )
      av_push( p_vals, newSVnv(ins_P[i]) );

  free(ins_P);

  return newRV_noinc( p_vals );

ERROR:
  if (ins_P!= NULL) free(ins_P);

  croak("error getting insert P ");
  return NULL;

}

SV*
inline_get_insertLengths (SV *hmm_p) {
  int           status;
  int           i;
  P7_HMM       *hmm        = c_obj(hmm_p, P7_HMM);
  float        *ins_lens   = NULL;
  AV           *lengths    = newAV();

  ESL_ALLOC(ins_lens, (hmm->M+1) * sizeof(float));
  hmmlogo_IndelValues(hmm, NULL, ins_lens, NULL);

   for ( i = 1; i <= hmm->M; ++i )
      av_push( lengths, newSVnv(ins_lens[i]) );

  free(ins_lens);

  return newRV_noinc( lengths );

ERROR:
  if (ins_lens!= NULL) free(ins_lens);

  croak("error getting insert lengths ");
  return NULL;

}

SV*
inline_get_deleteP (SV *hmm_p) {
  int           status;
  int           i;
  P7_HMM       *hmm       = c_obj(hmm_p, P7_HMM);
  float        *del_P     = NULL;
  AV           *p_vals    = newAV();

  ESL_ALLOC(del_P, (hmm->M+1) * sizeof(float));
  hmmlogo_IndelValues(hmm, NULL, NULL, del_P);

   for ( i = 1; i <= hmm->M; ++i )
      av_push( p_vals, newSVnv(del_P[i]) );

  free(del_P);

  return newRV_noinc( p_vals );

ERROR:
  if (del_P!= NULL) free(del_P);

  croak("error getting delete P ");
  return NULL;

}



SV*
inline_read_hmm (char *filename) {
  P7_HMMFILE      *hfp      = NULL;              // open input HMM file
  P7_HMM          *hmm      = NULL;              // one HMM query
  ESL_ALPHABET    *abc      = NULL;              // digital alphabet

  char   errbuf[eslERRBUFSIZE];
  int status;


  // Open the query profile HMM file
  status = p7_hmmfile_OpenE(filename, NULL, &hfp, errbuf);
  if (status != eslOK)  croak("unable to open HMM file");

  // read the file
  status = p7_hmmfile_Read(hfp, &abc, &hmm);
  if (status != eslOK)  croak("error reading HMM file");

  p7_hmmfile_Close(hfp);

  return perl_obj(hmm,"P7_HMM");
}



float
inline_hmmlogo_maxHeight( SV* abc_p) {
    ESL_ALPHABET *abc = c_obj(abc_p, ESL_ALPHABET);
    P7_BG        *bg  = p7_bg_Create(abc);
    return hmmlogo_maxHeight(bg);
}


SV*
inline_get_abc ( SV *hmm_p) {
  P7_HMM      *hmm      = c_obj(hmm_p, P7_HMM);
  return perl_obj(hmm->abc,"ESL_ALPHABET");
}


int inline_get_abc_type ( SV *hmm_p) {
  P7_HMM      *hmm      = c_obj(hmm_p, P7_HMM);
  return hmm->abc->type;
}

SV*
inline_destroy_abc ( SV *abc_p) {
  ESL_ALPHABET  *abc      = c_obj(abc_p, ESL_ALPHABET);
  esl_alphabet_Destroy(abc);
}


SV*
inline_destroy_hmm ( SV *hmm_p) {
  P7_HMM  *hmm      = c_obj(hmm_p, P7_HMM);
  p7_hmm_Destroy(hmm);
}


SV*
inline_get_alphabet_string (SV *abc_p) {
  ESL_ALPHABET   *abc   = c_obj(abc_p, ESL_ALPHABET);
  int i;
  return newSVpv(abc->sym, abc->K);
}



