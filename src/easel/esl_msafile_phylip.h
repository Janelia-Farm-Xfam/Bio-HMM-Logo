/* I/O of multiple sequence alignments in PHYLIP format
 */
#ifndef eslMSAFILE_PHYLIP_INCLUDED
#define eslMSAFILE_PHYLIP_INCLUDED

#include "esl_msa.h"
#include "esl_msafile.h"

extern int esl_msafile_phylip_SetInmap     (ESLX_MSAFILE *afp);
extern int esl_msafile_phylip_GuessAlphabet(ESLX_MSAFILE *afp, int *ret_type);
extern int esl_msafile_phylip_Read         (ESLX_MSAFILE *afp, ESL_MSA **ret_msa);
extern int esl_msafile_phylip_Write        (FILE *fp, const ESL_MSA *msa, int format, ESLX_MSAFILE_FMTDATA *opt_fmtd);

extern int esl_msafile_phylip_CheckFileFormat(ESL_BUFFER *bf, int *ret_format, int *ret_namewidth);

#endif /* eslMSAFILE_PHYLIP_INCLUDED */
/*****************************************************************
 * @LICENSE@
 *
 * SVN $Id$
 * SVN $URL$
 *****************************************************************/
