/*
Copyright ESIEE (2009) 

m.couprie@esiee.fr

This software is an image processing library whose purpose is to be
used primarily for research and teaching.

This software is governed by the CeCILL  license under French law and
abiding by the rules of distribution of free software. You can  use, 
modify and/ or redistribute the software under the terms of the CeCILL
license as circulated by CEA, CNRS and INRIA at the following URL
"http://www.cecill.info". 

As a counterpart to the access to the source code and  rights to copy,
modify and redistribute granted by the license, users are provided only
with a limited warranty  and the software's author,  the holder of the
economic rights,  and the successive licensors  have only  limited
liability. 

In this respect, the user's attention is drawn to the risks associated
with loading,  using,  modifying and/or developing or reproducing the
software by the user in light of its specific status of free software,
that may mean  that it is complicated to manipulate,  and  that  also
therefore means  that it is reserved for developers  and  experienced
professionals having in-depth computer knowledge. Users are therefore
encouraged to load and test the software's suitability as regards their
requirements in conditions enabling the security of their systems and/or 
data to be ensured and,  more generally, to use and operate it in the 
same conditions as regards security. 

The fact that you are presently reading this means that you have had
knowledge of the CeCILL license and that you accept its terms.
*/

/* authors : J. Cousty - L. Najman and M. Couprie */


/* $Id: lminimgraph */
#define LABMAX   0
#define LABMIN   1


/* ============== */
/* prototype for lminimgraph.c */
/* ============== */

extern int32_t lminimgraph(
        graphe *g,
        int32_t minimum,         
        int32_t *LABELe,
        int32_t *LABELv, 
        int32_t *nlabels
);

extern int32_t lminimvertex(
        graphe *g,
        double *Fv
        );

extern int32_t lmborderwshed(
         graphe *g, /* graphe de depart */                
        double *Fv, /*  */
        int32_t *LABELe,
	 int32_t *LABELv//,
	 //	double *depth,
	 //int32_t *lowedge
	);

extern int32_t lmbordercut(
         graphe *g, /* graphe de depart */                
        double *Fv, /*  */
        int32_t *LABELe,
	int32_t *LABELv,
	int32_t *cut
	);

extern u_int32_t mBorderWshed2drapide(graphe *g, u_int32_t ** Res);

extern void bKernelToRootedMSF(graphe *g, u_int32_t *Lab);
