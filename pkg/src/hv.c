/**
 * hv.c - Calculate dominated hypervolume
 *
 * Authors:
 *  Heike Trautmann  <trautmann@statistik.uni-dortmund.de>
 *  Detlef Steuer    <detlef.steuer@hsu-hamburg.de>
 *  Olaf Mersmann    <olafm@statistik.uni-dortmund.de>
 *
 * Based on C code written by
 *
 *                       Copyright (c) 2005, 2006
 *                 Carlos Fonseca <cmfonsec@ualg.pt>
 *            Manuel Lopez-Ibanez <m.lopez-ibanez@napier.ac.uk>
 *                   Luis Paquete <lpaquete@ualg.pt>
 * 
 * and AVL tree code by
 *
 *   Copyright (C) 1998  Michael H. Buselli <cosine@cosine.org>
 *   Copyright (C) 2000-2002  Wessel Dankers <wsl@nl.linux.org>
 *
 * LICENSE:
 * 
 * This program is free software (software libre); you can redistribute
 * it and/or modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2 of the
 * License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, you can obtain a copy of the GNU
 * General Public License at:
 *                 http://www.gnu.org/copyleft/gpl.html
 * or by writing to:
 *           Free Software Foundation, Inc., 59 Temple Place,
 *                 Suite 330, Boston, MA 02111-1307 USA
 *
 *
 * Relevant literature:
 *
 * [1]  C. M. Fonseca, L. Paquete, and M. Lopez-Ibanez. An
 *      improved dimension-sweep algorithm for the hypervolume
 *      indicator. In IEEE Congress on Evolutionary Computation,
 *      pages 1157-1163, Vancouver, Canada, July 2006.
 *
 * [2]  L. Paquete, C. M. Fonseca and M. Lopez-Ibanez. An optimal
 *      algorithm for a special case of Klee's measure problem in three
 *      dimensions. Technical Report CSI-RT-I-01/2006, CSI, Universidade
 *      do Algarve, 2006.
 *
 * CHANGES:
 *  * Merged avl.[ch] and hv.[ch] (OME)
 *  * Removed VARIANT stuff, only variant 4 used. (OME)
 *  * Use R memory functions
 */

#include <stdlib.h>
#include <stdio.h>
#include <limits.h>
#include <float.h>

#include <R.h>
#include <Rmath.h>
#include <Rinternals.h>
#include <R_ext/Applic.h>

/**
 * avl.c - AVL tree
 */
#include <errno.h>
#define AVL_DEPTH
#define AVL_COUNT

typedef int (*avl_compare_t)(const void *, const void *);
typedef void (*avl_freeitem_t)(void *);
typedef struct avl_node_t {
  struct avl_node_t *next;
  struct avl_node_t *prev;
  struct avl_node_t *parent;
  struct avl_node_t *left;
  struct avl_node_t *right;
  void *item;
  unsigned int count;
  unsigned char depth;
} avl_node_t;

typedef struct avl_tree_t {
  avl_node_t *head;
  avl_node_t *tail;
  avl_node_t *top;
  avl_compare_t cmp;
  avl_freeitem_t freeitem;
} avl_tree_t;

static void avl_rebalance(avl_tree_t *, avl_node_t *);

#define NODE_COUNT(n)  ((n) ? (n)->count : 0)
#define L_COUNT(n)     (NODE_COUNT((n)->left))
#define R_COUNT(n)     (NODE_COUNT((n)->right))
#define CALC_COUNT(n)  (L_COUNT(n) + R_COUNT(n) + 1)
#define NODE_DEPTH(n)  ((n) ? (n)->depth : 0)
#define L_DEPTH(n)     (NODE_DEPTH((n)->left))
#define R_DEPTH(n)     (NODE_DEPTH((n)->right))
#define CALC_DEPTH(n)  ((L_DEPTH(n)>R_DEPTH(n)?L_DEPTH(n):R_DEPTH(n)) + 1)

static int avl_check_balance(avl_node_t *avlnode) {
  int d;
  d = R_DEPTH(avlnode) - L_DEPTH(avlnode);
  return d<-1?-1:d>1?1:0;
}

static int avl_search_closest(const avl_tree_t *avltree, const void *item, avl_node_t **avlnode) {
  avl_node_t *node;
  avl_compare_t cmp;
  int c;

  if(!avlnode)
    avlnode = &node;

  node = avltree->top;

  if(!node)
    return *avlnode = NULL, 0;

  cmp = avltree->cmp;

  for(;;) {
    c = cmp(item, node->item);

    if(c < 0) {
      if(node->left)
	node = node->left;
      else
	return *avlnode = node, -1;
    } else if(c > 0) {
      if(node->right)
	node = node->right;
      else
	return *avlnode = node, 1;
    } else {
      return *avlnode = node, 0;
    }
  }
}

static avl_tree_t *avl_init_tree(avl_tree_t *rc, avl_compare_t cmp, avl_freeitem_t freeitem) {
  if(rc) {
    rc->head = NULL;
    rc->tail = NULL;
    rc->top = NULL;
    rc->cmp = cmp;
    rc->freeitem = freeitem;
  }
  return rc;
}

static avl_tree_t *avl_alloc_tree(avl_compare_t cmp, avl_freeitem_t freeitem) {
    return avl_init_tree((avl_tree_t *)R_alloc(1, sizeof(avl_tree_t)), 
			 cmp, 
			 freeitem);
}

static void avl_clear_tree(avl_tree_t *avltree) {
  avltree->top = avltree->head = avltree->tail = NULL;
}

static void avl_free_nodes(avl_tree_t *avltree) {
  avl_node_t *node, *next;
  avl_freeitem_t freeitem;
  
  freeitem = avltree->freeitem;
  
  for(node = avltree->head; node; node = next) {
    next = node->next;
    if(freeitem)
      freeitem(node->item);
    //free(node);
  }
  avl_clear_tree(avltree);
}

static void avl_free_tree(avl_tree_t *avltree) {
  avl_free_nodes(avltree);
  //free(avltree);
}

static void avl_clear_node(avl_node_t *newnode) {
  newnode->left = newnode->right = NULL;
  newnode->count = 1;
  newnode->depth = 1;
}

static avl_node_t *avl_init_node(avl_node_t *newnode, void *item) {
  if(newnode) {
    avl_clear_node(newnode); 
    newnode->item = item;
  }
  return newnode;
}

static avl_node_t *avl_insert_top(avl_tree_t *avltree, avl_node_t *newnode) {
  avl_clear_node(newnode);
  newnode->prev = newnode->next = newnode->parent = NULL;
  avltree->head = avltree->tail = avltree->top = newnode;
  return newnode;
}

static avl_node_t *avl_insert_after(avl_tree_t *avltree, avl_node_t *node, avl_node_t *newnode);
static avl_node_t *avl_insert_before(avl_tree_t *avltree, avl_node_t *node, avl_node_t *newnode) {
  if(!node) {
    if (NULL == avltree->tail)
      return avl_insert_after(avltree, avltree->tail, newnode);
    else
      return avl_insert_top(avltree, newnode);
  }
  if(node->left)
    return avl_insert_after(avltree, node->prev, newnode);

  avl_clear_node(newnode);

  newnode->next = node;
  newnode->parent = node;

  newnode->prev = node->prev;
  if(node->prev)
    node->prev->next = newnode;
  else
    avltree->head = newnode;
  node->prev = newnode;

  node->left = newnode;
  avl_rebalance(avltree, node);
  return newnode;
}

static avl_node_t *avl_insert_after(avl_tree_t *avltree, avl_node_t *node, avl_node_t *newnode) {
  if(!node)
    return avltree->head
      ? avl_insert_before(avltree, avltree->head, newnode)
      : avl_insert_top(avltree, newnode);

  if(node->right)
    return avl_insert_before(avltree, node->next, newnode);

  avl_clear_node(newnode);

  newnode->prev = node;
  newnode->parent = node;

  newnode->next = node->next;
  if(node->next)
    node->next->prev = newnode;
  else
    avltree->tail = newnode;
  node->next = newnode;

  node->right = newnode;
  avl_rebalance(avltree, node);
  return newnode;
}

static void avl_unlink_node(avl_tree_t *avltree, avl_node_t *avlnode) {
  avl_node_t *parent;
  avl_node_t **superparent;
  avl_node_t *subst, *left, *right;
  avl_node_t *balnode;

  if(avlnode->prev)
    avlnode->prev->next = avlnode->next;
  else
    avltree->head = avlnode->next;

  if(avlnode->next)
    avlnode->next->prev = avlnode->prev;
  else
    avltree->tail = avlnode->prev;

  parent = avlnode->parent;

  superparent = parent
    ? avlnode == parent->left ? &parent->left : &parent->right
    : &avltree->top;

  left = avlnode->left;
  right = avlnode->right;
  if(!left) {
    *superparent = right;
    if(right)
      right->parent = parent;
    balnode = parent;
  } else if(!right) {
    *superparent = left;
    left->parent = parent;
    balnode = parent;
  } else {
    subst = avlnode->prev;
    if(subst == left) {
      balnode = subst;
    } else {
      balnode = subst->parent;
      balnode->right = subst->left;
      if(balnode->right)
	balnode->right->parent = balnode;
      subst->left = left;
      left->parent = subst;
    }
    subst->right = right;
    subst->parent = parent;
    right->parent = subst;
    *superparent = subst;
  }

  avl_rebalance(avltree, balnode);
}

static void avl_rebalance(avl_tree_t *avltree, avl_node_t *avlnode) {
  avl_node_t *child;
  avl_node_t *gchild;
  avl_node_t *parent;
  avl_node_t **superparent;

  parent = avlnode;

  while(avlnode) {
    parent = avlnode->parent;

    superparent = parent
      ? avlnode == parent->left ? &parent->left : &parent->right
      : &avltree->top;

    switch(avl_check_balance(avlnode)) {
    case -1:
      child = avlnode->left;
      if(L_DEPTH(child) >= R_DEPTH(child)) {
	avlnode->left = child->right;
	if(avlnode->left)
	  avlnode->left->parent = avlnode;
	child->right = avlnode;
	avlnode->parent = child;
	*superparent = child;
	child->parent = parent;
	avlnode->count = CALC_COUNT(avlnode);
	child->count = CALC_COUNT(child);
	avlnode->depth = CALC_DEPTH(avlnode);
	child->depth = CALC_DEPTH(child);
      } else {
	gchild = child->right;
	avlnode->left = gchild->right;
	if(avlnode->left)
	  avlnode->left->parent = avlnode;
	child->right = gchild->left;
	if(child->right)
	  child->right->parent = child;
	gchild->right = avlnode;
	if(gchild->right)
	  gchild->right->parent = gchild;
	gchild->left = child;
	if(gchild->left)
	  gchild->left->parent = gchild;
	*superparent = gchild;
	gchild->parent = parent;
	avlnode->count = CALC_COUNT(avlnode);
	child->count = CALC_COUNT(child);
	gchild->count = CALC_COUNT(gchild);
	avlnode->depth = CALC_DEPTH(avlnode);
	child->depth = CALC_DEPTH(child);
	gchild->depth = CALC_DEPTH(gchild);
      }
      break;
    case 1:
      child = avlnode->right;
      if(R_DEPTH(child) >= L_DEPTH(child)) {
	avlnode->right = child->left;
	if(avlnode->right)
	  avlnode->right->parent = avlnode;
	child->left = avlnode;
	avlnode->parent = child;
	*superparent = child;
	child->parent = parent;
	avlnode->count = CALC_COUNT(avlnode);
	child->count = CALC_COUNT(child);
	avlnode->depth = CALC_DEPTH(avlnode);
	child->depth = CALC_DEPTH(child);
      } else {
	gchild = child->left;
	avlnode->right = gchild->left;
	if(avlnode->right)
	  avlnode->right->parent = avlnode;
	child->left = gchild->right;
	if(child->left)
	  child->left->parent = child;
	gchild->left = avlnode;
	if(gchild->left)
	  gchild->left->parent = gchild;
	gchild->right = child;
	if(gchild->right)
	  gchild->right->parent = gchild;
	*superparent = gchild;
	gchild->parent = parent;
	avlnode->count = CALC_COUNT(avlnode);
	child->count = CALC_COUNT(child);
	gchild->count = CALC_COUNT(gchild);
	avlnode->depth = CALC_DEPTH(avlnode);
	child->depth = CALC_DEPTH(child);
	gchild->depth = CALC_DEPTH(gchild);
      }
      break;
    default:
      avlnode->count = CALC_COUNT(avlnode);
      avlnode->depth = CALC_DEPTH(avlnode);
    }
    avlnode = parent;
  }
}
/** END AVL tree */

typedef struct dlnode {
  double *x;                    /* The data vector              */
  struct dlnode **next;         /* Next-node vector             */
  struct dlnode **prev;         /* Previous-node vector         */
  struct avl_node_t * tnode;
  int ignore;
  double *area;                 /* Area */
  double *vol;                  /* Volume */
} dlnode_t;


static avl_tree_t *tree;
static const int stop_dimension = 2; /* Use special case for last two dimensions. */

static int compare_node( const void *p1, const void* p2) {
  const double x1 = *((*(const dlnode_t **)p1)->x);
  const double x2 = *((*(const dlnode_t **)p2)->x);

  if (x1 == x2)
    return 0;
  else
    return ( x1 < x2 ) ? -1 : 1;
}

static int compare_tree_asc( const void *p1, const void *p2) {
  const double x1 = *((const double *)p1+1);
  const double x2 = *((const double *)p2+1);

  if (x1 != x2)
    return (x1 > x2) ? -1 : 1;
  else
    return 0;
}

/*
 * Setup circular double-linked list in each dimension
 */
static dlnode_t * setup_cdllist(double *data, int d, int n) {
  dlnode_t *head;
  dlnode_t **scratch;
  int i, j;

  head = R_alloc (n+1, sizeof(dlnode_t));

  head->x = data;
  head->ignore = 0;  /* should never get used */
  head->next   = R_alloc(d * (n+1), sizeof(dlnode_t*));
  head->prev   = R_alloc(d * (n+1), sizeof(dlnode_t*));
  head->tnode  = (dlnode_t *)R_alloc(1, sizeof(avl_node_t));
  head->area   = R_alloc(d * (n+1), sizeof(double));
  head->vol    = R_alloc(d * (n+1), sizeof(double));

  for (i = 1; i <= n; i++) {
      head[i].x = head[i-1].x + d ;/* this will be fixed a few lines below... */
      head[i].ignore = 0;
      head[i].next = head[i-1].next + d;
      head[i].prev = head[i-1].prev + d;
      head[i].tnode = (dlnode_t *)R_alloc(1, sizeof(avl_node_t));
      head[i].area = head[i-1].area + d;
      head[i].vol = head[i-1].vol + d;
  }
  head->x = NULL; /* head contains no data */
  
  scratch = Calloc(n, dlnode_t*);
  for (i = 0; i < n; i++)
      scratch[i] = head + i + 1;
  
  for (j = d-1; j >= 0; j--) {
      for (i = 0; i < n; i++)
	  scratch[i]->x--;
      qsort(scratch, n, sizeof(dlnode_t*), compare_node);
      head->next[j] = scratch[0];
      scratch[0]->prev[j] = head;
      for (i = 1; i < n; i++) {
	  scratch[i-1]->next[j] = scratch[i];
	  scratch[i]->prev[j] = scratch[i-1];
      }
      scratch[n-1]->next[j] = head;
      head->prev[j] = scratch[n-1];
  }
  Free(scratch);
  return head;
}

static void delete(dlnode_t *nodep, int dim, double *bound) {
    int i;
    
    for (i = 0; i < dim; i++) {
	nodep->prev[i]->next[i] = nodep->next[i];
	nodep->next[i]->prev[i] = nodep->prev[i];
	if (bound[i] > nodep->x[i])
	    bound[i] = nodep->x[i];
    }
}

static void reinsert (dlnode_t *nodep, int dim, double * bound) {
    int i;
    
    for (i = 0; i < dim; i++) {
	nodep->prev[i]->next[i] = nodep;
	nodep->next[i]->prev[i] = nodep;
	if (bound[i] > nodep->x[i])
	    bound[i] = nodep->x[i];
    }
}

static double hv_recursive( dlnode_t *list, int dim, int c, const double * ref, double * bound) {
    dlnode_t *p0,*p1,*pp;
    double hypera,hyperv=0;
    double height;
    

    /* ------------------------------------------------------
       General case for dimensions higher than stop_dimension
       ------------------------------------------------------ */
    if ( dim > stop_dimension ) {
	p0 = list;
	for (p1 = p0->prev[dim]; p1->x; p1 = p1->prev[dim]) {
	    if (p1->ignore < dim)
		p1->ignore = 0;
	}

	p1 = p0->prev[dim];
	while (c > 1 && (p1->x[dim] > bound[dim] || p1->prev[dim]->x[dim] >= bound[dim])) {
	    p0 = p1;
	    delete(p0, dim, bound);
	    p1 = p0->prev[dim];
	    c--;
	}
	
	if (c > 1)
	    hyperv = p1->prev[dim]->vol[dim] + p1->prev[dim]->area[dim]
		* (p1->x[dim] - p1->prev[dim]->x[dim]);
	else {
	    p1->area[0] = 1;
	    int i;
	    for (i = 1; i <= dim; i++)
		p1->area[i] = p1->area[i-1] * (ref[i-1] - p1->x[i-1]);
	}
	p1->vol[dim] = hyperv;
	if (p1->ignore >= dim)
	    p1->area[dim] = p1->prev[dim]->area[dim];
	else {
	    p1->area[dim] = hv_recursive(list, dim-1, c, ref, bound);
	    if (p1->area[dim] <= p1->prev[dim]->area[dim])
		p1->ignore = dim;
	}
	
	while (p0->x != NULL) {
	    hyperv += p1->area[dim] * (p0->x[dim] - p1->x[dim]);
	    bound[dim] = p0->x[dim];
	    reinsert(p0, dim, bound);
	    c++;
	    p1 = p0;
	    p0 = p0->next[dim];
	    p1->vol[dim] = hyperv;
	    if (p1->ignore >= dim)
		p1->area[dim] = p1->prev[dim]->area[dim];
	    else {
		p1->area[dim] = hv_recursive(list, dim-1, c, ref, bound);
		if (p1->area[dim] <= p1->prev[dim]->area[dim])
		    p1->ignore = dim;
	    }
	}
	hyperv +=  p1->area[dim] * (ref[dim] - p1->x[dim]);
	return hyperv;
    } else if (dim == 2) { // 3D case:
	pp = list->next[2];
	avl_init_node(pp->tnode,pp->x);
	avl_insert_top(tree,pp->tnode);
	
	hypera = (ref[0] - pp->x[0])*(ref[1] - pp->x[1]);
	
	if (c == 1)
	    height = ref[2] - pp->x[2];
	else
	    height = pp->next[2]->x[2] - pp->x[2];
	
	hyperv = hypera * height;
	
	while ((pp=pp->next[2])->x) {
	    height = (pp==list->prev[2])
		? ref[2] - pp->x[2]
		: pp->next[2]->x[2] - pp->x[2];
	    if (pp->ignore>=2)
		hyperv += hypera * height;
	    else {
		const double * prv_ip, * nxt_ip;
		avl_node_t *tnode;
		
		avl_init_node(pp->tnode, pp->x);
		
		if (avl_search_closest(tree, pp->x, &tnode) <= 0) {
		    nxt_ip = (double *)(tnode->item);
		    tnode = tnode->prev;
		} else {
		    nxt_ip = (tnode->next!=NULL)
			? (double *)(tnode->next->item)
			: ref;
		}
		
		if (nxt_ip[0] > pp->x[0]) {		    
		    avl_insert_after(tree, tnode, pp->tnode);
		    if (tnode !=NULL) {
			prv_ip = (double *)(tnode->item);
			if (prv_ip[0] > pp->x[0]) {
			    const double * cur_ip;
			    
			    tnode = pp->tnode->prev;
			    // cur_ip = point dominated by pp with highest [0]-coordinate
			    cur_ip = (double *)(tnode->item);
			    while (tnode->prev) {
				prv_ip = (double *)(tnode->prev->item);
				hypera -= (prv_ip[1] - cur_ip[1])*(nxt_ip[0] - cur_ip[0]);
				if (prv_ip[0] < pp->x[0])
				    break; // prv is not dominated by pp
				cur_ip = prv_ip;
				avl_unlink_node(tree,tnode);
				tnode = tnode->prev;
			    }			    
			    avl_unlink_node(tree,tnode);			    
			    if (!tnode->prev) {
				hypera -= (ref[1] - cur_ip[1])*(nxt_ip[0] - cur_ip[0]);
				prv_ip = ref;
			    }
			}
		    } else
			prv_ip = ref;		    
		    hypera += (prv_ip[1] - pp->x[1])*(nxt_ip[0] - pp->x[0]);		    
		} else
		    pp->ignore = 2;		
		if (height > 0)
		    hyperv += hypera * height;
	    }
	}
	avl_clear_tree(tree);
	return hyperv;
    } else if (dim == 1) { /* special case of dimension 2 */
	p1 = list->next[1];
	hypera = p1->x[0];
	while ((p0 = p1->next[1])->x) {
	    hyperv += (ref[0] - hypera) * (p0->x[1] - p1->x[1]);
	    if (p0->x[0] < hypera)
		hypera = p0->x[0];
	    p1 = p0;
	}
	hyperv += (ref[0] - hypera) * (ref[1] - p1->x[1]);
	return hyperv;
    } else if (dim == 0) {     /* special case of dimension 1 */
	return (ref[0] - list->next[0]->x[0]);
    } else {
	fprintf(stderr, "%s:%d: unreachable condition! \n"
		"This is a bug, please report it to "
		"m.lopez-ibanez@napier.ac.uk\n", __FILE__, __LINE__);
	exit(EXIT_FAILURE);
    }
}

#define UNPACK_REAL_VECTOR(S, D, N)             \
  double *D = REAL(S);                          \
  const R_len_t N = length(S);

#define UNPACK_REAL_MATRIX(S, D, N, K)		\
  double *D = REAL(S);                          \
  const R_len_t N = nrows(S);			\
  const R_len_t K = ncols(S);


SEXP do_hv(SEXP s_data, SEXP s_ref) {
  SEXP s_res;
  dlnode_t *list;
  double *bound;
  int i;

  if (!isReal(s_data))  error("Argument 's_data' is not a real matrix.");
  if (!isReal(s_ref))   error("Argument 's_ref' is not a real vector.");

  /* Unpack arguments */
  UNPACK_REAL_MATRIX(s_data, data, k_data, n_data);
  UNPACK_REAL_VECTOR(s_ref, ref, n_ref);

  if (n_ref != k_data)
    error("ref and data must be same dimension.");

  /* Allocate result */
  PROTECT(s_res = allocVector(REALSXP, 1));
  double *res = REAL(s_res);

  bound = Calloc(k_data, double);
  for (i = 0; i < k_data; i++) 
    bound[i] = -DBL_MAX;

  tree = avl_alloc_tree((avl_compare_t) compare_tree_asc, (avl_freeitem_t) free); 
  list = setup_cdllist(data, k_data, n_data); 
  res[0] = hv_recursive(list, k_data - 1, n_data, ref, bound);
  Free(bound);
  avl_free_tree(tree);
  UNPROTECT(1); /* s_res */
  return s_res;
}
