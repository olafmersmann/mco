#ifndef MCO_EXTERN_H
#define MCO_EXTERN_H
#include <R.h>
#include <Rmath.h>
#include <Rinternals.h>

SEXP do_hv(SEXP s_data, SEXP s_ref);
SEXP do_eps_ind(SEXP s_data, SEXP s_ref);
SEXP do_nsga2(SEXP s_function,
              SEXP s_constraint,
              SEXP s_env,
              SEXP s_obj_dim,
              SEXP s_constr_dim,
              SEXP s_input_dim,
              SEXP s_lower_bound,
              SEXP s_upper_bound,
              SEXP s_popsize,
              SEXP s_generations,
              SEXP s_crossing_prob,
              SEXP s_crossing_dist,
              SEXP s_mutation_prob,
              SEXP s_mutation_dist);

#endif
