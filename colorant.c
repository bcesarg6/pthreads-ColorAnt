/*************************************************************************
 * Created: Qua 27 Jan 2011 15:17:39 BRST
 *
 * Author: Carla N. Lintzmayer, carla0negri@gmail.com
 *
 * Modificado em 2017
 * 
 * Autor: Bruno Cesar Puli Dala Rosa, bcesar.g6@gmail.com
 *
 *************************************************************************/

#include <stdio.h>
#include <string.h>
#include <math.h>
#include <limits.h>
#include <pthread.h>

#include "colorant.h"
#include "tabucol.h"
#include "ant_fixed_k.h"
#include "util.h"

/* Global data */
static double **pheromone;
static double ***phero_var; //new dimension = thread-level
static gcp_solution_t *ant_k;
static gcp_solution_t *best_colony;
static gcp_solution_t **best_colonies;
static gcp_solution_t *best_ant;

static aco_memory_t *memory = NULL;
static gcp_solution_t *ant_memory_remove;
static gcp_solution_t *ant_memory_insert;

void colorant_printbanner(void) {

  char *schemes[] = {"All ants + Best ant + Best colony", "Best ant + Best colony", "Best ant + Best colony (gap)"};

  fprintf(problem->fileout, "COLORANT\n");
  fprintf(problem->fileout, "-------------------------------------------------\n");
  fprintf(problem->fileout, "Parameters:\n");
  fprintf(problem->fileout, "  Alpha..................................: %.2f\n", aco_info->alpha);
  fprintf(problem->fileout, "  Beta...................................: %.2f\n", aco_info->beta);
  fprintf(problem->fileout, "  Rho....................................: %.2f\n", aco_info->rho);
  if (!(get_flag(problem->flags, FLAG_ANTS_RATIO)))
    fprintf(problem->fileout, "  Ants...................................: %i\n", aco_info->ants);
  else
    fprintf(problem->fileout, "  Ants...................................: %i (%i of %i - vertices)\n", aco_info->ants, aco_info->ratio, problem->nof_vertices);
  fprintf(problem->fileout, "  Pheromone scheme.......................: %s\n", schemes[aco_info->pheromone_scheme-1]);

  if (get_flag(problem->flags, FLAG_CHANGE_PHEROMONE_SCHEME))
    fprintf(problem->fileout, "  Change pheromone scheme after %i iterations.\n", aco_info->change_phero_scheme_iterations);

  if (get_flag(problem->flags, FLAG_MEMORY)) {
      if (!(get_flag(problem->flags, FLAG_MEMORY_RATIO))) {
         fprintf(problem->fileout, "  Memory Usage:\n\tMemory size......................: %i\n", aco_info->memory_size);
      }
      else {
         fprintf(problem->fileout, "  Memory Usage:\n\tMemory size......................: %i (%i of %i - ants)\n", aco_info->memory_size, aco_info->memory_ratio, aco_info->ants);
      }
    fprintf(problem->fileout, "\tDelta............................: %.2f\n", aco_info->delta);
  }

  if (aco_info->pheromone_scheme == PHEROMONE_SCHEME_3)
    fprintf(problem->fileout, "  Pheromone Scheme 3:\n\tGap..............................: %i\n", aco_info->gap);

  if (get_flag(problem->flags, FLAG_REUSE_COLOR)) {
    fprintf(problem->fileout, "  Try to reuse colors:\n\tX................................: %.2f\n", aco_info->x);
    fprintf(problem->fileout, "\tY................................: %.2f\n", aco_info->y);
  }

  if (get_flag(problem->flags, FLAG_CHANGE_ALPHA_BETA)) {
    fprintf(problem->fileout, "  Change alpha and beta:\n\tGamma............................: %.2f\n", aco_info->gamma);
    fprintf(problem->fileout, "\tOmega............................: %.2f\n", aco_info->omega);
    fprintf(problem->fileout, "\tChange alpha and beta after %i iterations.\n", aco_info->iterations_alpha_beta);
  }

}

void colorant_malloc(void) {

  aco_info              = malloc_(sizeof(aco_t));
  aco_info->alpha       = COLORANT_ALPHA;
  aco_info->beta        = COLORANT_BETA;
  aco_info->rho         = COLORANT_RHO;
  aco_info->ants        = COLORANT_ANTS;
  aco_info->ratio       = 0;
  aco_info->pheromone_scheme = PHEROMONE_SCHEME_1;
  aco_info->change_phero_scheme_iterations = COLORANT_CHANGE_PHERO_SCHEME_ITERATIONS;
  aco_info->memory_size = COLORANT_MEMORY_SIZE;
  aco_info->delta       = COLORANT_DELTA;
  aco_info->gap         = COLORANT_GAP;
  aco_info->gamma       = COLORANT_GAMMA;
  aco_info->omega       = COLORANT_OMEGA;
  aco_info->x           = COLORANT_X;
  aco_info->y           = COLORANT_Y;
  aco_info->iterations_alpha_beta = COLORANT_ITERATIONS_ALPHA_BETA;

}

void colorant_initialization(void) {

  aco_info->alpha_base = aco_info->alpha;
  aco_info->beta_base  = aco_info->beta;

  if (get_flag(problem->flags, FLAG_ANTS_RATIO)) {
        aco_info->ratio = aco_info->ants;
	aco_info->ants = (problem->nof_vertices * aco_info->ants) / 100;
  }

  if (get_flag(problem->flags, FLAG_MEMORY_RATIO)) {
        aco_info->memory_ratio = aco_info->memory_size;
        aco_info->memory_size = (aco_info->memory_size * aco_info->ants) / 100;
	aco_info->memory_size = aco_info->memory_size < 1 ? 1 : aco_info->memory_size;
  }

}

void colorant_show_solution(void) {
  if (get_flag(problem->flags, FLAG_CHANGE_ALPHA_BETA)) {
    fprintf(problem->fileout, "-------------------------------------------------\n");
    fprintf(problem->fileout, "Alpha..............................: %.2f\n", aco_info->alpha);
    fprintf(problem->fileout, "Beta...............................: %.2f\n", aco_info->beta);
    fprintf(problem->fileout, "Rho................................: %.2f\n", aco_info->rho);
  }
}


static int memory_length() {
  aco_memory_t *item = memory;
  int length = 0;
  for ( ; item != NULL; item = item->tail, length++);
  return length;
}

static void print_memory() {
  aco_memory_t *lmemory = memory;
  int v, count = 1;
  for (; lmemory; lmemory = lmemory->tail) {
    fprintf(problem->fileout, "Item: %i\n", count++);
    fprintf(problem->fileout, "No. of conflicting edges: %d\n", lmemory->head->nof_confl_edges);
    fprintf(problem->fileout, "No. of conflicting vertices: %d\n", lmemory->head->nof_confl_vertices);
    fprintf(problem->fileout, "Color:\n");
    for (v = 0; v < problem->nof_vertices; v++)
      fprintf(problem->fileout, "%i, ", lmemory->head->color_of[v]);
    fprintf(problem->fileout, "\n\n");
  }

  fprintf(problem->fileout, "Removed:\n");
  fprintf(problem->fileout, "No. of conflicting edges: %d\n", ant_memory_remove->nof_confl_edges);
  fprintf(problem->fileout, "No. of conflicting vertices: %d\n", ant_memory_remove->nof_confl_vertices);
  fprintf(problem->fileout, "Color:\n");
  for (v = 0; v < problem->nof_vertices; v++)
    fprintf(problem->fileout, "%i, ", ant_memory_remove->color_of[v]);
  fprintf(problem->fileout, "\n\n");

}

static void insert_into_memory(gcp_solution_t *sol) {
  aco_memory_t *item = malloc_(sizeof(aco_memory_t*));
  aco_memory_t *last, *previous;

  item->head = sol;
  item->tail = memory;
  memory = item;

  cpy_solution(sol, ant_memory_insert);

  if (aco_info->memory_size < memory_length()) {
    last = memory->tail;
    previous = memory;
    for ( ; last && last->tail != NULL; last = last->tail, previous = previous->tail);
    previous->tail = NULL;
    cpy_solution(last->head, ant_memory_remove);
    set_flag(problem->flags, FLAG_MEMORY_REMOVE);
    free(last);
  }

  //print_memory();

}

static void initialize_data(void) {

  int i, j, t;

  pheromone = malloc_(sizeof(double*) * problem->nof_vertices);

  phero_var = malloc_(sizeof(double**) * n_threads);

  for(t = 0; t < n_threads; t++){
      phero_var[t] = malloc_(sizeof(double*) * problem->nof_vertices);
      for (i = 0; i < problem->nof_vertices; i++) {
          phero_var[t][i] = malloc_(sizeof(double) * problem->nof_vertices);
          for (j = 0; j < problem->nof_vertices; j++) {
              phero_var[t][i][j] = 0;
          }
      }
  }

  for (i = 0; i < problem->nof_vertices; i++) {
    pheromone[i] = malloc_(sizeof(double) * problem->nof_vertices);
    for (j = 0; j < problem->nof_vertices; j++) {
      pheromone[i][j] = 0;
      if (!problem->adj_matrix[i][j]) {
	pheromone[i][j] = 1;
      }
    }
  }

  best_ant = malloc_(sizeof(gcp_solution_t));
  best_ant->color_of = malloc_(sizeof(int) * problem->nof_vertices);
  best_ant->nof_confl_vertices = INT_MAX;
  best_ant->nof_colors = problem->colors;
  best_ant->spent_time_ls = 0;

  best_colony = malloc_(sizeof(gcp_solution_t));
  best_colony->color_of = malloc_(sizeof(int) * problem->nof_vertices);
  best_colony->nof_confl_vertices = INT_MAX;
  best_colony->nof_colors = problem->colors;
  best_colony->spent_time_ls = 0;

  best_colonies = malloc_(sizeof(gcp_solution_t*) * n_threads);

  for(t = 0; t < n_threads; t++){
    best_colonies[t] = malloc_(sizeof(gcp_solution_t));
    best_colonies[t]->color_of = malloc_(sizeof(int) * problem->nof_vertices);
    best_colonies[t]->nof_confl_vertices = INT_MAX;
    best_colonies[t]->nof_colors = problem->colors;
    best_colonies[t]->spent_time_ls = 0;
    best_colonies[t]->id = t;
  }

  ant_k = malloc_(sizeof(gcp_solution_t));
  ant_k->color_of = malloc_(sizeof(int) * problem->nof_vertices);
  ant_k->nof_colors = problem->colors;
  ant_k->spent_time_ls = 0;

  if (get_flag(problem->flags, FLAG_MEMORY)) {
    ant_memory_remove = malloc_(sizeof(gcp_solution_t));
    ant_memory_remove->color_of = malloc_(sizeof(int) * problem->nof_vertices);
    ant_memory_remove->nof_colors = problem->colors;
    ant_memory_remove->spent_time_ls = 0;
    ant_memory_insert = malloc_(sizeof(gcp_solution_t));
    ant_memory_insert->color_of = malloc_(sizeof(int) * problem->nof_vertices);
    ant_memory_insert->nof_colors = problem->colors;
    ant_memory_insert->spent_time_ls = 0;

  }

  afk_initialize_data(aco_info->alpha, aco_info->beta);

}


/* Functions to help updating pheromone */
// as duas funcoes sao iguais
// head eh ant_memory_insert
// lembre que sempre inseri na cabeça (memory)
// LOGO eh a cabeça que sera usada no feromonio
//
// ant_memory_insert eh uma copia de head
static void update_pheromone_trails_memoryOLD(void) {
 int i, j;

 for (i = 0; i < problem->nof_vertices; i++) {
   for (j = 0; j < problem->nof_vertices; j++) {

     if (!problem->adj_matrix[i][j]) {
       if (memory->head->color_of[i] == memory->head->color_of[j])
	 pheromone[i][j] *=  1 + aco_info->delta;

       if ((get_flag(problem->flags, FLAG_MEMORY_REMOVE)) &&
	   (ant_memory_remove->color_of[i] == ant_memory_remove->color_of[j]))
	 pheromone[i][j] *=  1 - aco_info->delta;
     }
   }
 }
}

static void update_pheromone_trails_memory(void) {
 int i, j;

 for (i = 0; i < problem->nof_vertices; i++) {
   for (j = 0; j < problem->nof_vertices; j++) {

     if (!problem->adj_matrix[i][j]) {
       if (ant_memory_insert->color_of[i] == ant_memory_insert->color_of[j])
         pheromone[i][j] *=  1 + aco_info->delta;

       if ((get_flag(problem->flags, FLAG_MEMORY_REMOVE)) &&
           (ant_memory_remove->color_of[i] == ant_memory_remove->color_of[j]))
         pheromone[i][j] *=  1 - aco_info->delta;
     }
   }
 }
}

static void update_var_phero(gcp_solution_t *solution, int thread_id) {

  int i, j;

  for (i = 0; i < problem->nof_vertices; i++) {
    for (j = 0; j < problem->nof_vertices; j++) {
      if (!problem->adj_matrix[i][j] &&
	  (solution->color_of[i] == solution->color_of[j])) {
	phero_var[thread_id][i][j] += (solution->nof_confl_vertices == 0) ?
	  1 : 1.0/solution->nof_confl_vertices;
      }
    }
  }

}

static void update_pheromone_trails_scheme_1(void) {

  int i, j,t;

#if defined DEBUG
  fprintf(stderr, "update_pheromone_trails_scheme_1.\n");
#endif

  for(t = 0; t < n_threads; t++){
      for (i = 0; i < problem->nof_vertices; i++) {
          for (j = 0; j < problem->nof_vertices; j++) {
              if(t == 0) pheromone[i][j] *= aco_info->rho;
              pheromone[i][j] += phero_var[t][i][j];

              if(t == 0){
                  if (!problem->adj_matrix[i][j]) {
                      if (best_ant->color_of[i] == best_ant->color_of[j])
                      pheromone[i][j] += (best_ant->nof_confl_vertices == 0) ?
                      1 : 1.0/best_ant->nof_confl_vertices;
                      if (best_colony->color_of[i] == best_colony->color_of[j])
                      pheromone[i][j] += (best_colony->nof_confl_vertices == 0) ?
                      1 : 1.0/best_colony->nof_confl_vertices;
                  }
              }
              phero_var[t][i][j] = 0;
          }
      }
  }

}

static void update_pheromone_trails_scheme_2(void) {

  int i, j;

#if defined DEBUG
  fprintf(stderr, "update_pheromone_trails_scheme_2.\n");
#endif

  /*
  fprintf(problem->fileout, "ANTES\n");
  for (i = 0; i < problem->nof_vertices; i++) {
    for (j = 0; j < problem->nof_vertices; j++) {
      fprintf(problem->fileout, "%.2f, ", pheromone[i][j]);
    }
    fprintf(problem->fileout, "\n");
  }
  */
  for (i = 0; i < problem->nof_vertices; i++) {
    for (j = 0; j < problem->nof_vertices; j++) {
      pheromone[i][j] *= aco_info->rho;
      if (!problem->adj_matrix[i][j]) {
	if (best_ant->color_of[i] == best_ant->color_of[j]) {
	  pheromone[i][j] += (best_ant->nof_confl_vertices == 0) ?
	    1 : 1.0/best_ant->nof_confl_vertices;
	}
	if (best_colony->color_of[i] == best_colony->color_of[j]) {
	  pheromone[i][j] += (best_colony->nof_confl_vertices == 0) ?
	    1 : 1.0/best_colony->nof_confl_vertices;
	}
      }
    }
  }
  /*
  fprintf(problem->fileout, "DEPOIS\n");
  for (i = 0; i < problem->nof_vertices; i++) {
    for (j = 0; j < problem->nof_vertices; j++) {
      fprintf(problem->fileout, "%.2f, ", pheromone[i][j]);
    }
    fprintf(problem->fileout, "\n");
  }
  */
}

static void update_pheromone_trails_scheme_3() {

  int i, j;
  static int gap;

#if defined DEBUG
  fprintf(stderr, "update_pheromone_trails_scheme_3.\n");
#endif

  if (cycle % aco_info->gap == 0) gap = cycle / aco_info->gap;

  for (i = 0; i < problem->nof_vertices; i++) {
    for (j = i; j < problem->nof_vertices; j++) {
      pheromone[i][j] *= aco_info->rho;
      if (!problem->adj_matrix[i][j]) {
        if (gap) {
          if (best_ant->color_of[i] == best_ant->color_of[j])
            pheromone[i][j] += (best_ant->nof_confl_vertices == 0) ? 1 : 1.0/best_ant->nof_confl_vertices;
        }
        else {
          if (best_colony->color_of[i] == best_colony->color_of[j])
            pheromone[i][j] += (best_colony->nof_confl_vertices == 0) ? 1 : 1.0/best_colony->nof_confl_vertices;
        }
      }
      pheromone[j][i] = pheromone[i][j];
    }
  }
  gap--;

}

/* END Functions to help updating pheromone */


/* Vai ser transformado em função pthread */
static void *construct_solutions(void *r_solution){

  /* Inicialização dos dados locais da thread */
  gcp_solution_t *best_local_colony = (gcp_solution_t*) r_solution; //retorno da função
  gcp_solution_t *local_ant_k;
  gcp_solution_t *ant_memory;

  local_ant_k = malloc_(sizeof(gcp_solution_t));
  local_ant_k->color_of = malloc_(sizeof(int) * problem->nof_vertices);
  local_ant_k->nof_colors = problem->colors;
  local_ant_k->spent_time_ls = 0;

  int k;
  int ants = aco_info->ants;
  best_local_colony->f1 = INT_MAX * 1.0;

  for (k = 0; k < ants / n_threads; k++) {

    local_ant_k->total_cycles = cycle;

    ant_fixed_k(local_ant_k, pheromone, best_local_colony->id);

    /* Apply local search in all ants */
    if ( get_flag(problem->flags, FLAG_TABUCOL_ALL_ANTS) && (local_ant_k->nof_confl_vertices != 0) && (tabucol_info->cycles > 0) ) {
#if defined DEBUG
      //fprintf(stderr, "FLAG_TABUCOL_ALL_ANTS\n");
#endif
      tabucol(local_ant_k, tabucol_info->cycles, tabucol_info->tl_style, best_local_colony->id);
    }

    local_ant_k->spent_time = current_time_secs(TIME_FINAL, time_initial);

    if (local_ant_k->f1 < best_local_colony->f1) {
      cpy_solution(local_ant_k, best_local_colony);
      best_local_colony->cycles_to_best = cycle;
      best_local_colony->time_to_best = local_ant_k->spent_time;
    }

    if (aco_info->pheromone_scheme == PHEROMONE_SCHEME_1) update_var_phero(local_ant_k, best_local_colony->id);


    if ((best_local_colony->nof_confl_vertices == 0) || ((get_flag(problem->flags, FLAG_TIME)) && (problem->time <= current_time_secs(TIME_FINAL, time_initial)))) break;

  }

  if ( !(get_flag(problem->flags, FLAG_TABUCOL_ALL_ANTS)) && (best_local_colony->nof_confl_vertices != 0) && (tabucol_info->cycles > 0)) {
#if defined DEBUG
    fprintf(stderr, "FLAG_TABUCOL_BEST_ANT\n");
#endif
    tabucol(best_local_colony, tabucol_info->cycles, tabucol_info->tl_style, best_local_colony->id);
    best_local_colony->spent_time = current_time_secs(TIME_FINAL, time_initial);
  }

  if (get_flag(problem->flags, FLAG_MEMORY)) {
    ant_memory = malloc_(sizeof(gcp_solution_t));
    ant_memory->color_of = malloc_(sizeof(int) * problem->nof_vertices);
    cpy_solution(best_local_colony, ant_memory);
    insert_into_memory(ant_memory);
    }

  pthread_exit((void*) r_solution);
}

void printThreadsColonies(gcp_solution_t** colonies){
    int t;

    fprintf(problem->fileout, "\nThreads results:\n");
    for(t = 0; t < n_threads; t++){
        fprintf(problem->fileout, "T%d > Conflitos : %d\n", t, colonies[t]->nof_confl_vertices);
    }

    fprintf(problem->fileout, "\n");
}


gcp_solution_t* colorant(void) {
  void* pth_status;
  int tc; //Thread check
  int t = 0;
  int converg = 0;
  int change = 0;
  int cycle_phero = 0;

  initialize_data();
  best_ant->stop_criterion = 0;

  //threads
  pthread_t threads[n_threads];

  while (!terminate_conditions(best_ant, converg)) {

    cycle++;
    converg++;
    cycle_phero++;

    for(t = 0; t < n_threads; t++){
        tc = pthread_create(&threads[t], NULL, construct_solutions, best_colonies[t]);
        if (tc){
            printf("ERROR: return code from pthread_create() is %d\n", tc);
            exit(-1);
        }
    }

    for(t = 0; t < n_threads; t++){
        tc = pthread_join(threads[t], &pth_status);
        if (tc){
            printf("ERROR: return code from pthread_join() is %d\n", tc);
            exit(-1);
        }
    }

    sortSolutions(best_colonies);
    //if (get_flag(problem->flags, FLAG_VERBOSE)) printThreadsColonies(best_colonies);

    //Considera apenas a melhor solução de todas como best_colony
    cpy_solution(best_colonies[0], best_colony);

    if (best_colony->nof_confl_vertices < best_ant->nof_confl_vertices) {
      cpy_solution(best_colony, best_ant);
      best_ant->cycles_to_best = cycle;
      best_ant->time_to_best = best_colony->spent_time;
      converg = 0;
      change = 1;
    }

    switch (aco_info->pheromone_scheme){
    case PHEROMONE_SCHEME_1:
      update_pheromone_trails_scheme_1();
      break;
    case PHEROMONE_SCHEME_2:
      update_pheromone_trails_scheme_2();
      break;
    case PHEROMONE_SCHEME_3:
      update_pheromone_trails_scheme_3(cycle);
      break;
    }

    if (get_flag(problem->flags, FLAG_MEMORY)) {
      update_pheromone_trails_memory();
    }

    if (get_flag(problem->flags, FLAG_VERBOSE)) {
  fprintf(problem->fileout, "Cycle %d - Conflicts found: %d (edges), %d (vertices)\n", cycle, best_ant->nof_confl_edges, best_ant->nof_confl_vertices);
	fflush(stdout);
    }

    if (best_ant->nof_confl_vertices == 0) {
      best_ant->stop_criterion = STOP_BEST;
      break;
    }

    if ( get_flag(problem->flags, FLAG_CHANGE_ALPHA_BETA) && ((cycle % aco_info->iterations_alpha_beta)==0) ) {

      aco_info->gamma = change ? (1 - aco_info->omega) * aco_info->gamma:
                                (1 + aco_info->omega) * aco_info->gamma;

      aco_info->alpha = aco_info->alpha_base * aco_info->gamma;
      aco_info->beta  = aco_info->beta_base  * (1 - aco_info->gamma);

      if ( (aco_info->alpha <= 0) && ( aco_info->beta > 0) )
         aco_info->gamma = (1 + aco_info->omega) * aco_info->gamma;
      else
         if ( (aco_info->alpha > 0) && (aco_info->beta <= 0) )
            aco_info->gamma = (1 - aco_info->omega) * aco_info->gamma;

      aco_info->alpha = aco_info->alpha_base * aco_info->gamma;
      aco_info->beta  = aco_info->beta_base  * (1 - aco_info->gamma);

      change = 0;

      //printf("2 alfa: %.2f beta:%.2f gama:%.2f omega:%.2f\n", aco_info->alpha, aco_info->beta, aco_info->gamma, aco_info->omega);

    }

    if ( (get_flag(problem->flags, FLAG_DIFF_TABUCOL_SCHEME)) && ((cycle % tabucol_info->diff_scheme_iterations)==0) ) {
      tabucol_info->tl_style = tabucol_info->tl_style == TABUCOL_REACTIVE ? TABUCOL_DYNAMIC : TABUCOL_REACTIVE;
#if defined DEBUG
      printf("ACO LS: %i\n", tabucol_info->tl_style);
#endif
    }

    //printf("1 trocou: phero: %i cycle_phero: %i converg: %i \n", aco_info->pheromone_scheme, cycle_phero, converg);

    if ( (get_flag(problem->flags, FLAG_CHANGE_PHEROMONE_SCHEME)) && (converg >= aco_info->change_phero_scheme_iterations) && (cycle_phero >= aco_info->change_phero_scheme_iterations) ) {

      aco_info->pheromone_scheme++;
      if (aco_info->pheromone_scheme > 3)
	aco_info->pheromone_scheme = 1;

      cycle_phero = 0;
      //printf("trocou: phero: %i cycle: %i converg: %i \n", aco_info->pheromone_scheme, cycle, converg);

    }
  }

  best_ant->spent_time = current_time_secs(TIME_FINAL, time_initial);
  best_ant->total_cycles = cycle;

  return best_ant;

}
