/***********************************************************
 * Created: Seg 08 Ago 2011 17:01:23 BRT
 *
 * Author: Carla N. Lintzmayer, carla0negri@gmail.com
 *
 ***********************************************************
 *
 * Local search Tabucol [1987, Hertz]
 * Works with infeasible k-colorings.
 * Can be used with either a react or a dynamic tabu tenure.
 *
 * This program is an adaptation of the code available in
 * http://rose.epfl.ch/~bloechli/coloring/
 * explained below.
 *
 ***********************************************************/
/******************************************************************************/
//
//  ReactPartialCol, PartialCol, ReactTabucol and Tabucol graph coloring
//  heuristics. Reference code for the paper
//  "A Reactive Tabu Search Using Partial Solutions for the
//  Graph Coloring Problem" by Ivo Bloechliger and Nicolas Zuffery.
//
//  Copyright (C) 2003 - Ecole Poyltechnique Federale de Lausanne - EPFL,
//  Laboratory of Operations Research South Est-ROSE, 1015 Lausanne, Switzerland
//  Written by Ivo Bloechliger, Ivo.Bloechliger@epfl.ch
//  http://rose.epfl.ch/~bloechli/coloring/
//
/******************************************************************************/
//
// This program is distributed under the terms of the GNU General Public License
// as published by the Free Software Foundation. In paticular, this program is
// distributed WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. The EPFL shall in no
// case be liable for any damage of any kind in connection with the use of this
// program.  See the GNU General Public License for more details
// (http://www.gnu.org/copyleft/gpl.html#SEC1).
//
/******************************************************************************/

#include <stdio.h>
#include "color.h"
#include "tabucol.h"
#include "util.h"

void tabucol_printbanner(void) {/*{{{*/

  fprintf(problem->fileout, "TABUCOL\n");
  fprintf(problem->fileout, "-------------------------------------------------\n");
  if (tabucol_info->cycles > 0) {
    fprintf(problem->fileout, "Parameters:\n");
    fprintf(problem->fileout, "  Iterations.............................: %i\n", tabucol_info->cycles);
    if (get_flag(problem->flags, FLAG_TABUCOL_CONV)) fprintf(problem->fileout, "  Stop tabucol after %i cycles without improvement.\n", tabucol_info->convergence_cycles);

    if (get_flag(problem->flags, FLAG_CHANGE_TABUCOL_SCHEME)) {
      fprintf(problem->fileout, "  Change tabucol scheme after %d cycles.\n", tabucol_info->change_scheme_iterations);
      fprintf(problem->fileout, "\tInitial scheme: %s\n", tabucol_info->tl_style == TABUCOL_REACTIVE ? "Reactive" : "Dynamic");
    }

    if (get_flag(problem->flags, FLAG_DIFF_TABUCOL_SCHEME)) {
      fprintf(problem->fileout, "  Different tabucol scheme after %d cycles.\n", tabucol_info->change_scheme_iterations);
      fprintf(problem->fileout, "\tInitial scheme: %s\n", tabucol_info->tl_style == TABUCOL_REACTIVE ? "Reactive" : "Dynamic");
    }

    if ( (!(get_flag(problem->flags, FLAG_CHANGE_TABUCOL_SCHEME))) && (!(get_flag(problem->flags, FLAG_DIFF_TABUCOL_SCHEME))) ){
      if (tabucol_info->tl_style == TABUCOL_REACTIVE)
	fprintf(problem->fileout, "  Reactive scheme for tabu tenure.\n");
      else
	fprintf(problem->fileout, "  Dynamic scheme for tabu tenure.\n");
    }

    if (get_flag(problem->flags, FLAG_TABUCOL_ALL_ANTS))
      fprintf(problem->fileout, "  Apply tabu search on all ants.\n");
    else
      fprintf(problem->fileout, "  Apply tabu search only on the best ant.\n");

  }
  else
    fprintf(problem->fileout, "  Do not use tabu search.\n");

}/*}}}*/

void tabucol_malloc(void) {/*{{{*/

  tabucol_info = malloc_(sizeof(tabucol_t));
  tabucol_info->tl_style = TABUCOL_DYNAMIC;
  tabucol_info->cycles = TABUCOL_CYCLES;
  tabucol_info->convergence_cycles = TABUCOL_CONVERGENCE_CYCLES;
  tabucol_info->change_scheme_iterations = TABUCOL_CHANGE_SCHEME_ITERATIONS;
  tabucol_info->diff_scheme_iterations = TABUCOL_DIFF_SCHEME_ITERATIONS;
  tabucol_info->spent_time = 0;
}/*}}}*/

/* Number of conflicts for each color and node: (k+1)x(n+1) */
static int ***conflicts = NULL;
/* Tabu status for each color and node: (k+1)x(n) */
static int ***tabu_status = NULL;
/* Nodes that have conflicts. Position 0 keeps quantity: (k+1)x(n+1) */
static int **nodes_in_conflict = NULL;
/* Index of nodes in array nodes_in_conflicts */
static int **conf_position = NULL;

/* Aloca espaço para n_threads estruturas da busca tabu */
void tabucol_initialization(int n_threads) {
    int t;

    conflicts = malloc(sizeof(int**) * n_threads);
    tabu_status = malloc(sizeof(int**) * n_threads);
    nodes_in_conflict = malloc(sizeof(int*) * n_threads);
    conf_position = malloc(sizeof(int*) * n_threads);

    for(t = 0; t < n_threads; t++){
        nodes_in_conflict[t] = NULL;
    }
    //*conflicts = NULL;
    //*tabu_status = NULL;
    //*nodes_in_conflict = NULL;
    //*conf_position = NULL;
}

void tabucol_show_solution(void) {
  return;
}

static void initialize_arrays(gcp_solution_t *solution, int thread_id) {/*{{{*/

  int i, j, n;

  if (nodes_in_conflict[thread_id] == NULL) {
    nodes_in_conflict[thread_id] = malloc_(sizeof(int) * (problem->nof_vertices+1));
    conf_position[thread_id] = malloc_(sizeof(int) * problem->nof_vertices);

    conflicts[thread_id] = malloc_(sizeof(int*) * (problem->colors+1));
    tabu_status[thread_id] = malloc_(sizeof(int*) * (problem->colors+1));

    for (i = 0; i <= problem->colors; i++) {
      conflicts[thread_id][i] = malloc_(sizeof(int) * (problem->nof_vertices+1));
      tabu_status[thread_id][i] = malloc_(sizeof(int) * problem->nof_vertices);
    }
  }

  for (i = 0; i <= problem->colors; i++) {
    for (j = 0; j < problem->nof_vertices; j++) {
      conflicts[thread_id][i][j] = 0;
      tabu_status[thread_id][i][j] = 0;
    }
    conflicts[thread_id][i][problem->nof_vertices] = 0;
  }

  /* Initializing the conflicts array */
  if (get_flag(problem->flags, FLAG_ADJ_LIST)) {
    for (i = 0; i < problem->nof_vertices; i++) {
      for (j = 1; j <= problem->adj_list[i][0]; j++) {
	n = problem->adj_list[i][j];
	conflicts[thread_id][solution->color_of[n]][i]++;
      }
    }
  }
  else if (get_flag(problem->flags, FLAG_ADJ_MATRIX)) {
    for (i = 0; i < problem->nof_vertices; i++) {
      for (j = 0; j < problem->nof_vertices; j++) {
	if (problem->adj_matrix[i][j])
	  conflicts[thread_id][solution->color_of[j]][i]++;
      }
    }
  }

}/*}}}*/

static void neighbor_solution(int best_node, int best_color, /*{{{*/
			      gcp_solution_t *solution, int total_it, int t_tenure, int thread_id) {

  int j, last, n;
  int old_color;

  /* Change color of best_node */
  old_color = solution->color_of[best_node];
  solution->color_of[best_node] = best_color;
  if (get_flag(problem->flags, FLAG_S_PARTITION)) {
    for (j = 1; j <= solution->class_color[old_color][0]; j++) {
      if (solution->class_color[old_color][j] == best_node) {
	solution->class_color[old_color][j] =
	  solution->class_color[old_color][solution->class_color[old_color][0]];
	solution->class_color[old_color][0]--;
	break;
      }
    }
    solution->class_color[best_color][0]++;
    solution->class_color[best_color][solution->class_color[best_color][0]] = best_node;
  }

  /* If <old_color> was a conflict and <best_color> is not, remove <best_node>
   * from the list of conflicting nodes */
  if (conflicts[thread_id][old_color][best_node] && !(conflicts[thread_id][best_color][best_node])){
    last = nodes_in_conflict[thread_id][nodes_in_conflict[thread_id][0]];
    conf_position[thread_id][last] = conf_position[thread_id][best_node];
    nodes_in_conflict[thread_id][conf_position[thread_id][best_node]] =
      nodes_in_conflict[thread_id][nodes_in_conflict[thread_id][0]--];
  }
  else {
    /* If <old_color> was not a conflict and <best_color> is, put
     * <best_node> in the list */
    if (!(conflicts[thread_id][old_color][best_node]) &&
	conflicts[thread_id][best_color][best_node]) {
      nodes_in_conflict[thread_id][0]++;
      conf_position[thread_id][best_node] = nodes_in_conflict[thread_id][0];
      nodes_in_conflict[thread_id][conf_position[thread_id][best_node]] = best_node;
    }
  }

  /* Update the conflicts of the neighbors */
  if (get_flag(problem->flags, FLAG_ADJ_LIST)) {
    for (j = 1; j <= problem->adj_list[best_node][0]; j++) {

      n = problem->adj_list[best_node][j];

      /* Decrease the number of conflicts in the old color */
      conflicts[thread_id][old_color][n]--;

      if (conflicts[thread_id][old_color][n]==0 && solution->color_of[n]==old_color) {
	/* Remove <n> from the list of conflicting nodes if there are 0
	 * conflicts in its own color */
	last = nodes_in_conflict[thread_id][nodes_in_conflict[thread_id][0]];
	conf_position[thread_id][last] = conf_position[thread_id][n];
	nodes_in_conflict[thread_id][conf_position[thread_id][n]] =
	  nodes_in_conflict[thread_id][nodes_in_conflict[thread_id][0]--];
      }

      /* Increase the number of conflicts in the new color */
      conflicts[thread_id][best_color][n]++;

      if (conflicts[thread_id][best_color][n]==1 && solution->color_of[n]==best_color) {
	/* Add <n> in the list conflicting nodes if there is a new
	 * conflict in its own color */
	nodes_in_conflict[thread_id][0]++;
	conf_position[thread_id][n] = nodes_in_conflict[thread_id][0];
	nodes_in_conflict[thread_id][conf_position[thread_id][n]] = n;
      }
    }
  }
  else if (get_flag(problem->flags, FLAG_ADJ_MATRIX)) {
    for (n = 0; n < problem->nof_vertices; n++) {
      if (problem->adj_matrix[best_node][n]) {

	/* Decrease the number of conflicts in the old color */
	conflicts[thread_id][old_color][n]--;

	if (conflicts[thread_id][old_color][n]==0 && solution->color_of[n]==old_color) {
	  /* Remove <n> from the list of conflicting nodes if there are 0
	   * conflicts in its own color */
	  last = nodes_in_conflict[thread_id][nodes_in_conflict[thread_id][0]];
	  conf_position[thread_id][last] = conf_position[thread_id][n];
	  nodes_in_conflict[thread_id][conf_position[thread_id][n]] =
	    nodes_in_conflict[thread_id][nodes_in_conflict[thread_id][0]--];
	}

	/* Increase the number of conflicts in the new color */
	conflicts[thread_id][best_color][n]++;

	if (conflicts[thread_id][best_color][n]==1 && solution->color_of[n]==best_color) {
	  /* Add <n> in the list conflicting nodes if there is a new
	   * conflict in its own color */
	  nodes_in_conflict[thread_id][0]++;
	  conf_position[thread_id][n] = nodes_in_conflict[thread_id][0];
	  nodes_in_conflict[thread_id][conf_position[thread_id][n]] = n;
	}
      }
    }
  }

  /* Set the tabu status */
  tabu_status[thread_id][old_color][best_node] = total_it + t_tenure;

}/*}}}*/

static void free_(int thread_id) {/*{{{*/
  int i;
  for (i = 0; i <= problem->colors; i++) {
    free(conflicts[thread_id][i]);
    free(tabu_status[thread_id][i]);
  }

  free(nodes_in_conflict[thread_id]);
  free(conf_position[thread_id]);
  free(conflicts[thread_id]);
  free(tabu_status[thread_id]);

  conflicts[thread_id] = NULL;
  tabu_status[thread_id] = NULL;
  nodes_in_conflict[thread_id] = NULL;
  conf_position[thread_id] = NULL;

}/*}}}*/

void tabucol(gcp_solution_t *solution, int max_cycles, int type_of_tl, int thread_id) {/*{{{*/

#if defined DEBUG
  fprintf(stderr, "[INICIO] tabucol_info->tl_style: %i. %i\n", tabucol_info->tl_style,type_of_tl);
#endif

  /* declaring variables *//*{{{*/
  int i, c;
  int pairs[][3] = {
    {10000,10,5}, {10000,15,3}, {10000,5,10}, {5000,15,10},
    {5000,10,15}, {5000,5,20}, {1000,15,30}, {1000,10,50},
    {1000,5,100}, {500,5,100}, {500,10,150}, {500,15,200}
  };

  int num_pairs = sizeof(pairs)/sizeof(int)/3;

  int pair_cycles = 0;
  int frequency = pairs[0][0];
  int increment = pairs[0][1];
  int next_pair = pairs[0][2];

  int total_it = 0;
  int iteration = 0;
  int total_conflicts = 0;
  int best_solution_value;

  int tabu_tenure = problem->nof_vertices/10;

  int min_solution_value = problem->nof_vertices;
  int max_solution_value = 0;

  int nc, best_node, best_color, best_value, node;
  int max_min, p, new_value;

  double time_initial_ls = 0;

  //printf("\tTABUCOL: %i\n", type_of_tl);

  /*}}}*/

  initialize_arrays(solution, thread_id);

  /* Count the number of conflicts and set up the list nodes_in_conflict
   * with the associated list conf_position */
  nodes_in_conflict[thread_id][0] = 0;
  for (i = 0; i < problem->nof_vertices; i++) {
    if (conflicts[thread_id][solution->color_of[i]][i] > 0) {
      total_conflicts += conflicts[thread_id][solution->color_of[i]][i];
      nodes_in_conflict[thread_id][0]++;
      conf_position[thread_id][i] = nodes_in_conflict[thread_id][0];
      nodes_in_conflict[thread_id][conf_position[thread_id][i]] = i;
    }
  }
  total_conflicts /= 2;

  best_solution_value = total_conflicts;
  if (best_solution_value == 0) {
    return;
  }

  if (get_flag(problem->flags, FLAG_VERBOSE) && max_cycles < 0) {
    fprintf(problem->fileout, "Tabucol: arrays initialized; edges in conflict = %d\n, vertices in conflict = %d\n",
	    total_conflicts, nodes_in_conflict[thread_id][0]);
  }

  time_initial_ls = current_time_secs(TIME_INITIAL, 0);

  while (TRUE) {

    //printf("\t\tTABUCOL: %i - %i (%i)\n", type_of_tl, total_it, max_cycles);

    nc = nodes_in_conflict[thread_id][0];
    total_it++;
    iteration++;

    best_node = -1;
    best_color = -1;
    best_value = problem->nof_vertices * problem->nof_vertices;

    for (i = 1; i <= nodes_in_conflict[thread_id][0]; i++) {
      node = nodes_in_conflict[thread_id][i];

      /* Try for every node in conflict to move it to every color */
      for (c = 0; c < problem->colors; c++) {
	if (c != solution->color_of[node]) {
	  new_value = total_conflicts + conflicts[thread_id][c][node]
	    - conflicts[thread_id][solution->color_of[node]][node];
	  /* Take the new move with minimum value of new conflicts */
	  if (new_value <= best_value) {
	    if ((tabu_status[thread_id][c][node] < total_it) ||
		(new_value < best_solution_value)) {
	      best_node = node;
	      best_color = c;
	      best_value = new_value;
	    }
	  }
	}
      }
    }

    /* If no non-tabu moves have been found, take any random move */
    if (best_node == -1) {
#if defined LRAND
      //best_node = RANDOM(problem->nof_vertices);
      RANDOM(problem->buffer, best_node, int, problem->nof_vertices);
#elif defined NRAND
      //best_node = RANDOM(problem->seed, problem->nof_vertices);
      RANDOM(problem->seed, problem->buffer, best_node, int, problem->nof_vertices);
#endif
      while (conf_position[thread_id][best_node] == 0) {
#if defined LRAND
	//best_node = RANDOM(problem->nof_vertices);
	RANDOM(problem->buffer, best_node, int, problem->nof_vertices);
#elif defined NRAND
        //best_node = RANDOM(problem->seed, problem->nof_vertices);
        RANDOM(problem->seed, problem->buffer, best_node, int, problem->nof_vertices);
#endif
      }

#if defined LRAND
      //best_color = RANDOM(problem->colors);
      RANDOM(problem->buffer, best_color, int, problem->nof_vertices);
#elif defined NRAND
      //best_color = RANDOM(problem->seed, problem->colors);
      RANDOM(problem->seed, problem->buffer, best_color, int, problem->nof_vertices);
#endif

      /* Choose a color that is not the current color of <best_node> but
       * try this for just <|V|> number of tries */
      while (best_color == solution->color_of[best_node]) {
#if defined LRAND
	//best_color = RANDOM(problem->colors);
	RANDOM(problem->buffer, best_color, int, problem->nof_vertices);
#elif defined NRAND
        //best_color = RANDOM(problem->seed, problem->colors);
	RANDOM(problem->seed, problem->buffer, best_color, int, problem->nof_vertices);
#endif
      }
      best_value = total_conflicts + conflicts[thread_id][best_color][best_node]
	- conflicts[thread_id][solution->color_of[best_node]][best_node];
    }

    /* Execute the move */
    neighbor_solution(best_node, best_color, solution, total_it, tabu_tenure, thread_id);

    total_conflicts = best_value;

    if (get_flag(problem->flags, FLAG_TABUCOL_VERBOSE)) {
      fprintf(problem->fileout, "Tabucol: move executed = (%d,%d); edges in conflict = %d; vertices in conflict = %d\n",
	      best_node+1, best_color, total_conflicts, nodes_in_conflict[thread_id][0]);
      fprintf(problem->fileout, "\tCycle = %d; best solution so far = %d\n", total_it, nodes_in_conflict[thread_id][0] /*best_solution_value*/);

    }

    /* Calculating tabu_tenure: *//*{{{*/
    if (type_of_tl == TABUCOL_REACTIVE) {
      max_min = 0;
      /* Update the min and max objective function value */
      if (total_conflicts > max_solution_value)
	max_solution_value = total_conflicts;
      if (total_conflicts < min_solution_value)
	min_solution_value = total_conflicts;

      max_min = max_solution_value - min_solution_value;

      if (iteration % frequency == 0) {
	/* Adjust tabu_tenure every <frequency> iterations */
	if ((max_min < (4+total_conflicts/80)) || tabu_tenure == 0) {
	  tabu_tenure += increment;
	  if (pair_cycles == next_pair) {
#if defined LRAND
	    //p = (int) RANDOM(num_pairs);
            RANDOM(problem->buffer, p, int, num_pairs);
#elif defined NRAND
            //p = (int) RANDOM(problem->seed, num_pairs);
            RANDOM(problem->seed, problem->buffer, p, int, num_pairs);
#endif
	    frequency = pairs[p][0];
	    increment = pairs[p][1];
	    pair_cycles = 0;
	    next_pair = pairs[p][2];
	  }
	}
	else if (tabu_tenure) {
	  tabu_tenure--;
	}

	min_solution_value = problem->nof_vertices * problem->nof_vertices;
	max_solution_value = 0;

	if (pair_cycles == next_pair) {
#if defined LRAND
	  //p = RANDOM(num_pairs);
	  RANDOM(problem->buffer, p, int, num_pairs);
#elif defined NRAND
          //p = RANDOM(problem->seed, num_pairs);
          RANDOM(problem->seed, problem->buffer, p, int, num_pairs);
#endif
	  frequency = pairs[p][0];
	  increment = pairs[p][1];
	  pair_cycles = 0;
	  next_pair = pairs[p][2];
	}
	else {
	  pair_cycles++;
	}
      }
    }
    else {
      /* Using a dynamic tabu tenure */
#if defined LRAND
      //tabu_tenure = (int) (0.6 * nc) + RANDOM(10);
      RANDOM(problem->buffer, tabu_tenure, int, 10);
      tabu_tenure = (int) (0.6 * nc) + tabu_tenure;
#elif defined NRAND
      //tabu_tenure = (int) (0.6 * nc) + RANDOM(problem->seed, 10);
      RANDOM(problem->seed, problem->buffer, tabu_tenure, int, 10);
      tabu_tenure = (int) (0.6 * nc) + tabu_tenure;
#endif
    }/*}}}*/

#if defined DEBUG
    //fprintf(stderr, "total: %i best: %i.\n", total_conflicts,best_solution_value);
#endif

    /* Test if there is a new globally best solution */
    if (total_conflicts <= best_solution_value) {

      best_solution_value = total_conflicts;
      solution->time_to_best = current_time_secs(TIME_FINAL, time_initial);

      /* If all nodes are colored, stop iterating */
      if (best_solution_value == 0) {
	break;
      }

      /* Reinitialize some values */
      min_solution_value = problem->nof_vertices * problem->nof_vertices;
      max_solution_value = 0;
      iteration = 0;
    }

    if ( (get_flag(problem->flags, FLAG_CHANGE_TABUCOL_SCHEME)) && ((total_it % tabucol_info->change_scheme_iterations)==0) ) {
      type_of_tl = type_of_tl == TABUCOL_DYNAMIC ? TABUCOL_REACTIVE : TABUCOL_DYNAMIC;

#if defined DEBUG
      fprintf(stderr, "type_of_tl: %i. %i\n", type_of_tl,total_it);
#endif
    }

    if ((get_flag(problem->flags, FLAG_TIME)) && (problem->time <= current_time_secs(TIME_FINAL, time_initial))) {
       break;
    }
    else {
#if defined DEBUG
      //fprintf(stderr, "FLAG_TABUCOL_CONV: break. %i %i\n", iteration, tabucol_info->convergence_cycles);
#endif
       if ((get_flag(problem->flags,  FLAG_TABUCOL_CONV)) && (iteration >= tabucol_info->convergence_cycles)) {
#if defined DEBUG
	 fprintf(stderr, "FLAG_TABUCOL_CONV: break. %i\n", iteration);
#endif
	 break;
       }
       else {
          if (total_it >= max_cycles) {
             break;
          }
       }
    }

  }

  solution->nof_confl_edges = total_conflicts;
  solution->nof_confl_vertices = nodes_in_conflict[thread_id][0];
  solution->spent_time = current_time_secs(TIME_FINAL, time_initial);
  tabucol_info->spent_time += current_time_secs(TIME_FINAL, time_initial_ls) / n_threads;
  solution->spent_time_ls = tabucol_info->spent_time;
  solution->total_cycles = total_it;

  free_(thread_id);

}/*}}}*/
