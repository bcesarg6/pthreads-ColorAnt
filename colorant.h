
/*************************************************************************
 * Created: Qua 27 Jan 2011 15:17:39 BRST
 *
 * Author: Carla N. Lintzmayer, carla0negri@gmail.com
 *
 * VERSAO 1:
 * * Todas as formigas reforçam o feromônio
 * * A melhor formiga da colônia reforça o feromônio
 * * A melhor formiga global reforça o feromônio
 * * A busca local é a React-Tabucol
 * * A busca local é aplicada apenas na melhor formiga da colônia
 *
 * VERSAO 2:
 * * A melhor formiga da colônia reforça o feromônio
 * * A melhor formiga global reforça o feromônio
 * * A busca local é a React-Tabucol
 * * A busca local é aplicada apenas na melhor formiga da colônia
 *
 * VERSAO 3:
 * * A melhor formiga da colônia reforça o feromônio
 * * A melhor formiga global reforça o feromônio
 * * Em cada ciclo, ou a melhor formiga da colônia ou a melhor formiga global
 * * reforçam o feromônio; nunca as duas juntas
 * * A busca local é a React-Tabucol
 * * A busca local é aplicada em todas as formigas da colônia
 *
 * VERSAO 4:
 * * A melhor formiga da colônia reforça o feromônio
 * * A melhor formiga global reforça o feromônio
 * * Em cada ciclo, ou a melhor formiga da colônia ou a melhor formiga global
 * * reforçam o feromônio; nunca as duas juntas
 * * A busca local é a React-Tabucol
 * * A busca local é aplicada em todas as formigas da colônia
 * * O algoritmo tenta reutilizar cores
 *
 * VERSAO 5:
 * * A melhor formiga da colônia reforça o feromônio
 * * A melhor formiga global reforça o feromônio
 * * Em cada ciclo, ou a melhor formiga da colônia ou a melhor formiga global
 * * reforçam o feromônio; nunca as duas juntas
 * * A busca local é a React-Tabucol
 * * A busca local é aplicada em todas as formigas da colônia
 * * O algoritmo ajusta alfa e beta
 *
 * VERSAO 6: (VERSAO 1 + VERSAO 4)
 * * Todas as formigas reforçam o feromônio
 * * A melhor formiga da colônia reforça o feromônio
 * * A melhor formiga global reforça o feromônio
 * * A busca local é a React-Tabucol
 * * A busca local é aplicada apenas na melhor formiga da colônia
 * * O algoritmo tenta reutilizar cores
 * *
 *
 * VERSAO 7: (VERSAO 3, EXCETO QUE A BUSCA LOCAO EH APENAS APLICADA A MELHOR FORMIGA)
 * * A melhor formiga da colônia reforça o feromônio
 * * A melhor formiga global reforça o feromônio
 * * Em cada ciclo, ou a melhor formiga da colônia ou a melhor formiga global
 * * reforçam o feromônio; nunca as duas juntas
 * * A busca local é a React-Tabucol
 * * A busca local é aplicada apenas na melhor formiga da colônia
 *
*************************************************************************/

#ifndef __COLORANT_H
#define __COLORANT_H

#include "color.h"

#define COLORANT_ALPHA            2
#define COLORANT_BETA             8
#define COLORANT_RHO           0.60
#define COLORANT_ANTS          200
#define COLORANT_MEMORY_SIZE           25
#define COLORANT_DELTA          0.5

#define COLORANT_GAP             10
#define COLORANT_X              1.0
#define COLORANT_Y              2.0

#define COLORANT_GAMMA          0.4
#define COLORANT_OMEGA          0.2
#define COLORANT_ITERATIONS_ALPHA_BETA      25
#define COLORANT_PHEROMONE_SCHEME 1
#define COLORANT_CHANGE_PHERO_SCHEME_ITERATIONS 10


/* pheromone scheme 1 */
#define PHEROMONE_SCHEME_1 1

/* pheromone scheme 2 */
#define PHEROMONE_SCHEME_2 2

/* pheromone scheme 3 */
#define PHEROMONE_SCHEME_3 3

struct aco_t {

  float alpha_base;
  float beta_base;
  float rho_base;
  int iterations_alpha_beta;
  int gap;
  float x;
  float y;
  float gamma;
  float omega;
  int memory_size;
  int memory_ratio;
  float delta;
  int ants;
  int ratio;
  float alpha;
  float beta;
  float rho;
  int pheromone_scheme;
  int change_phero_scheme_iterations;
};

typedef struct aco_memory_t aco_memory_t;

struct aco_memory_t {
  gcp_solution_t *head;
  aco_memory_t *tail;
};


typedef struct aco_t aco_t;

aco_t *aco_info;

void colorant_printbanner(void);
void colorant_malloc(void);
void colorant_initialization(void);
void colorant_show_solution(void);
gcp_solution_t* colorant(void);

#endif /* __COLORANT1_H */

