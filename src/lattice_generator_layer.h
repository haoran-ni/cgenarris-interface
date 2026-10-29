#ifndef LATTICE_GENERATOR_LAYER_H
#define LATTICE_GENERATOR_LAYER_H

#include "crystal.h"

//layer lattice types; layer groups 1-7 are oblique, 8-48 rectangular,
//49-64 square and 65-80 hexagonal
enum
{
	LAYER_OBLIQUE,
	LAYER_RECTANGULAR,
	LAYER_SQUARE,
	LAYER_HEXAGONAL,
	NUM_LAYER_LATTICE_TYPES
};

//A film cell F' and its common cell with the substrate cell S,
//C = H_s S = H_f F'. Matrices are row-major and act on row lattice vectors.
typedef struct
{
	int h_s[4];
	int h_f[4];
	double film[2][3];	//F' in the Cartesian frame of S
} LAYER_LATTICE_COMBO;

//film cells to sample from, one list per layer lattice type
typedef struct
{
	LAYER_LATTICE_COMBO *combos[NUM_LAYER_LATTICE_TYPES];
	int num_combos[NUM_LAYER_LATTICE_TYPES];
	int capacity[NUM_LAYER_LATTICE_TYPES];
	int max_ds;		//bounds used, after defaults
	int max_df;
	double substrate_area;
	int num_film_lattices;	//distinct film lattices F with at least one cell
} LAYER_LATTICE_LIST;


//void gen_triclinic_lattice(float lattice_vector[3][3], 
//	float target_volume, float max_angle, float min_angle);
//void gen_monoclinic_lattice(float lattice_vector[3][3], 
//	float target_volume, float max_angle, float min_angle);
/*
void gen_orthorhombic_lattice(float lattice_vector[3][3], 
	float target_volume);
void gen_tetragonal_lattice(float lattice_vector[3][3],
	float target_volume);
void gen_hexagonal_lattice(float lattice_vector[3][3], 
	float target_volume);
void gen_cubic_lattice(float lattice_vector[3][3],
	float target_volume);

*/
void generate_oblique_triclinic(float lattice_vector[3][3], 
	float target_volume, float max_angle, float min_angle,float tmp_lattice_vec_a [3],
        float tmp_lattice_vec_b [3],float tmp_lattice_vec_a_norm,float tmp_lattice_vec_b_norm,float gamma);

void generate_oblique_monoclinic(float lattice_vector[3][3], 
	float target_volume, float max_angle, float min_angle,float tmp_lattice_vec_a [3],
        float tmp_lattice_vec_b [3],float tmp_lattice_vec_a_norm,float tmp_lattice_vec_b_norm,float gamma);

void generate_rectangle_monoclinic(float lattice_vector[3][3], 
	float target_volume,float max_angle, float min_angle,float tmp_lattice_vec_a [3],
        float tmp_lattice_vec_b [3],float tmp_lattice_vec_a_norm,float tmp_lattice_vec_b_norm,float gamma);

void generate_rectangle_orthorhombic(float lattice_vector[3][3], 
	float target_volume,float max_angle, float min_angle,float tmp_lattice_vec_a [3],
                                                        float tmp_lattice_vec_b [3],float tmp_lattice_vec_a_norm,float tmp_lattice_vec_b_norm,float gamma);

void generate_hexagonal_hexagonal(float lattice_vector[3][3], 
	float target_volume,float max_angle, float min_angle,float tmp_lattice_vec_a [3],
                                                       float tmp_lattice_vec_b [3],float tmp_lattice_vec_a_norm,float tmp_lattice_vec_b_norm,float gamma);

void generate_square_tetragonal(float lattice_vector[3][3], 
	float target_volume,float max_angle, float min_angle,float tmp_lattice_vec_a [3],
                                                        float tmp_lattice_vec_b [3],float tmp_lattice_vec_a_norm,float tmp_lattice_vec_b_norm,float gamma);


void generate_layer_lattice(const LAYER_LATTICE_LIST *list, crystal* random_crystal, int spg,
	float max_angle, float min_angle, float target_volume, int volume_multiplier);

void generate_fake_layer_lattice(float lattice_vector[3][3], int spg);

void standardise_layer_lattice(float lattice[3][3], int spg);

int build_layer_lattice_list(LAYER_LATTICE_LIST *list, float lattice_vector_2d[2][3],
	int max_ds, int max_df, double area_min, double area_max,
	int max_interface_atoms, int num_substrate_atoms, int num_film_atoms,
	float max_angle, float min_angle);

void free_layer_lattice_list(LAYER_LATTICE_LIST *list);

int layer_lattice_type(int spg);


void get_lower_triangle(float lattice_vector[3][3],float lower_lattice_vector[3][3]);
#endif

