#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include "crystal.h"
#include "lattice_generator_layer.h"
#include "randomgen.h"
#include "algebra.h"
//#include "spg_generation.h"
//#include "read_input.h"
#define LOWB 2.5  //lower bound for length of lattice vector
#define PI 3.141592653
#define epsilon_length 0.001
#define epsilon 0.1
#define MAX_LENGTH_RATIO 3.0	//legacy shape rule on film cell bases
#define SHAPE_TOL 1e-6		//keeps bases exactly on the ratio and angle limits
#define AREA_TOL 1e-4		//roundoff allowance on the film area window (A^2)


//min and max angles are in radians

static void fill_layer_lattice_vectors(float lv[3][3], float a[3], float b[3],
                                        float c, float c2x, float c2y)
{
	lv[0][0] = a[0];
	lv[1][1] = b[1];
	lv[2][2] = c;

	lv[0][1] = a[1];
	lv[0][2] = a[2];
	lv[1][2] = b[2];

	lv[1][0] = b[0];
	lv[2][0] = c2x;
	lv[2][1] = c2y;
}


void generate_oblique_triclinic(float lattice_vector[3][3],
	float target_volume, float max_angle, float min_angle,float tmp_lattice_vec_a [3],
        float tmp_lattice_vec_b [3],float tmp_lattice_vec_a_norm,float tmp_lattice_vec_b_norm,
	float gamma)


{	

	//printf("I am in generate_oblique_triclinic");
	int found =0;
	float cz = target_volume / (tmp_lattice_vec_a_norm * tmp_lattice_vec_b_norm * sin(gamma));
	float c;
	float coeff_cx;
	float coeff_cy;

	float random02 = uniform_dist_01();
	float beta = random02*(max_angle - min_angle) + min_angle;

	while (found ==0)
	{
		float random01 = uniform_dist_01();
		float alpha = random01*(max_angle - min_angle) + min_angle;
		float cz_square = pow(cz,2);
		coeff_cx=cos(beta);
		coeff_cy=(cos(alpha)-cos(beta)*cos(gamma))/sin(gamma);    //cosbeta2
		float coeff_c=1-pow(coeff_cx,2)-pow(coeff_cy,2);

		if (coeff_c > 0.2 && fabs(coeff_cy) < 1)   //set coeff_c > 0.2 to avoid lattice too oblique
		{
			c = sqrt(cz_square/coeff_c);
			found = 1;
		}

	}

	float cx = coeff_cx * c;
	float cy = coeff_cy * c;

	fill_layer_lattice_vectors(lattice_vector, tmp_lattice_vec_a, tmp_lattice_vec_b, cz, cx, cy);


	

}


void generate_oblique_monoclinic(float lattice_vector[3][3], 
	float target_volume, float max_angle, float min_angle,float tmp_lattice_vec_a [3],
        float tmp_lattice_vec_b [3],float tmp_lattice_vec_a_norm,float tmp_lattice_vec_b_norm,        float gamma)


{
	float c = target_volume / (tmp_lattice_vec_a_norm * tmp_lattice_vec_b_norm * sin(gamma));
	//printf("I am in generate_oblique_monoclinic and c is %f",c);
	
	fill_layer_lattice_vectors(lattice_vector, tmp_lattice_vec_a, tmp_lattice_vec_b, c, 0, 0);
	//printf("I am at end of generate_oblique_monoclinic");


}

void generate_rectangle_monoclinic(float lattice_vector[3][3], 
	float target_volume,float max_angle, float min_angle,float tmp_lattice_vec_a [3],
        float tmp_lattice_vec_b [3],float tmp_lattice_vec_a_norm,float tmp_lattice_vec_b_norm,	      float gamma)


{

	float random01 = uniform_dist_01();
	float alpha = random01*(max_angle - min_angle) + min_angle;
	float c = target_volume / (tmp_lattice_vec_a_norm * tmp_lattice_vec_b_norm * sin(alpha));
	float cy = c * cos(alpha);
	float cz = c * sin(alpha);

	fill_layer_lattice_vectors(lattice_vector, tmp_lattice_vec_a, tmp_lattice_vec_b, cz, 0, cy);
}


void generate_rectangle_orthorhombic(float lattice_vector[3][3], 
	float target_volume,float max_angle, float min_angle,float tmp_lattice_vec_a [3],
        float tmp_lattice_vec_b [3],float tmp_lattice_vec_a_norm,float tmp_lattice_vec_b_norm,	      float gamma)
{


	float c = target_volume / (tmp_lattice_vec_a_norm * tmp_lattice_vec_b_norm);

	fill_layer_lattice_vectors(lattice_vector, tmp_lattice_vec_a, tmp_lattice_vec_b, c, 0, 0);
}

void generate_hexagonal_hexagonal(float lattice_vector[3][3], 
	float target_volume,float max_angle, float min_angle,float tmp_lattice_vec_a [3],
        float tmp_lattice_vec_b [3],float tmp_lattice_vec_a_norm,float tmp_lattice_vec_b_norm,        float gamma)
{
	float c = target_volume / (sqrt(3)/2 * pow(tmp_lattice_vec_a_norm,2));
	//printf("c generated is %f\n",c);

	fill_layer_lattice_vectors(lattice_vector, tmp_lattice_vec_a, tmp_lattice_vec_b, c, 0, 0);
}

void generate_square_tetragonal(float lattice_vector[3][3], 
	float target_volume,float max_angle, float min_angle,float tmp_lattice_vec_a [3],
        float tmp_lattice_vec_b [3],float tmp_lattice_vec_a_norm,float tmp_lattice_vec_b_norm,        float gamma)
{
	//printf("I am in square tetragonal");
	//fflush(stdout);
	float c = target_volume / (pow(tmp_lattice_vec_a_norm,2));

	fill_layer_lattice_vectors(lattice_vector, tmp_lattice_vec_a, tmp_lattice_vec_b, c, 0, 0);
}


//metric tests of the layer lattice types, with the legacy tolerances.
//Returns the lattice type whose layer groups accept a basis with these
//lengths and angle (degrees), or -1 if none does.
static int metric_lattice_type(double norm_a, double norm_b, double gamma_deg)
{
	int equal = fabs(norm_a - norm_b) < epsilon_length;

	if (fabs(gamma_deg - 90) <= epsilon)
		return equal ? LAYER_SQUARE : LAYER_RECTANGULAR;
	if (fabs(gamma_deg - 120) <= epsilon)
		return equal ? LAYER_HEXAGONAL : -1;
	return equal ? -1 : LAYER_OBLIQUE;
}


int layer_lattice_type(int spg)
{
	if (spg < 1 || spg > 80)
		return -1;
	if (spg <= 7)
		return LAYER_OBLIQUE;
	if (spg <= 48)
		return LAYER_RECTANGULAR;
	if (spg <= 64)
		return LAYER_SQUARE;
	return LAYER_HEXAGONAL;
}


static double dot_dvec3(const double a[3], const double b[3])
{
	return a[0]*b[0] + a[1]*b[1] + a[2]*b[2];
}


static double cross_norm_dvec3(const double a[3], const double b[3])
{
	double c0 = a[1]*b[2] - a[2]*b[1];
	double c1 = a[2]*b[0] - a[0]*b[2];
	double c2 = a[0]*b[1] - a[1]*b[0];
	return sqrt(c0*c0 + c1*c1 + c2*c2);
}


//2x2 integer matrices are stored row-major {m00, m01, m10, m11} and act
//on row lattice vectors.
static void imat2_multiply(const int a[4], const int b[4], int c[4])
{
	int c0 = a[0]*b[0] + a[1]*b[2];
	int c1 = a[0]*b[1] + a[1]*b[3];
	int c2 = a[2]*b[0] + a[3]*b[2];
	int c3 = a[2]*b[1] + a[3]*b[3];
	c[0] = c0;
	c[1] = c1;
	c[2] = c2;
	c[3] = c3;
}


//inverse of a determinant +1 matrix
static void imat2_unimodular_inverse(const int a[4], int inv[4])
{
	inv[0] = a[3];
	inv[1] = -a[1];
	inv[2] = -a[2];
	inv[3] = a[0];
}


//rows of c = m (rows of v)
static void imat2_dvec_multiply(const int m[4], double v[2][3], double c[2][3])
{
	for (int k = 0; k < 3; k++)
	{
		c[0][k] = m[0]*v[0][k] + m[1]*v[1][k];
		c[1][k] = m[2]*v[0][k] + m[3]*v[1][k];
	}
}


//Row Hermite normal forms {{p, q}, {0, r}} with p r = det and 0 <= q < r,
//one per sublattice of index det. Writes them to out (4 ints each) unless
//out is NULL, and returns how many there are.
static int row_hnf(int det, int *out)
{
	int n = 0;
	for (int p = 1; p <= det; p++)
	{
		if (det % p)
			continue;
		int r = det / p;
		for (int q = 0; q < r; q++)
		{
			if (out)
			{
				out[4*n + 0] = p;
				out[4*n + 1] = q;
				out[4*n + 2] = 0;
				out[4*n + 3] = r;
			}
			n++;
		}
	}
	return n;
}


static int gcd_int(int a, int b)
{
	a = abs(a);
	b = abs(b);
	while (b)
	{
		int t = a % b;
		a = b;
		b = t;
	}
	return a;
}


//gcd of the 2x2 minors of [H_s | H_f]. It is 1 exactly when C = H_s S is
//the smallest common cell of S and F, and it does not change under basis
//changes of either cell.
static int gcd_of_minors(const int h_s[4], const int h_f[4])
{
	int row0[4] = {h_s[0], h_s[1], h_f[0], h_f[1]};
	int row1[4] = {h_s[2], h_s[3], h_f[2], h_f[3]};
	int g = 0;
	for (int i = 0; i < 4; i++)
		for (int j = i + 1; j < 4; j++)
			g = gcd_int(g, row0[i]*row1[j] - row0[j]*row1[i]);
	return g;
}


//Lagrange-Gauss reduction of the lattice spanned by the rows of v. Only
//determinant +1 steps are used, so the handedness is kept. Afterwards
//|v0| <= |v1| and |v0 . v1| <= |v0|^2 / 2, and t holds the integer matrix
//with v(reduced) = t v(input).
static void gauss_reduce(double v[2][3], int t[4])
{
	t[0] = 1;
	t[1] = 0;
	t[2] = 0;
	t[3] = 1;
	for (;;)
	{
		if (dot_dvec3(v[0], v[0]) > dot_dvec3(v[1], v[1]))
		{
			//(v0, v1) <- (v1, -v0)
			for (int k = 0; k < 3; k++)
			{
				double tmp = v[0][k];
				v[0][k] = v[1][k];
				v[1][k] = -tmp;
			}
			int t0 = t[0];
			int t1 = t[1];
			t[0] = t[2];
			t[1] = t[3];
			t[2] = -t0;
			t[3] = -t1;
		}
		int mu = (int) lround(dot_dvec3(v[0], v[1]) / dot_dvec3(v[0], v[0]));
		for (int k = 0; k < 3; k++)
			v[1][k] -= mu * v[0][k];
		t[2] -= mu * t[0];
		t[3] -= mu * t[1];
		if (dot_dvec3(v[1], v[1]) >= dot_dvec3(v[0], v[0]))
			break;
	}
}


typedef struct
{
	int i;			//vector = i a + j b in the reduced basis (a, b)
	int j;
	double vec[3];
	double norm;
} LATTICE_POINT;


//nonzero vectors of the lattice with reduced basis red, up to length radius.
//*points is allocated here; returns how many there are.
static int short_lattice_vectors(double red[2][3], double radius, LATTICE_POINT **points)
{
	double aa = dot_dvec3(red[0], red[0]);
	double ab = dot_dvec3(red[0], red[1]);
	double norm_a = sqrt(aa);
	double area = cross_norm_dvec3(red[0], red[1]);
	int capacity = 64;
	int n = 0;
	LATTICE_POINT *p = (LATTICE_POINT *) malloc(capacity * sizeof(LATTICE_POINT));

	//|j| times the height of b over a cannot exceed radius
	int jmax = (int) floor(radius * norm_a / area);
	for (int j = -jmax; j <= jmax; j++)
	{
		//the component along a, (i aa + j ab) / |a|, lies in [-radius, radius]
		int imin = (int) ceil((-radius * norm_a - j * ab) / aa);
		int imax = (int) floor((radius * norm_a - j * ab) / aa);
		for (int i = imin; i <= imax; i++)
		{
			if (i == 0 && j == 0)
				continue;
			double vec[3];
			for (int k = 0; k < 3; k++)
				vec[k] = i * red[0][k] + j * red[1][k];
			double norm = sqrt(dot_dvec3(vec, vec));
			if (norm > radius)
				continue;
			if (n == capacity)
			{
				capacity *= 2;
				p = (LATTICE_POINT *) realloc(p, capacity * sizeof(LATTICE_POINT));
			}
			p[n].i = i;
			p[n].j = j;
			for (int k = 0; k < 3; k++)
				p[n].vec[k] = vec[k];
			p[n].norm = norm;
			n++;
		}
	}
	*points = p;
	return n;
}


static void append_layer_lattice_combo(LAYER_LATTICE_LIST *list, int type,
	const LAYER_LATTICE_COMBO *combo)
{
	if (list->num_combos[type] == list->capacity[type])
	{
		list->capacity[type] = list->capacity[type] ? 2 * list->capacity[type] : 256;
		list->combos[type] = (LAYER_LATTICE_COMBO *) realloc(list->combos[type],
			list->capacity[type] * sizeof(LAYER_LATTICE_COMBO));
		if (!list->combos[type])
		{
			printf("***ERROR: build_layer_lattice_list: out of memory***\n");
			exit(1);
		}
	}
	list->combos[type][list->num_combos[type]] = *combo;
	list->num_combos[type]++;
}


//Adds to the list every basis F' of the film lattice F (rows of film) that
//the legacy shape rule accepts: right-handed, length ratio <=
//MAX_LENGTH_RATIO, angle in [min_angle, max_angle], and a basis some layer
//lattice type accepts. Every such basis has lengths <= 2 lambda_2 (lambda_1 <=
//lambda_2 the reduced lengths), since |a||b| = area / sin(angle) <= 2 area <=
//2 lambda_1 lambda_2. common = H_s S = H_f F. Returns how many were added.
static int add_film_bases(LAYER_LATTICE_LIST *list, const int h_s[4], const int h_f[4],
	double film[2][3], double common[2][3], double max_angle, double min_angle)
{
	double red[2][3];
	for (int r = 0; r < 2; r++)
		for (int k = 0; k < 3; k++)
			red[r][k] = film[r][k];
	int w[4];
	gauss_reduce(red, w);

	//the shorter vector of a basis is at most sqrt(2 lambda_1 lambda_2) long
	//and the longer one at least lambda_2, so no basis passes the ratio rule
	double lambda_1 = sqrt(dot_dvec3(red[0], red[0]));
	double lambda_2 = sqrt(dot_dvec3(red[1], red[1]));
	if (lambda_2 > 2 * MAX_LENGTH_RATIO * MAX_LENGTH_RATIO * lambda_1)
		return 0;

	LATTICE_POINT *points;
	int num_points = short_lattice_vectors(red, 2 * lambda_2 * (1 + 1e-9), &points);
	int num_added = 0;
	for (int k = 0; k < num_points; k++)
	{
		for (int l = 0; l < num_points; l++)
		{
			LATTICE_POINT *a = points + k;
			LATTICE_POINT *b = points + l;
			//a basis with the handedness of the reduced basis
			if (a->i * b->j - a->j * b->i != 1)
				continue;
			//limits are inclusive; hexagonal lattices have bases exactly at
			//30 and 150 degrees and length ratio 3
			double max_ratio = MAX_LENGTH_RATIO * (1 + SHAPE_TOL);
			if (a->norm > max_ratio * b->norm || b->norm > max_ratio * a->norm)
				continue;
			double gamma = atan2(cross_norm_dvec3(a->vec, b->vec), dot_dvec3(a->vec, b->vec));
			if (gamma < min_angle - SHAPE_TOL || gamma > max_angle + SHAPE_TOL)
				continue;
			int type = metric_lattice_type(a->norm, b->norm, gamma * 180 / PI);
			if (type < 0)
				continue;

			//F' = u F with u = (basis coefficients) w, so H_f' = H_f u^-1
			LAYER_LATTICE_COMBO combo;
			int coeff[4] = {a->i, a->j, b->i, b->j};
			int u[4];
			int u_inv[4];
			imat2_multiply(coeff, w, u);
			imat2_unimodular_inverse(u, u_inv);
			imat2_multiply(h_f, u_inv, combo.h_f);
			for (int i = 0; i < 4; i++)
				combo.h_s[i] = h_s[i];
			for (int i = 0; i < 3; i++)
			{
				combo.film[0][i] = a->vec[i];
				combo.film[1][i] = b->vec[i];
			}

			//H_f' F' must reproduce the common cell
			double check[2][3];
			imat2_dvec_multiply(combo.h_f, combo.film, check);
			for (int r = 0; r < 2; r++)
				for (int i = 0; i < 3; i++)
					if (fabs(check[r][i] - common[r][i]) > 1e-6 * (1 + fabs(common[r][i])))
					{
						printf("***ERROR: build_layer_lattice_list: H_f F != H_s S***\n");
						exit(1);
					}

			append_layer_lattice_combo(list, type, &combo);
			num_added++;
		}
	}
	free(points);
	return num_added;
}


int build_layer_lattice_list(LAYER_LATTICE_LIST *list, float lattice_vector_2d[2][3],
	int max_ds, int max_df, double area_min, double area_max,
	int max_interface_atoms, int num_substrate_atoms, int num_film_atoms,
	float max_angle, float min_angle)
{
	double substrate[2][3];
	for (int r = 0; r < 2; r++)
		for (int k = 0; k < 3; k++)
			substrate[r][k] = lattice_vector_2d[r][k];
	double substrate_area = cross_norm_dvec3(substrate[0], substrate[1]);

	//the defaults reproduce the legacy search F = M S with A_f <= area_max
	if (max_df <= 0)
		max_df = 1;
	if (max_ds <= 0)
	{
		max_ds = (int) floor(area_max / substrate_area);
		if (max_ds < 1)
			max_ds = 1;
	}

	for (int type = 0; type < NUM_LAYER_LATTICE_TYPES; type++)
	{
		list->combos[type] = NULL;
		list->num_combos[type] = 0;
		list->capacity[type] = 0;
	}
	list->max_ds = max_ds;
	list->max_df = max_df;
	list->substrate_area = substrate_area;
	list->num_film_lattices = 0;

	for (int d_s = 1; d_s <= max_ds; d_s++)
	{
		int num_hnf_s = row_hnf(d_s, NULL);
		int *hnf_s = (int *) malloc(4 * num_hnf_s * sizeof(int));
		row_hnf(d_s, hnf_s);

		for (int d_f = 1; d_f <= max_df; d_f++)
		{
			//the area window applies to the film cell, never to the common cell
			double film_area = d_s * substrate_area / d_f;
			if (film_area < area_min - AREA_TOL || film_area > area_max + AREA_TOL)
				continue;
			if ((long) d_s * num_substrate_atoms + (long) d_f * num_film_atoms > max_interface_atoms)
				continue;

			int num_hnf_f = row_hnf(d_f, NULL);
			int *hnf_f = (int *) malloc(4 * num_hnf_f * sizeof(int));
			row_hnf(d_f, hnf_f);

			for (int s = 0; s < num_hnf_s; s++)
			{
				for (int f = 0; f < num_hnf_f; f++)
				{
					//film side: transposed row HNF, so each film lattice
					//containing C at index d_f appears once
					int *h_s0 = hnf_s + 4*s;
					int h_f0[4] = {hnf_f[4*f + 0], hnf_f[4*f + 2], hnf_f[4*f + 1], hnf_f[4*f + 3]};
					if (gcd_of_minors(h_s0, h_f0) != 1)
						continue;

					//reduce C and carry the basis change to both matrices
					double common[2][3];
					int v[4];
					int h_s[4];
					int h_f[4];
					imat2_dvec_multiply(h_s0, substrate, common);
					gauss_reduce(common, v);
					imat2_multiply(v, h_s0, h_s);
					imat2_multiply(v, h_f0, h_f);

					//F = H_f^-1 C = adj(H_f) C / d_f
					int adj[4] = {h_f[3], -h_f[1], -h_f[2], h_f[0]};
					double film[2][3];
					imat2_dvec_multiply(adj, common, film);
					for (int r = 0; r < 2; r++)
						for (int k = 0; k < 3; k++)
							film[r][k] /= d_f;

					if (add_film_bases(list, h_s, h_f, film, common, max_angle, min_angle))
						list->num_film_lattices++;
				}
			}
			free(hnf_f);
		}
		free(hnf_s);
	}

	int total = 0;
	for (int type = 0; type < NUM_LAYER_LATTICE_TYPES; type++)
		total += list->num_combos[type];
	return total;
}


void free_layer_lattice_list(LAYER_LATTICE_LIST *list)
{
	for (int type = 0; type < NUM_LAYER_LATTICE_TYPES; type++)
	{
		free(list->combos[type]);
		list->combos[type] = NULL;
		list->num_combos[type] = 0;
		list->capacity[type] = 0;
	}
}


//film basis rotated so that a is along x and b lies in the xy plane
static void film_lower_triangle(const double film[2][3], float a[3], float b[3],
	float *norm_a, float *norm_b, float *gamma)
{
	double la = sqrt(dot_dvec3(film[0], film[0]));
	double lb = sqrt(dot_dvec3(film[1], film[1]));
	double g = atan2(cross_norm_dvec3(film[0], film[1]), dot_dvec3(film[0], film[1]));

	a[0] = la;
	a[1] = 0;
	a[2] = 0;
	b[0] = lb * cos(g);
	b[1] = lb * sin(g);
	b[2] = 0;
	*norm_a = la;
	*norm_b = lb;
	*gamma = g;
}


void generate_layer_lattice(const LAYER_LATTICE_LIST *list, crystal* random_crystal, int spg,
	float max_angle, float min_angle, float target_volume, int volume_multiplier)
{
	//no film cell of this lattice type on this rank: return a zero lattice,
	//which the caller treats as a failed attempt
	int type = layer_lattice_type(spg);
	if (type < 0 || list->num_combos[type] == 0)
	{
		for (int i = 0; i < 3; i++)
			for (int j = 0; j < 3; j++)
				random_crystal->lattice_vectors[i][j] = 0;
		return;
	}

	const LAYER_LATTICE_COMBO *combo =
		&list->combos[type][genrand_int32() % (unsigned long) list->num_combos[type]];
	for (int i = 0; i < 4; i++)
	{
		random_crystal->epitaxy_matrix[i] = combo->h_s[i];
		random_crystal->film_matrix[i] = combo->h_f[i];
	}

	float tmp_lattice_vec_a[3];
	float tmp_lattice_vec_b[3];
	float norm_tmp_lattice_vec_a;
	float norm_tmp_lattice_vec_b;
	float gamma;
	film_lower_triangle(combo->film, tmp_lattice_vec_a, tmp_lattice_vec_b,
		&norm_tmp_lattice_vec_a, &norm_tmp_lattice_vec_b, &gamma);
	float volume = target_volume * volume_multiplier;

	if (spg <= 2)  //oblique_triclinic
		generate_oblique_triclinic(random_crystal->lattice_vectors, volume, max_angle, min_angle,
			tmp_lattice_vec_a, tmp_lattice_vec_b, norm_tmp_lattice_vec_a, norm_tmp_lattice_vec_b, gamma);
	else if (spg <= 7)  //oblique_monoclinic
		generate_oblique_monoclinic(random_crystal->lattice_vectors, volume, max_angle, min_angle,
			tmp_lattice_vec_a, tmp_lattice_vec_b, norm_tmp_lattice_vec_a, norm_tmp_lattice_vec_b, gamma);
	else if (spg <= 18)  //rectangle_monoclinic
		generate_rectangle_monoclinic(random_crystal->lattice_vectors, volume, max_angle, min_angle,
			tmp_lattice_vec_a, tmp_lattice_vec_b, norm_tmp_lattice_vec_a, norm_tmp_lattice_vec_b, gamma);
	else if (spg <= 48)  //rectangle_orthorombic
		generate_rectangle_orthorhombic(random_crystal->lattice_vectors, volume, max_angle, min_angle,
			tmp_lattice_vec_a, tmp_lattice_vec_b, norm_tmp_lattice_vec_a, norm_tmp_lattice_vec_b, gamma);
	else if (spg <= 64)  //square_tetragonal
		generate_square_tetragonal(random_crystal->lattice_vectors, volume, max_angle, min_angle,
			tmp_lattice_vec_a, tmp_lattice_vec_b, norm_tmp_lattice_vec_a, norm_tmp_lattice_vec_b, gamma);
	else  //hexagonal_trigonal, hexagonal_hexagonal
		generate_hexagonal_hexagonal(random_crystal->lattice_vectors, volume, max_angle, min_angle,
			tmp_lattice_vec_a, tmp_lattice_vec_b, norm_tmp_lattice_vec_a, norm_tmp_lattice_vec_b, gamma);
}

//create a large volume lattice for testing compatiility
void generate_fake_layer_lattice(float lattice_vector[3][3], int spg)
{
	const float ax = 15;
	
	
	if(spg < 1 || spg > 80)
	printf("***ERROR: generate_lattice: spg out of bounds***");
	
	else if (spg <= 2)      //oblique_triclinic
	{	
		lattice_vector[0][0] = ax;
		lattice_vector[1][1] = 0.8*ax;
		lattice_vector[2][2] = 0.5*ax;

		lattice_vector[0][1] = 0;
		lattice_vector[0][2] = 0;
		lattice_vector[1][2] = 0;

		lattice_vector[1][0] = 0.8*ax*tan(10*PI/180);
		lattice_vector[2][0] = 0.5*ax*cos(85*PI/180);
		lattice_vector[2][1] = 0.5*ax*cos(70*PI/180);
	}
	
	else if ( spg <= 7  )    //oblique_monoclinic, alpha=beta=90
        {
                lattice_vector[0][0] = ax;
                lattice_vector[1][1] = 0.8*ax;
                lattice_vector[2][2] = 0.5*ax;

                lattice_vector[0][1] = ax * tan(70*PI/180);
                lattice_vector[0][2] = 0;
                lattice_vector[1][2] = 0;

                lattice_vector[1][0] = 0;
                lattice_vector[2][0] = 0;
                lattice_vector[2][1] = 0;
        }


	else if (  spg <= 18 )    //rectangle_monoclinic, gamma = beta = 90
	{
		lattice_vector[0][0] = ax;
		lattice_vector[1][1] = 0.8*ax;
		lattice_vector[2][2] = 0.5*ax;

		lattice_vector[0][1] = 0;
		lattice_vector[0][2] = 0;
		lattice_vector[1][2] = 0;

		lattice_vector[1][0] = 0;
		lattice_vector[2][0] = 0;
		lattice_vector[2][1] = 0.5*ax/tan(70*PI/180);
	}
	
	else if (spg <= 48)       //rectangle_orthorombic
	{
		lattice_vector[0][0] = ax;
		lattice_vector[1][1] = 0.8*ax;
		lattice_vector[2][2] = 0.5*ax;

		lattice_vector[0][1] = 0;
		lattice_vector[0][2] = 0;
		lattice_vector[1][2] = 0;

		lattice_vector[1][0] = 0;
		lattice_vector[2][0] = 0;
		lattice_vector[2][1] = 0;
	}
	
	else if (spg <= 64)     //square_tetragonal
	{
		lattice_vector[0][0] = ax;
		lattice_vector[1][1] = ax;
		lattice_vector[2][2] = 0.5*ax;

		lattice_vector[0][1] = 0;
		lattice_vector[0][2] = 0;
		lattice_vector[1][2] = 0;

		lattice_vector[1][0] = 0;
		lattice_vector[2][0] = 0;
		lattice_vector[2][1] = 0;
	}
	
	else if (spg <= 80)   //hexagonal_trigonnal, hexagonal_hexagonal
	{
		lattice_vector[0][0] = ax;
		lattice_vector[1][1] = ax*sin(120*PI/180);
		lattice_vector[2][2] = 0.5*ax;

		lattice_vector[0][1] = 0;
		lattice_vector[0][2] = 0;
		lattice_vector[1][2] = 0;

		lattice_vector[1][0] = ax*cos(120*PI/180);
		lattice_vector[2][0] = 0;
		lattice_vector[2][1] = 0;
	}
	
}

static inline float fmodulo (float n, float d)
{
	long q =n/d;
	float r = n - q*d ;
	return r;	
}
//need to change standardise lattice!!!!!
void standardise_layer_lattice( float lattice[3][3], int spg)
{
	//triclinic
	if (spg == 1 || spg == 2)
	{
		float bx = lattice[1][0];
		float ax = lattice[0][0];
		float bx_new = fmodulo (bx, ax);
		lattice[1][0] = bx_new;
		
		float by = lattice[1][1];
		float cy = lattice[2][1];
		float cy_new = fmodulo(cy, by);
		lattice[2][1] = cy_new;
		
		float cx = lattice[2][0];
		float cx_new = fmodulo(cx, ax);
		lattice[2][0] = cx_new;
	}
	
	/*
	else if (spg > 2 && spg < 16)
	{
		float bx = lattice[1][0];
		float ax = lattice[0][0];
		float bx_new = fmodulo (bx, ax);
		lattice[1][0] = bx_new;
		
		float cx = lattice[2][0];
		float cx_new = fmodulo(cx, ax);
		lattice[2][0] = cx_new;
		
	}
	*/
	
	else if (spg <= 0 || spg > 80)
		printf("***ERROR: lattice_generator: standardise_lattice invalid spg");
	
}

void get_lower_triangle(float lattice_vector[3][3],float lower_lattice_vector[3][3])
{
	
	float norm_lattice_vectors_a = sqrt(pow(lattice_vector[0][0],2)+pow(lattice_vector[0][1],2)+pow(lattice_vector[0][2],2));
	float norm_lattice_vectors_b = sqrt(pow(lattice_vector[1][0],2)+pow(lattice_vector[1][1],2)+pow(lattice_vector[1][2],2));
	//float norm_lattice_vectors_c = sqrt(pow(tmp_lattice_vectors[2][0],2)+pow(tmp_lattice_vectors[2][1],2)+
	//																		pow(tmp_lattice_vectors[2][2],2));
	float dot_product_ab = lattice_vector[0][0] * lattice_vector[1][0] + lattice_vector[0][1] * lattice_vector[1][1]+lattice_vector[0][2] * lattice_vector[1][2];
	//float dot_product_bc = tmp_lattice_vectors[1][0] * tmp_lattice_vectors[2][0] + tmp_lattice_vectors[1][1] * 
	//										tmp_lattice_vectors[2][1]+tmp_lattice_vectors[1][2] * tmp_lattice_vectors[2][2];
	//float dot_product_ac = tmp_lattice_vectors[0][0] * tmp_lattice_vectors[2][0] + tmp_lattice_vectors[0][1] * 
	//										tmp_lattice_vectors[2][1]+tmp_lattice_vectors[0][2] * tmp_lattice_vectors[2][2];

	//float alpha = acos(dot_product_bc/(norm_lattice_vectors_b * norm_lattice_vectors_c));
	//float beta = acos(dot_product_ac/(norm_lattice_vectors_a * norm_lattice_vectors_c));
	float gamma = acos(dot_product_ab/(norm_lattice_vectors_a * norm_lattice_vectors_b));

	lower_lattice_vector[0][0] = norm_lattice_vectors_a;
	lower_lattice_vector[0][1] = 0;
	lower_lattice_vector[0][2] = 0;
	lower_lattice_vector[1][0] = cos(gamma) * norm_lattice_vectors_b;
	lower_lattice_vector[1][1] = sin(gamma) * norm_lattice_vectors_b;
	lower_lattice_vector[1][2] = 0;
	//random_crystal->lattice_vectors[2][0] = cos(beta) * norm_lattice_vectors_c;
	//random_crystal->lattice_vectors[2][1] = (norm_lattice_vectors_b*norm_lattice_vectors_c*cos(alpha)-
	//	random_crystal->lattice_vectors[1][0] * random_crystal->lattice_vectors[2][0])/random_crystal->lattice_vectors[1][1];
    //random_crystal->lattice_vectors[2][2] = pow((pow(norm_lattice_vectors_c,2)-pow(random_crystal->lattice_vectors[2][0],2)- 
	//													pow(random_crystal->lattice_vectors[2][1],2)),0.5);


}











