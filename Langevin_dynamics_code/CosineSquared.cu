/*
GALAMOST - GPU-Accelerated Large-Scale Molecular Simulation Toolkit
COPYRIGHT
	GALAMOST Copyright (c) (2013) The group of Prof. Zhong-Yuan Lu
LICENSE
	This program is a free software: you can redistribute it and/or
	modify it under the terms of the GNU General Public License.
	This program is distributed in the hope that it will be useful,
	but WITHOUT ANY WARRANTY; without even the implied warranty of
	MERCHANT ABILITY or FITNESS FOR A PARTICULAR PURPOSE.
	See the General Public License v3 for more details.
	You should have received a copy of the GNU General Public License
	along with this program. If not, see <http://www.gnu.org/licenses/>.
DISCLAIMER
	The authors of GALAMOST do not guarantee that this program and its
	derivatives are free from error. In no event shall the copyright
	holder or contributors be liable for any indirect, incidental,
	special, exemplary, or consequential loss or damage that results
	from its use. We also have no responsibility for providing the
	service of functional extension of this program to general users.
USER OBLIGATION
	If any results obtained with GALAMOST are published in the scientific
	literature, the users have an obligation to distribute this program
	and acknowledge our efforts by citing the paper "Y.-L. Zhu, H. Liu,
	Z.-W. Li, H.-J. Qian, G. Milano, and Z.-Y. Lu, J. Comput. Chem. 2013,
	34, 2197-2211" in their article.
CORRESPONDENCE
	State Key Laboratory of Polymer Physics and Chemistry,
	Changchun Institute of Applied Chemistry, Chinese Academy of Sciences, China,
	Dr. You-Liang Zhu,
	Email: youliangzhu@ciac.ac.cn
*/
//	Maintainer: You-Liang Zhu

#include "CosineSquared.cuh"


Real4_tex_t pos_tex;
__global__ void gpu_compute_cosine_squared_kernel(Real4* d_force,
	ForceLog force_log,
	Real4* d_pos,
	BoxSize box,
	const unsigned int* d_n_neigh,
	const unsigned int* d_nlist,
	Index2D nli,
	Real4* d_params,
	int coeff_width,
	unsigned int Np,
	bool energy_shift)
{
	extern __shared__ Real4 s_params[];
	for (unsigned int cur_offset = 0; cur_offset < coeff_width * coeff_width; cur_offset += blockDim.x)
	{
		if (cur_offset + threadIdx.x < coeff_width * coeff_width)
			s_params[cur_offset + threadIdx.x] = d_params[cur_offset + threadIdx.x];
	}
	__syncthreads();

	unsigned int idx = blockIdx.x * blockDim.x + threadIdx.x;
	if (idx >= Np)
		return;

	unsigned int n_neigh = d_n_neigh[idx];
	Real4 pos = texFetchReal4(d_pos, pos_tex, idx);

	Real4 force = d_force[idx];
	Real virial = Real(0.0);
	Real6 virial_matrix = ToReal6(0.0, 0.0, 0.0, 0.0, 0.0, 0.0);
	if (force_log.virial)
		virial = force_log.d_virial[idx];
	if (force_log.virial_matrix)
		virial_matrix = force_log.d_virial_matrix[idx];

	unsigned int cur_neigh = 0;
	unsigned int next_neigh = d_nlist[nli(idx, 0)];

	for (int neigh_idx = 0; neigh_idx < n_neigh; neigh_idx++)
	{
		cur_neigh = next_neigh;
		next_neigh = d_nlist[nli(idx, neigh_idx + 1)];
		Real4 neigh_pos = texFetchReal4(d_pos, pos_tex, cur_neigh);

		Real dx = pos.x - neigh_pos.x;
		Real dy = pos.y - neigh_pos.y;
		Real dz = pos.z - neigh_pos.z;

		box.minDisImage(dx, dy, dz);
		Real rsq = dx * dx + dy * dy + dz * dz;
        Real r = sqrt_gala(rsq);
		int typ_pair = __real_as_int(neigh_pos.w) * coeff_width + __real_as_int(pos.w);

		Real epsilon = s_params[typ_pair].x;
		Real rcut = s_params[typ_pair].y;
		Real rcutsq = s_params[typ_pair].z;
		Real wc = s_params[typ_pair].w;
		Real wcsq = wc * wc;

		if (r < rcut)
		{

			Real pair_eng = -epsilon;
			force.w += pair_eng * Real(0.5);
		}

		if (r >= rcut && r <= (rcut+wc)  )
		{
			Real rinv = Real(1.0) / r;
			// Real r6inv = r2inv * r2inv * r2inv;
			// Real force_divr = r2inv * r6inv * (Real(12.0) * lj1 * r6inv - Real(6.0) * lj2);

			if (energy_shift)
            {
                
            }
			Real piw = Real(M_PI) / Real(2.0) / wc;
			Real deltarpiw = piw * (r - rcut);
            Real cosdeltarpiw = cos_gala(deltarpiw);
			Real force_divr = - 2 * epsilon * cosdeltarpiw * sin_gala(deltarpiw) * piw * rinv;

			//Real pair_eng = r6inv * (lj1 * r6inv - lj2);
			Real pair_eng = - epsilon * cosdeltarpiw * cosdeltarpiw;


			if (force_log.virial)
				virial += Real(1.0) / Real(6.0) * rsq * force_divr;
			if (force_log.virial_matrix)
			{
                Real force_div2r = Real(0.5) * force_divr;
				virial_matrix.x += dx * dx * force_div2r;   // press_tensor_xx
				virial_matrix.y += dx * dy * force_div2r;   // press_tensor_xy
				virial_matrix.z += dx * dz * force_div2r;   // press_tensor_xz
				virial_matrix.w += dy * dy * force_div2r;   // press_tensor_yy
				virial_matrix.m += dy * dz * force_div2r;   // press_tensor_yz
				virial_matrix.n += dz * dz * force_div2r;   // press_tensor_zz			
			}

			force.x += dx * force_divr;
			force.y += dy * force_divr;
			force.z += dz * force_divr;
			// energy is double counted: multiply by 0.5		
			force.w += pair_eng * Real(0.5);
		}
	}

	d_force[idx] = force;
	if (force_log.virial)
		force_log.d_virial[idx] = virial;
	if (force_log.virial_matrix)
		force_log.d_virial_matrix[idx] = virial_matrix;
}


cudaError_t gpu_compute_cosine_squared(Real4* d_force,
	ForceLog& force_log,
	Real4* d_pos,
	const BoxSize& box,
	const unsigned int* d_n_neigh,
	const unsigned int* d_nlist,
	const Index2D& nli,
	Real4* d_params,
	int coeff_width,
	int blocksize,
	unsigned int Np,
	unsigned int Ntot,
	bool energy_shift,
	unsigned int compute_capability)
{
	dim3 grid((int)ceil((Real)Np / (Real)blocksize), 1, 1);
	dim3 threads(blocksize, 1, 1);

	if (compute_capability < 350)
	{
		pos_tex.normalized = false;
		pos_tex.filterMode = cudaFilterModePoint;
		cudaError_t error = cudaBindTexture(0, pos_tex, d_pos, sizeof(Real4) * Ntot);
		if (error != cudaSuccess)
			return error;
	}

	gpu_compute_cosine_squared_kernel << < grid, threads, sizeof(Real4)* coeff_width* coeff_width >> > (d_force,
		force_log,
		d_pos,
		box,
		d_n_neigh,
		d_nlist,
		nli,
		d_params,
		coeff_width,
		Np,
		energy_shift);


	return cudaSuccess;
}


