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

#include "NeighborList.cuh"
#include "Info.cuh"


#ifndef __COSINESQUARED_CUH__
#define __COSINESQUARED_CUH__



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
	unsigned int compute_capability);

#endif


