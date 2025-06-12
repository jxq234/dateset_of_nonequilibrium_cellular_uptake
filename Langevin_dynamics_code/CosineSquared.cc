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

#include <boost/python.hpp>
using namespace boost::python;

#include "CosineSquared.h"
#include <stdexcept>


using namespace std;

CosineSquared::CosineSquared(boost::shared_ptr<AllInfo> all_info, boost::shared_ptr<NeighborList> nlist, Real r_cut)
	: Force(all_info), m_nlist(nlist), m_r_cut(r_cut), m_energy_shift(false)
{
	m_block_size = 320;
	Real nlistRcut = m_nlist->getRcut();

	if (r_cut < 0.0 || r_cut > nlistRcut)
	{
		cerr << endl << "***Error! The rcut is " << r_cut << " !" << endl << endl;
		throw runtime_error("Error building LjForce, rcut is negative or larger than the rcut of list");
	}
	m_c6 = 0.0;
	m_nc6 = 0;
	m_params = boost::shared_ptr<Array<Real4> >(new Array<Real4>(m_nkinds * m_nkinds, location::host));
	m_correct_press = false;
	m_first_compute = false;
	m_np_interacted = 0;
	m_pair_param_set.resize(m_nkinds * m_nkinds);
	m_pair_param_check = false;
	m_ObjectName = "CosineSquared";
	if (m_perf_conf->isRoot())
		cout << "INFO : " << m_ObjectName << " has been created" << endl;
}

CosineSquared::~CosineSquared()
{
}

void CosineSquared::computeSlowForce(unsigned int timestep)
{
	computeForce(timestep);
}

void CosineSquared::addInteractionType(unsigned int type)
{
	for (unsigned int i = 0; i < m_interacted_types.size(); i++)
	{
		if (m_interacted_types[i] == type)
			return;
	}
	m_interacted_types.push_back(type);
}

void CosineSquared::setParams(const std::string& name1, const std::string& name2, Real epsilon, Real wc)
{
	unsigned int typ1 = m_basic_info->switchNameToIndex(name1);
	unsigned int typ2 = m_basic_info->switchNameToIndex(name2);

	if (typ1 >= m_nkinds || typ2 >= m_nkinds)
	{
		cerr << endl << "***Error! Trying to set CosineSquared parameters for a non existed type! " << typ1 << "," << typ2 << endl << endl;
		throw runtime_error("CosineSquared::setParams argument error");
	}

	addInteractionType(typ1);
	addInteractionType(typ2);

	Real4* h_params = m_params->getArray(location::host, access::readwrite);

	h_params[typ1 * m_nkinds + typ2] = ToReal4(epsilon, m_r_cut, m_r_cut * m_r_cut, wc);
	h_params[typ2 * m_nkinds + typ1] = ToReal4(epsilon, m_r_cut, m_r_cut * m_r_cut, wc);
	m_pair_param_set[typ1 * m_nkinds + typ2] = true;
	m_pair_param_set[typ2 * m_nkinds + typ1] = true;
	m_pair_param_check = false;
}

void CosineSquared::setParams(const std::string& name1, const std::string& name2, Real epsilon, Real wc, Real r_cut)
{

	unsigned int typ1 = m_basic_info->switchNameToIndex(name1);
	unsigned int typ2 = m_basic_info->switchNameToIndex(name2);

	if (typ1 >= m_nkinds || typ2 >= m_nkinds)
	{
		cerr << endl << "***Error! Trying to set CosineSquared params for a non existant type! " << typ1 << "," << typ2 << endl << endl;
		throw runtime_error("CosineSquared::setParams argument error");
	}

	Real nlistRcut = m_nlist->getRcut();

	if (r_cut < 0.0 || r_cut > nlistRcut)
		throw runtime_error("Error CosineSquared setParams, negative rcut or larger than rcut of list");

    addInteractionType(typ1);
    addInteractionType(typ2);

	Real4* h_params = m_params->getArray(location::host, access::readwrite);

    h_params[typ1 * m_nkinds + typ2] = ToReal4(epsilon, r_cut, r_cut * r_cut, wc);
	h_params[typ2 * m_nkinds + typ1] = ToReal4(epsilon, r_cut, r_cut * r_cut, wc);
	m_pair_param_set[typ1 * m_nkinds + typ2] = true;
	m_pair_param_set[typ2 * m_nkinds + typ1] = true;
	m_pair_param_check = false;
}

void CosineSquared::setEnergy_shift()
{
	m_energy_shift = true;
}

void CosineSquared::setDispVirialCorr(bool dvc)
{
	m_correct_press = dvc;
}

void CosineSquared::computeForce(unsigned int timestep)
{
	if (!m_pair_param_check)
	{
		for (unsigned int i = 0; i < m_nkinds; i++)
			for (unsigned int j = i; j < m_nkinds; j++)
				if (!m_pair_param_set[i * m_nkinds + j])
					cerr << endl << "***Warning! CosineSquared, pair '" << m_basic_info->switchIndexToName(i) << "' and '" << m_basic_info->switchIndexToName(j) << "' has not been given parameters!" << endl << endl;
		m_pair_param_check = true;
	}
	m_nlist->compute(timestep);

	unsigned int Np = m_basic_info->getN();
	Real4* d_pos = m_basic_info->getPos()->getArray(location::device, access::read);
	const BoxSize& box = m_basic_info->getBox();
	Real4* d_force = m_basic_info->getForce()->getArray(location::device, access::readwrite);
	Real4* d_params = m_params->getArray(location::device, access::read);
	unsigned int compute_capability = m_perf_conf->getComputeCapability();

	ForceLog force_log;
	force_log.virial = m_all_info->getLogFlags()[log_flag::virial];
	force_log.potential = m_all_info->getLogFlags()[log_flag::potential];
	force_log.virial_matrix = m_all_info->getLogFlags()[log_flag::virial_matrix] || m_all_info->getLogFlags()[log_flag::press_tensor];
	force_log.d_virial = m_basic_info->getVirial()->getArray(location::device, access::readwrite);
	force_log.d_virial_matrix = m_basic_info->getVirialMatrix()->getArray(location::device, access::readwrite);

	gpu_compute_cosine_squared(d_force,
		force_log,
		d_pos,
		box,
		m_nlist->getGPUNNeigh(),
		m_nlist->getGPUNList(),
		m_nlist->getNListIndexer(),
		d_params,
		m_basic_info->getNKinds(),
		m_block_size,
		Np,
		m_basic_info->getN() + m_basic_info->getNGhosts(),
		m_energy_shift,
		compute_capability);

	CHECK_CUDA_ERROR();

}


void export_CosineSquared()
{
	scope in_lj = class_<CosineSquared, boost::shared_ptr<CosineSquared>, bases<Force>, boost::noncopyable >
		("CosineSquared", init< boost::shared_ptr<AllInfo>, boost::shared_ptr<NeighborList>, Real >())
		.def("setParams", static_cast<void (CosineSquared::*)(const std::string&, const std::string&, Real, Real)>(&CosineSquared::setParams))
		.def("setParams", static_cast<void (CosineSquared::*)(const std::string&, const std::string&, Real, Real, Real)>(&CosineSquared::setParams))
		.def("setEnergy_shift", &CosineSquared::setEnergy_shift)
		.def("setDispVirialCorr", &CosineSquared::setDispVirialCorr)
		;

}


