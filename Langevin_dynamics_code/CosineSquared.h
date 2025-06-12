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

#include <boost/shared_ptr.hpp>

#include "Force.h"
#include "NeighborList.h"
#include "CosineSquared.cuh"

#ifndef __CONSINESQUARED_H__
#define __CONSINESQUARED_H__

class CosineSquared : public Force
{
public:
	CosineSquared(boost::shared_ptr<AllInfo> all_info, boost::shared_ptr<NeighborList> nlist, Real r_cut);

	virtual ~CosineSquared();

	virtual void setParams(const std::string& name1, const std::string& name2, Real epsilon,  Real wc);
	virtual void setParams(const std::string& name1, const std::string& name2, Real epsilon,  Real wc, Real r_cut);
	void setEnergy_shift();
	void setDispVirialCorr(bool dvc);
	void addInteractionType(unsigned int type);
protected:
	boost::shared_ptr<NeighborList> m_nlist;
	Real m_r_cut;
	bool m_energy_shift;
	boost::shared_ptr<Array<Real4> >  m_params;
	Real m_c6;
	unsigned int m_nc6;
	bool m_correct_press;
	bool m_first_compute;
	unsigned int m_np_interacted;
	std::vector<unsigned int> m_interacted_types;
	std::vector< bool > m_pair_param_set;
	bool m_pair_param_check;
	virtual void computeForce(unsigned int timestep);
	virtual void computeSlowForce(unsigned int timestep);
};

void export_CosineSquared();

#endif