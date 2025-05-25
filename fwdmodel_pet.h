/**
 * fwdmodel_pet.h
 *
 * Base class for PET models. Provides basic functions for dealing with AIF and
 * Moss Zhao - Center for Advanced Functional Neuroimaging (CAFN), Stanford University
 */

#pragma once

#include <fabber_core/fwdmodel.h>

#include <armawrap/newmat.h>

#include <string>
#include <vector>

using namespace NEWMAT;

/**
 * Base class for PET models as they share options
 *
 */
class PETFwdModel : public FwdModel
{
public:
    virtual ~PETFwdModel()
    {
    }

    virtual std::string ModelVersion() const;
    virtual void GetOptions(std::vector<OptionSpec> &opts) const;

    virtual void Initialize(FabberRunData &rundata);
    virtual void GetParameterDefaults(std::vector<Parameter> &params) const;
    virtual void InitVoxelPosterior(MVNDist &posterior) const;
    
    virtual Matrix convolve_matrix(const ColumnVector &kernel) const;
    virtual Matrix interp_matrix(const ColumnVector &x, const ColumnVector &x_p) const;
    
protected:

    ColumnVector m_kernel_time;
    ColumnVector m_aif_pet;
    Matrix m_c_mat;   
    double m_init_vB;
    double m_density;

};
