/**
 * fwdmodel_pet.h
 *
 * Base class for PET models. Provides basic functions for dealing with AIF and
 * Moss Zhao - Center for Advanced Functional Neuroimaging (CAFN), Stanford University
 */

#pragma once

#include "fabber_core/fwdmodel.h"

#include <newmat.h>

#include <string>
#include <vector>

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
    
protected:
    // Mandatory PET configuration
    // Inference flags
    bool m_infer_delay, m_infer_t10, m_infer_sig0;
    
    // AIF as concentration curve
    std::string m_aif_type;
    NEWMAT::ColumnVector m_aif, m_time;

    // Other options
    mutable int m_sig0_idx;
    bool m_auto_init_delay;

    double m_dt;

};
