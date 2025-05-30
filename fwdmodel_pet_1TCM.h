/**
 * fwdmodel_pet_1TCM.h
 * 
 * Implementation of one tissue compartment model for PET
 * Moss Zhao - Center for Advanced Functional Neuroimaging (CAFN), Stanford University

 */

/*  CCOPYRIGHT */
#pragma once

#include "fwdmodel_pet.h"

#include <fabber_core/fwdmodel.h>

#include <armawrap/newmat.h>

#include <string>
#include <vector>


class PET_1TCM_FwdModel : public PETFwdModel
{
public:
    static FwdModel *NewInstance();

    PET_1TCM_FwdModel()
    {
    }

    std::string GetDescription() const;
    void GetOptions(std::vector<OptionSpec> &opts) const;
    void Initialize(FabberRunData &rundata);
    void GetParameterDefaults(std::vector<Parameter> &params) const;
    void EvaluateModel(const NEWMAT::ColumnVector &params, NEWMAT::ColumnVector &result, const std::string &key = "") const;
    void GetOutputs(std::vector<std::string> &outputs) const;

protected:
    void Evaluate(const NEWMAT::ColumnVector &params, NEWMAT::ColumnVector &result) const;

private:
    double m_init_K1;
    double m_init_k2;
    
    /** Auto-register with forward model factory. */
    static FactoryRegistration<FwdModelFactory, PET_1TCM_FwdModel> registration;
};
