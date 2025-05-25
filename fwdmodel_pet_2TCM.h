/**
 * fwdmodel_pet_2TCM.h
 */

#pragma once

#include "fwdmodel_pet.h"

#include <fabber_core/fwdmodel.h>

#include <armawrap/newmat.h>

#include <string>
#include <vector>


class PET_2TCM_FwdModel : public PETFwdModel
{
public:
    static FwdModel *NewInstance();

    PET_2TCM_FwdModel()
    {
    }

    std::string GetDescription() const;
    void GetOptions(std::vector<OptionSpec> &opts) const;
    void Initialize(FabberRunData &rundata);
    void GetParameterDefaults(std::vector<Parameter> &params) const;
    void EvaluateModel(const NEWMAT::ColumnVector &params, NEWMAT::ColumnVector &result, const std::string &key = "") const;
    void ConvertParams(const NEWMAT::ColumnVector &params, NEWMAT::ColumnVector &result, const std::string &key = "") const;
    void GetOutputs(std::vector<std::string> &outputs) const;

protected:
     void Evaluate(const NEWMAT::ColumnVector &params, NEWMAT::ColumnVector &result) const;

private:
    // Initial values of model parameters - always inferred
    double m_init_alpha_1;
    double m_init_alpha_2;
    double m_init_beta_1;
    double m_init_beta_2;
    double m_ca;
    double m_lc;

    /** Auto-register with forward model factory. */
    static FactoryRegistration<FwdModelFactory, PET_2TCM_FwdModel> registration;
};
