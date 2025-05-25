/**
 * fwdmodel_pet_1TCM.cc
 *
 * Implementation of one tissue compartment model for PET
 * Moss Zhao - Center for Advanced Functional Neuroimaging (CAFN), Stanford University
 */

/*  CCOPYRIGHT */

#include "fwdmodel_pet_1TCM.h"

#include <fabber_core/easylog.h>
#include <fabber_core/priors.h>

#include <miscmaths/miscprob.h>
#include <newimage/newimageall.h>
#include <armawrap/newmat.h>

#include <iostream>
#include <stdexcept>

using namespace std;
using namespace NEWMAT;

FactoryRegistration<FwdModelFactory, PET_1TCM_FwdModel> PET_1TCM_FwdModel::registration("pet_1TCM");

std::string PET_1TCM_FwdModel::GetDescription() const
{
    return "PET one tissue compartment model";
}

static OptionSpec OPTIONS[] = {
    {"init-K1", OPT_FLOAT, "Influx rate (mL/mL/s)", OPT_NONREQ,
     "0.05"},
    {"init-k2", OPT_FLOAT, "Efflux rate (1/s)", OPT_NONREQ,
     "0.05"},
    { "" },
};

void PET_1TCM_FwdModel::GetOptions(vector<OptionSpec> &opts) const
{
    PETFwdModel::GetOptions(opts);
    for (int i = 0; OPTIONS[i].name != ""; i++) {
        opts.push_back(OPTIONS[i]);
  }
}

void PET_1TCM_FwdModel::Initialize(FabberRunData &rundata)
{
    PETFwdModel::Initialize(rundata);
    
    // Initial values of model specific parameters
    m_init_K1 = rundata.GetDoubleDefault("init-K1", 0.05);
    m_init_k2 = rundata.GetDoubleDefault("init-k2", 0.05);
}

void PET_1TCM_FwdModel::GetParameterDefaults(std::vector<Parameter> &params) const
{
    PETFwdModel::GetParameterDefaults(params);
    
    int p = params.size();

    // specific parameters
    params.push_back(Parameter(p++, "K1", DistParams(m_init_K1, 100),
                               DistParams(m_init_K1, 100), PRIOR_NORMAL,
                               TRANSFORM_LOG()));
    params.push_back(Parameter(p++, "k2", DistParams(m_init_k2, 100),
                               DistParams(m_init_k2, 100), PRIOR_NORMAL,
                               TRANSFORM_LOG()));
}

void PET_1TCM_FwdModel::EvaluateModel(const ColumnVector &params, ColumnVector &result, const std::string &key) const
{
    if (key == ""){
        Evaluate(params, result);
    }
    else{
        result.ReSize(1);
        result(1) = params(2) * 6000.0 / m_density;
    } 
}


void PET_1TCM_FwdModel::Evaluate(const ColumnVector &params, ColumnVector &result) const
{
    // Parameters that are inferred - extract and give sensible names
    int p = 1;

    double vB = params(p++);
    double K1 = params(p++);
    double k2 = params(p++);


    ColumnVector exp_results = MISCMATHS::exp((-k2) * m_kernel_time);
    ColumnVector convolution_result = m_c_mat * exp_results;
   
    result = (1 - vB) * K1 * convolution_result + vB * m_aif_pet;

    for (int i = 1; i <= data.Nrows(); i++)
    {
        if (isnan(result(i)) || isinf(result(i)))
        {
            LOG << "Warning NaN or inf in result" << endl;
            LOG << "result: " << result.t() << endl;
            LOG << "params: " << params.t() << endl;

            result = 0.0;
            break;
        }
    }
}

void PET_1TCM_FwdModel::GetOutputs(std::vector<std::string> &outputs) const
{
    outputs.push_back("CBF");
}

FwdModel *PET_1TCM_FwdModel::NewInstance()
{
    return new PET_1TCM_FwdModel();
}
