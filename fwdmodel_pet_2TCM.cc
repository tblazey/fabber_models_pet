/**
 * fwdmodel_pet_2TCM.cc
*/

#include "fwdmodel_pet_2TCM.h"
#include <fabber_core/easylog.h>
#include <fabber_core/priors.h>

#include <armawrap/newmat.h>
#include <miscmaths/miscprob.h>
#include <newimage/newimageall.h>

#include <iostream>
#include <stdexcept>

using namespace std;
using namespace NEWMAT;

FactoryRegistration<FwdModelFactory, PET_2TCM_FwdModel>
    PET_2TCM_FwdModel::registration("pet_2TCM");

std::string PET_2TCM_FwdModel::GetDescription() const {
  return "PET irreversible two compartment models";
}

static OptionSpec OPTIONS[] = {
    {"init-alpha-1", OPT_FLOAT, "Influx transfer rate (mL/s/mL)", OPT_NONREQ,
     "0.001"},
    {"init-alpha-2", OPT_FLOAT, "Influx transfer rate (mL/s/mL)", OPT_NONREQ,
     "0.001"},
    {"init-beta-1", OPT_FLOAT, "Influx transfer rate (mL/s/mL)", OPT_NONREQ,
     "0.0001"},
    {"init-beta-2", OPT_FLOAT, "Influx transfer rate (mL/s/mL)", OPT_NONREQ,
     "0.001"},
    {"ca", OPT_FLOAT, "Plasma blood glucose (mg/dL)", OPT_NONREQ, ""},
    {"lc", OPT_FLOAT, "Lumped constant", OPT_NONREQ, "0.81"},
    {""},
};

void PET_2TCM_FwdModel::GetOptions(vector<OptionSpec> &opts) const {
  PETFwdModel::GetOptions(opts);
  for (int i = 0; OPTIONS[i].name != ""; i++) {
    opts.push_back(OPTIONS[i]);
  }
}

void PET_2TCM_FwdModel::Initialize(FabberRunData &rundata) {
  PETFwdModel::Initialize(rundata);

  // Initial values of model specific parameters
  m_init_alpha_1 = rundata.GetDoubleDefault("init-alpha-1", 0.001);
  m_init_alpha_2 = rundata.GetDoubleDefault("init-alpha-2", 0.001);
  m_init_beta_1 = rundata.GetDoubleDefault("init-beta-1", 0.0001);
  m_init_beta_2 = rundata.GetDoubleDefault("init-beta-2", 0.001);
  m_ca = rundata.GetDoubleDefault("ca", 0);
  m_lc = rundata.GetDoubleDefault("lc", 0.81);
}

void PET_2TCM_FwdModel::GetParameterDefaults(
    std::vector<Parameter> &params) const {

  // generic model parameters
  PETFwdModel::GetParameterDefaults(params);
  int p = params.size();

  // specific parameters
  params.push_back(Parameter(p++, "alpha_1", DistParams(m_init_alpha_1, 100),
                             DistParams(m_init_alpha_1, 100), PRIOR_NORMAL,
                             TRANSFORM_LOG()));
  params.push_back(Parameter(p++, "alpha_2", DistParams(m_init_alpha_2, 100),
                             DistParams(m_init_alpha_2, 100), PRIOR_NORMAL,
                             TRANSFORM_LOG()));
  params.push_back(Parameter(p++, "beta_1", DistParams(m_init_beta_1, 100),
                             DistParams(m_init_beta_1, 100), PRIOR_NORMAL,
                             TRANSFORM_LOG()));
  params.push_back(Parameter(p++, "beta_2", DistParams(m_init_beta_2, 100),
                             DistParams(m_init_beta_2, 100), PRIOR_NORMAL,
                             TRANSFORM_LOG()));                            
                             
}

void PET_2TCM_FwdModel::EvaluateModel(const ColumnVector &params,
                                         ColumnVector &result,
                                         const std::string &key) const {
  if (key == "") {
    Evaluate(params, result);
  } else {
    ConvertParams(params, result, key);
  }
}

void PET_2TCM_FwdModel::ConvertParams(const ColumnVector &params,
                                      ColumnVector &result,
                                      const std::string &key) const {
    
    // extract parameters
    int p = 1;
    double vB = params(p++);
    double alpha_1 = params(p++);
    double alpha_2 = params(p++);
    double beta_1 = params(p++);
    double beta_2 = params(p++);
    
    // convert from optimized parameters to rate constants (see Hong and Fryer, Neuroimage,s 2010)
    double K1 = (alpha_1 + alpha_2) / (1.0 - vB) * 6000.0 / m_density;
    double k2 = (alpha_1 * beta_1 + alpha_2 * beta_2) / (alpha_1 + alpha_2);
    double k4 = beta_1 * beta_2 / k2;
    double k3 = beta_1 + beta_2 - k2 - k4;
    
    // convert from 1/seconds to 1/minutes
    k2 *= 60;
    k3 *= 60;
    k4 *= 60;
    
    // Macro parameters
    double Ki = K1 * k3 / (k2 + k3);
    double Vt = K1 / k2 * (1.0 + k3 / k4);
    double cmr;
    if (m_ca != 0){
        result.ReSize(7);
        cmr = Ki * m_ca / 18.0156 / m_lc;
    } else{
        result.ReSize(6);
    }
    
    // Save parameters
    result(1) = K1;
    result(2) = k2;
    result(3) = k3;
    result(4) = k4;
    result(5) = Ki;
    result(6) = Vt;
    if (m_ca != 0){
        result(7) = cmr;
    }
                                 
}

void PET_2TCM_FwdModel::Evaluate(const ColumnVector &params,
                                    ColumnVector &result) const {
  // Parameters that are inferred - extract and give sensible names
  int p = 1;

  double vB = params(p++);
  double alpha_1 = params(p++);
  double alpha_2 = params(p++);
  double beta_1 = params(p++);
  double beta_2 = params(p++);

  ColumnVector c_1 = alpha_1 * m_c_mat * MISCMATHS::exp(-beta_1 * m_kernel_time);
  ColumnVector c_2 = alpha_2 * m_c_mat * MISCMATHS::exp(-beta_2 * m_kernel_time);

  result = c_1 + c_2 + vB * m_aif_pet;

  for (int i = 1; i <= data.Nrows(); i++) {
    if (isnan(result(i)) || isinf(result(i))) {
      LOG << "Warning NaN or inf in result" << endl;
      LOG << "result: " << result.t() << endl;
      LOG << "params: " << params.t() << endl;

      result = 0.0;
      break;
    }
  }
}

void PET_2TCM_FwdModel::GetOutputs(std::vector<std::string> &outputs) const {
  outputs.push_back("rates");
}

FwdModel *PET_2TCM_FwdModel::NewInstance() {
  return new PET_2TCM_FwdModel();
}
