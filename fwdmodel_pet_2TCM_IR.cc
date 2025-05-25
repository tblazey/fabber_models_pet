/**
 * fwdmodel_pet_2TCM_IR.cc
 */

#include "fwdmodel_pet_2TCM_IR.h"
#include <fabber_core/easylog.h>
#include <fabber_core/priors.h>

#include <armawrap/newmat.h>
#include <miscmaths/miscprob.h>
#include <newimage/newimageall.h>

#include <iostream>
#include <stdexcept>

using namespace std;
using namespace NEWMAT;

FactoryRegistration<FwdModelFactory, PET_2TCM_IR_FwdModel>
    PET_2TCM_IR_FwdModel::registration("pet_2TCM_IR");

std::string PET_2TCM_IR_FwdModel::GetDescription() const {
  return "PET irreversible two compartment models";
}

static OptionSpec OPTIONS[] = {
    {"init-K1", OPT_FLOAT, "Influx transfer rate (mL/s/mL)", OPT_NONREQ,
     "0.05"},
    {"init-Ki", OPT_FLOAT, "Net transfer rate (mL/s/mL)", OPT_NONREQ,
     "0.05"},
    {"init-k-sum", OPT_FLOAT, "Sum of k2 and k3 (1/s)", OPT_NONREQ,
     "0.05"},
    {"ca", OPT_FLOAT, "Plasma blood glucose (mg/dL)", OPT_NONREQ, ""},
    {"lc", OPT_FLOAT, "Lumped constant", OPT_NONREQ, "0.81"},
    {""},
};

void PET_2TCM_IR_FwdModel::GetOptions(vector<OptionSpec> &opts) const {
  PETFwdModel::GetOptions(opts);
  for (int i = 0; OPTIONS[i].name != ""; i++) {
    opts.push_back(OPTIONS[i]);
  }
}

void PET_2TCM_IR_FwdModel::Initialize(FabberRunData &rundata) {
  PETFwdModel::Initialize(rundata);

  // Initial values of model specific parameters
  m_init_K1 = rundata.GetDoubleDefault("init-K1", 0.05);
  m_init_Ki = rundata.GetDoubleDefault("init-Ki", 0.05);
  m_init_k_sum = rundata.GetDoubleDefault("init-k-sum", 0.05);

  m_ca = rundata.GetDoubleDefault("ca", 0);
  m_lc = rundata.GetDoubleDefault("lc", 0.81);
}

void PET_2TCM_IR_FwdModel::GetParameterDefaults(
    std::vector<Parameter> &params) const {

  // generic model parameters
  PETFwdModel::GetParameterDefaults(params);
  int p = params.size();

  // specific parameters
  params.push_back(Parameter(p++, "K1", DistParams(m_init_K1, 100),
                             DistParams(m_init_K1, 100), PRIOR_NORMAL,
                             TRANSFORM_LOG()));
  params.push_back(Parameter(p++, "Ki", DistParams(m_init_Ki, 100),
                             DistParams(m_init_Ki, 100), PRIOR_NORMAL,
                             TRANSFORM_LOG()));
  params.push_back(Parameter(p++, "k_sum", DistParams(m_init_k_sum, 100),
                             DistParams(m_init_k_sum, 100), PRIOR_NORMAL,
                             TRANSFORM_LOG()));
                             
}

void PET_2TCM_IR_FwdModel::EvaluateModel(const ColumnVector &params,
                                         ColumnVector &result,
                                         const std::string &key) const {
  if (key == "") {
    Evaluate(params, result);
  } else {
    result.ReSize(1);
    result(1) = params(3) * 6000 / m_density * m_ca / 18.0156 / m_lc;
  }
}

void PET_2TCM_IR_FwdModel::Evaluate(const ColumnVector &params,
                                    ColumnVector &result) const {
  // Parameters that are inferred - extract and give sensible names
  int p = 1;

  double vB = params(p++);
  double K1 = params(p++);
  double Ki = params(p++);
  double k_sum = params(p++);

  ColumnVector exp_results = MISCMATHS::exp(-k_sum * m_kernel_time);

  ColumnVector convolution_result_1 = m_c_mat * exp_results;
  ColumnVector convolution_result_2 = m_c_mat * (1 - exp_results);

  ColumnVector c_1 = K1 * convolution_result_1;
  ColumnVector c_2 = Ki * convolution_result_2;

  result = (1 - vB) * (c_1 + c_2) + vB * m_aif_pet;

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

void PET_2TCM_IR_FwdModel::GetOutputs(std::vector<std::string> &outputs) const {
  if (m_ca != 0) {
    outputs.push_back("CMRglc");
  }
}

FwdModel *PET_2TCM_IR_FwdModel::NewInstance() {
  return new PET_2TCM_IR_FwdModel();
}
