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

#include <newmatio.h>

#include <iostream>
#include <stdexcept>

using namespace NEWMAT;

FactoryRegistration<FwdModelFactory, PET_1TCM_FwdModel> PET_1TCM_FwdModel::registration("pet_1TCM");

std::string PET_1TCM_FwdModel::GetDescription() const
{
    return "PET one tissue compartment model";
}

static OptionSpec OPTIONS[] = {
    { "K1", OPT_FLOAT, "Flow in min-1", OPT_NONREQ, "0.01" },
    { "k2", OPT_FLOAT, "Flow in min-1", OPT_NONREQ, "0.1" },
    { "vB", OPT_FLOAT, "Blood volume in 1", OPT_NONREQ, "0.05" },
    //{ "fp", OPT_FLOAT, "Flow in min-1", OPT_NONREQ, "0.5" },
    //{ "ps", OPT_FLOAT, "Permeability surface area product in min-1", OPT_NONREQ, "0.5" },
    //{ "vp", OPT_FLOAT, "Plasma volume in decimal between zero and one", OPT_NONREQ, "0.5" },
    //{ "ve", OPT_FLOAT, "Extracellular space volume in decimal between zero and one", OPT_NONREQ, "0.05" },
    //{ "conv-method", OPT_STR, "Method to compute convolution, trapezium, matrix or iterative. Default is iterative", OPT_REQ, "iterative"},
    { "" },
};

void PET_1TCM_FwdModel::GetOptions(vector<OptionSpec> &opts) const
{
    PETFwdModel::GetOptions(opts);
    for (int i = 0; OPTIONS[i].name != ""; i++)
    {
        opts.push_back(OPTIONS[i]);
    }
}

void PET_1TCM_FwdModel::Initialize(FabberRunData &rundata)
{
    PETFwdModel::Initialize(rundata);

    // Initial values of main parameters
    m_K1 = rundata.GetDoubleDefault("K1", 0.05);
    m_k2 = rundata.GetDoubleDefault("k2", 0.05);
    m_vB = rundata.GetDoubleDefault("vB", 0.03);

    // Initial values of main parameters
    //m_fp = rundata.GetDoubleDefault("fp", 0.5);
    //m_ps = rundata.GetDoubleDefault("ps", 0.5);
    //m_ve = rundata.GetDoubleDefault("ve", 0.5);
    //m_vp = rundata.GetDoubleDefault("vp", 0.05);

    // Other model options
    //m_conv_method = rundata.GetStringDefault("conv-method", "iterative");
}

void PET_1TCM_FwdModel::GetParameterDefaults(std::vector<Parameter> &params) const
{
    params.clear();

    // Basic model parameters
    int p = 0;
    // This set of parameters worked
    //params.push_back(Parameter(p++, "K1", DistParams(m_K1, 1e5), DistParams(m_K1, 100), PRIOR_NORMAL, TRANSFORM_LOG()));
    //params.push_back(Parameter(p++, "k2", DistParams(m_k2, 1e5), DistParams(m_k2, 100), PRIOR_NORMAL, TRANSFORM_LOG()));
    //params.push_back(Parameter(p++, "vB", DistParams(m_vB, 1e5), DistParams(m_vB, 100), PRIOR_NORMAL, TRANSFORM_LOG()));
    

    params.push_back(Parameter(p++, "K1", DistParams(m_K1, 100), DistParams(m_K1, 100), PRIOR_NORMAL, TRANSFORM_LOG()));
    params.push_back(Parameter(p++, "k2", DistParams(m_k2, 100), DistParams(m_k2, 100), PRIOR_NORMAL, TRANSFORM_LOG()));
    params.push_back(Parameter(p++, "vB", DistParams(m_vB, 100), DistParams(m_vB, 100), PRIOR_NORMAL, TRANSFORM_LOG()));
    //params.push_back(Parameter(p++, "fp", DistParams(m_fp, 1e5), DistParams(m_fp, 100), PRIOR_NORMAL, TRANSFORM_LOG()));
    //params.push_back(Parameter(p++, "ps", DistParams(m_ps, 1e5), DistParams(m_ps, 100), PRIOR_NORMAL, TRANSFORM_LOG()));
    //params.push_back(Parameter(p++, "ve", DistParams(m_ve, 1e5), DistParams(m_ve, 1), PRIOR_NORMAL, TRANSFORM_FRACTIONAL()));
    //params.push_back(Parameter(p++, "vp", DistParams(m_vp, 1), DistParams(m_vp, 1), PRIOR_NORMAL, TRANSFORM_FRACTIONAL()));
    
    // Standard DCE parameters
    PETFwdModel::GetParameterDefaults(params);
}

ColumnVector PET_1TCM_FwdModel::compute_convolution(const ColumnVector &vector_1, const ColumnVector &vector_2, const ColumnVector &vector_time) const
{

    ColumnVector convolution_result(data.Nrows());

    int vector_length_1 = vector_1.Nrows();
    int vector_length_2 = vector_2.Nrows();

    //cout << data.Nrows() << endl;
    //getchar();

    // Length of the full convolution results
    int temp_length = vector_length_1 + vector_length_2 - 1;

    ColumnVector convolution_full(temp_length);

    int time_0 = 0;

    // Compute convolution
    // Source code from here
    // https://toto-share.com/2011/11/cc-convolution-source-code/
    for (int i = 0; i < temp_length; i++) {
        int i1 = i;
        double temp = 0.0;

        for (int j = 0; j < vector_length_2; j++) {
            if(i1 >= 0 && i1 < vector_length_1) {
                temp = temp + vector_1(i1 + 1) * vector_2(j + 1);
            }
            i1 = i1 - 1;
            convolution_full(i + 1) = temp;
        }  
    }

    // Do the 'same' operation for the final convolution results
    /*
    int start_point = floor(temp_length / 2) + 1;
    for (int i = 1; i <= vector_length_1; i++) {
        //cout << start_point << endl;
        convolution_result(i) = convolution_full(start_point);
        start_point++;
    }
    */

    // Here we only need the first x elements of the convolution result
    for (int i = 1; i <= vector_length_1; i++) {
        /*
        if(i == 1) {
            convolution_result(i) = (vector_time(i) - time_0) * convolution_full(i);
        }
        else {
            convolution_result(i) = (vector_time(i) - vector_time(i - 1)) * convolution_full(i);
        }
        */
        
        convolution_result(i) = convolution_full(i);
    }  
    
    return convolution_result;
}

void PET_1TCM_FwdModel::Evaluate(const ColumnVector &params, ColumnVector &result) const
{
    // Parameters that are inferred - extract and give sensible names
    int p = 1;

    double K1 = params(p++);
    double k2 = params(p++);
    double vB = params(p++);

    /*
    double fp = params(p++);
    double ps = params(p++);
    double ve = params(p++);
    double vp = params(p++);

    // Standard DCE parameters - may be inferred
    double sig0 = m_sig0;
    double delay = m_delay;
    double t10 = m_t10;
    if (m_infer_sig0)
    {
        sig0 = params(p++);
    }
    if (m_infer_delay)
    {
        delay = params(p++);
    }
    if (m_infer_t10)
    {
        t10 = params(p++);
    }
    */



    // Tissue concentration results
    //ColumnVector concentration_tissue = compute_concentration(K1, k2);

    //double pre_delay = 30.2;


    //int vector_length_1 = m_time.Nrows();
    //ColumnVector time_new(vector_length_1);
    //for (int i = 1; i <= vector_length_1; i++) {
    //    time_new(i) = m_time(i) + pre_delay;
    //    //cout << time_new(i) << endl;
    //}


    ColumnVector exp_results = exp((-k2) * m_time);

    ColumnVector convolution_result = compute_convolution(exp_results, m_aif, m_time);
    
    // Converts concentration back to DCE signal
    /*
    result.ReSize(data.Nrows());
    for (int i = 1; i <= data.Nrows(); i++)
    {   
        result(i) = SignalFromConcentration(concentration_tissue(i), t10, sig0);
    }
    */

    //result = K1 * convolution_result;
    result = (1 - vB) * K1 * convolution_result + vB * m_aif;
    //result = K1 * exp_results;

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

FwdModel *PET_1TCM_FwdModel::NewInstance()
{
    return new PET_1TCM_FwdModel();
}
