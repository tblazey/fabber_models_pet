/**
 * fwdmodel_pet.cc
 * 
 * Implementation of one tissue compartment model for PET
 * Moss Zhao - Center for Advanced Functional Neuroimaging (CAFN), Stanford University

 */

#include "fwdmodel_pet.h"

#include <fabber_core/easylog.h>
#include <fabber_core/priors.h>

#include <miscmaths/miscprob.h>
#include <newimage/newimageall.h>

#include <newmatio.h>

#include <iostream>
#include <stdexcept>

using namespace NEWMAT;

static OptionSpec OPTIONS[] = {
    { "K1", OPT_FLOAT, "K1 transfer rate", OPT_REQ, "" },
    { "k2", OPT_FLOAT, "k2 transfer rate", OPT_REQ, "" },
    { "aif", OPT_STR, "Source of AIF function: signal=User-supplied vascular signal", OPT_REQ, "none"},
    { "aif-data", OPT_MATRIX,
        "File containing single-column ASCII data defining the AIF. For aif=signal, this is the vascular signal curve.",
        OPT_NONREQ, "none" },
    { "time-data", OPT_MATRIX,
        "File containing single-column ASCII data defining the timing information.",
        OPT_NONREQ, "none" },
    /*    
    { "delt", OPT_FLOAT, "Time resolution between volumes, in minutes", OPT_REQ, "" },
    { "fa", OPT_FLOAT, "Flip angle in degrees.", OPT_REQ, "" },
    { "tr", OPT_FLOAT, "Repetition time (TR) In seconds.", OPT_REQ, "" },
    { "r1", OPT_FLOAT, "Relaxivity of contrast agent, In s^-1 mM^-1.", OPT_REQ, "" },

    { "t10", OPT_FLOAT, "Baseline T1 value in seconds. May be inferred.", OPT_NONREQ, "1" },
    { "sig0", OPT_FLOAT, "Baseline signal. This value is ignored if sig0 is inferred.", OPT_NONREQ, "1" },
    { "delay", OPT_FLOAT, "Injection time (or delay time when using measured AIF) in minutes. May be inferred.", OPT_NONREQ, "0" },
    { "infer-t10", OPT_BOOL, "Infer t10 value", OPT_NONREQ, "" },
    { "infer-sig0", OPT_BOOL, "Infer baseline signal", OPT_NONREQ, "" },
    { "infer-delay", OPT_BOOL, "Infer the delay parameter", OPT_NONREQ, "" },

    { "aif", OPT_STR, "Source of AIF function: signal=User-supplied vascular signal", OPT_REQ, "none"},
    */
    { "auto-init-delay", OPT_BOOL, "Automatically initialize posterior value of delay parameter", OPT_NONREQ, "" },
    { "" },
};

void PETFwdModel::GetOptions(vector<OptionSpec> &opts) const
{
    for (int i = 0; OPTIONS[i].name != ""; i++)
    {
        opts.push_back(OPTIONS[i]);
    }
}

string PETFwdModel::ModelVersion() const
{
    string version = "Fabber PET models: ";
#ifdef GIT_SHA1
    version += string(" Revision ") + GIT_SHA1;
#endif
#ifdef GIT_DATE
    version += string(" Last commit ") + GIT_DATE;
#endif
    return version;
}

void PETFwdModel::Initialize(FabberRunData &rundata)
{
    // Mandatory parameters
    // AIF
    m_aif_type = rundata.GetString("aif");
    if (m_aif_type == "signal")
    {   
        // Read in AIF signal from text file
        m_aif = read_ascii_matrix(rundata.GetString("aif-data"));
        m_time = read_ascii_matrix(rundata.GetString("time-data"));
        //cout << "Your AIf is..." << endl;
        //cout  << "C: " << m_aif.t() << endl;
        //cout << "Your time is..." << endl;
        //cout  << "C: " << m_time.t() << endl;
    }
    else
    {
        throw InvalidOptionValue("aif", m_aif_type, "Must be a signel file");
    }

    // Automatically initialise delay posterior
    m_auto_init_delay = rundata.ReadBool("auto-init-delay");
}

// This is to be implement in future releases
void PETFwdModel::GetParameterDefaults(std::vector<Parameter> &params) const
{
    // Add standard PET parameters to whatever's there
    /*
    m_sig0_idx = params.size();
    unsigned int p = m_sig0_idx;
    if (m_infer_sig0)
        params.push_back(Parameter(p++, "sig0", DistParams(m_sig0, 1e8), DistParams(m_sig0, 100), PRIOR_NORMAL, TRANSFORM_ABS()));
    if (m_infer_delay)
        params.push_back(Parameter(p++, "delay", DistParams(m_delay, 0.25), DistParams(m_delay, 0.25)));
    if (m_infer_t10)
        params.push_back(Parameter(p++, "t10", DistParams(m_t10, 1), DistParams(m_t10, 1), PRIOR_NORMAL, TRANSFORM_ABS()));
    */

}


static int fit_step(const ColumnVector &data)
{
    // MSC method for initializing delay parameter.
    // We fit a step function to the PET signal to estimate the delay. Otherwise the
    // delay can be problematic to infer.
    double best_ssq = -1;
    int step_pos = 0;
    for (int pos=1; pos<data.Nrows(); pos++) {
        double mean_left = data.Rows(1, pos).Sum() / pos;
        double mean_right = data.Rows(pos+1, data.Nrows()).Sum() / (data.Nrows() - pos);
        double ssq = 0;
        for (int t=1; t<=pos; t++) {
            ssq += (data(t) - mean_left)*(data(t) - mean_left);
        }
        for (int t=pos+1; t<=data.Nrows(); t++) {
            ssq += (data(t) - mean_right)*(data(t) - mean_right);
        }
        if ((ssq < best_ssq) || (best_ssq < 0)) {
            best_ssq = ssq;
            step_pos = pos;
        }
    }

    return step_pos-1;
}


void PETFwdModel::InitVoxelPosterior(MVNDist &posterior) const
{
    int delay_idx = m_sig0_idx;
    if (m_infer_sig0) {
        // Note that sig0 is the first PET signal
        double sig0_data = data(1);
        posterior.means(m_sig0_idx+1) = sig0_data;
        delay_idx++;
    }
    if (m_infer_delay && m_auto_init_delay) {
        posterior.means(delay_idx+1) = m_dt * fit_step(data);
    }
}
