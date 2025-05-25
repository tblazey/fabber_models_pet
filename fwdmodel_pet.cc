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
#include <armawrap/newmat.h>

#include <iostream>
#include <fstream>
#include <stdexcept>

using namespace std;
using namespace NEWMAT;
using MISCMATHS::read_ascii_matrix;

static OptionSpec OPTIONS[] = {
    { "aif-data", OPT_MATRIX,
        "File containing single-column ASCII data defining the AIF.",
        OPT_REQ, "" },
    { "pet-time-data", OPT_MATRIX,
        "File containing single-column ASCII data defining the timing information for PET data.",
        OPT_REQ, "" },
    { "aif-time-data", OPT_MATRIX,
        "File containing single-column ASCII data defining the aif timing information.",
        OPT_NONREQ, "" },
    { "init-vB", OPT_FLOAT, "Initial volume of blood (between 0 and 1)", OPT_NONREQ, "0.03" },
    { "density", OPT_FLOAT, "Density of brain (g/mL)", OPT_NONREQ, "1.05" },
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


Matrix PETFwdModel::convolve_matrix(const ColumnVector &kernel) const
{
    int n = kernel.Nrows();
    Matrix c_mat = Matrix(n, n);
    
    for (int i = 0; i < n; i++){
        for (int j = 0; j < n; j++){
            if (j < i + 1 && i - j < n ) {
                c_mat(i + 1, j + 1) = kernel(i - j + 1);
            } else{
                c_mat(i + 1, j + 1) = 0;
            }
        }       
    }
    return c_mat;
}

Matrix PETFwdModel::interp_matrix(const ColumnVector &x, const ColumnVector &x_p) const
{
    int n_x = x.Nrows();
    int n_p = x_p.Nrows();
    int r;
    int s;
    
    double mu;
    double x_one;
    double x_zero;
    
    Matrix i_mat = Matrix(n_p, n_x);
    
    // Loop through points to interpolate
    for (int i = 0; i < n_p; i++){
    
        // Loop through known points
        for (int j = 0; j < n_x; j++){
            
            // See if we found a boundary point
            if (x_p(i + 1) - x(j + 1) <= 0 || j + 1 == n_x){
                
                // Set boundry indicies
                if (j == 0){
                    r = j;
                    s = j + 1;
                } else{
                    r = j - 1;
                    s = j;
                }
                
                // Get boundry values
                x_zero = x(r + 1);
                x_one = x(s + 1);

                // Compute interpolation weight
                mu = (x_p(i + 1) - x_zero) / (x_one - x_zero);

                // Fill in interpolation matrix
                i_mat(i + 1, r + 1) = 1 - mu;
                i_mat(i + 1, s + 1) = mu;

                // Stop looping through reference points
                break;
            }
        }
    }
    
    return i_mat;
}

void PETFwdModel::Initialize(FabberRunData &rundata)
{

    // Initial values of main parameters
    m_init_vB = rundata.GetDoubleDefault("init-vB", 0.03);
    m_density = rundata.GetDoubleDefault("density", 1.05);

    // Read in AIF signal from text file
    ColumnVector aif = read_ascii_matrix(rundata.GetString("aif-data"));
    ColumnVector pet_time = read_ascii_matrix(rundata.GetString("pet-time-data"));

    // Load in aif time vector
    string aif_time_path = rundata.GetStringDefault("aif-time-data", "");
    if ( aif_time_path != "" ){
        Matrix aif_time = read_ascii_matrix(aif_time_path);
        
        // Figure out smallest sampling time
        double dt = aif_time(2) - aif_time(1);
        double dt_new;
        for (int i = 2; i < aif_time.Nrows(); i++){
            dt_new = aif_time(i + 1) - aif_time(i);
            if (dt_new < dt){
                dt = dt_new;
            }
        }
        
        // Make new array of timepoints
        double aif_min = aif_time.Minimum();
        double aif_max = aif_time.Maximum();
        int n_i = ceil((aif_max - aif_min) / dt);
        ColumnVector aif_time_i(n_i);
        for (int i = 1; i <= n_i; i++){
            aif_time_i(i) = aif_min + (i - 1) * dt;
        }
        m_kernel_time = aif_time_i - aif_min;
        
        // Interpolate aif to even sampling and to pet sampling
        ColumnVector aif_i = interp_matrix(aif_time, aif_time_i) * aif;
        m_aif_pet = interp_matrix(aif_time, pet_time) * aif;
        
        // Get matrix to interpolate + convolve
        m_c_mat = interp_matrix(aif_time_i, pet_time) * convolve_matrix(aif_i) * dt;
        
    } else{
        m_c_mat = convolve_matrix(aif) * (pet_time(2) - pet_time(1));
        m_kernel_time = pet_time - pet_time.Minimum();
        m_aif_pet = aif;
    }
       
}

void PETFwdModel::GetParameterDefaults(std::vector<Parameter> &params) const
{
    params.clear();

    // Basic model parameters
    int p = 0;
    params.push_back(Parameter(p++, "vB", DistParams(m_init_vB, 10), DistParams(m_init_vB, 10), PRIOR_NORMAL, TRANSFORM_LOG()));
}

void PETFwdModel::InitVoxelPosterior(MVNDist &posterior) const
{
}
