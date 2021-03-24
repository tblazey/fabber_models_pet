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

#include <newmat.h>

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
    void Evaluate(const NEWMAT::ColumnVector &params, NEWMAT::ColumnVector &result) const;
    
private:
    // Initial values of model parameters - always inferred
    double m_K1, m_k2, m_vB;
    //double m_fp, m_ps, m_ve, m_vp;

    // Convolution method
    //std::string m_conv_method;

    NEWMAT::ColumnVector compute_convolution(const NEWMAT::ColumnVector &vector_1, const NEWMAT::ColumnVector &vector_2, const NEWMAT::ColumnVector &vector_time) const;
    //NEWMAT::ColumnVector compute_convolution_matrix(const double delay, const double T, const double T_plus, const double T_minus) const;
    //NEWMAT::ColumnVector compute_convolution_iterative(const double delay, const double T_term) const;
    //NEWMAT::ColumnVector compute_convolution_trap(const double delay, const double T, const double T_plus, const double T_minus) const;
    
    /** Auto-register with forward model factory. */
    static FactoryRegistration<FwdModelFactory, PET_1TCM_FwdModel> registration;
};
