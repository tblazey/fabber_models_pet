/**
 * pet_models.cc
 * 
 * Implementation of differnet pharmacokinetic models for PET
 * Moss Zhao - Center for Advanced Functional Neuroimaging (CAFN), Stanford University

 */


#include "fabber_core/fwdmodel.h"

#include <algorithm>

#include "fwdmodel_pet_1TCM.h" // One Tissue Compartment Model

#include "fwdmodel_pet_2TCM.h" // reversible Two Compartment Model

#include "fwdmodel_pet_2TCM_IR.h" // Irreversible Two Compartment Model

#include "pet_models.h"

extern "C" {

// Number of available models
int CALL get_num_models()
{
    return 1;
}

const char *CALL get_model_name(int index)
{
    switch (index)
    {
    case 0:
        return "pet_1TCM";
        break;
    case 1:
        return "pet_2TCM";
        break;
    case 2
        return "pet_2TCM_IR";
        break;
    default:
        return NULL;
    }
}

// Create a new instance of the model
NewInstanceFptr CALL get_new_instance_func(const char *name)
{
    if (string(name) == "pet_1TCM")
    {
        return PET_1TCM_FwdModel::NewInstance;
    }
    else if (string(name) == "pet_2TCM")
    {
        return PET_2TCM_FwdModel::NewInstance;
    }
    else if (string(name) == "pet_2TCM_IR")
    {
        return PET_2TCM_IR_FwdModel::NewInstance;
    }
    else
    {
        return NULL;
    }
}
}
