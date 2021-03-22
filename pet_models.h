/**
 * pet_models.h
 * 
 * Header file for differnet pharmacokinetic models for PET
 * Moss Zhao - Center for Advanced Functional Neuroimaging (CAFN), Stanford University

 */

#pragma once

#ifdef _WIN32
#ifdef fabber_pet_EXPORTS
#define FABBER_PET_API __declspec(dllexport)
#else
#define FABBER_PET_API __declspec(dllimport)
#endif
#define CALL __stdcall
#else
#define FABBER_PET_API
#define CALL
#endif

#include "fabber_core/fwdmodel.h"

extern "C" {
FABBER_PET_API int CALL get_num_models();
FABBER_PET_API const char *CALL get_model_name(int index);
FABBER_PET_API NewInstanceFptr CALL get_new_instance_func(const char *name);
}
