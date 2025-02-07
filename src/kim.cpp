/*
 * Kim.cpp
 *
 *  Created on: Nov 18, 2019
 *      Author: Nikhil
 */

#include <iomanip>
#include <iostream>
#include "kim.h"
#include "typedef.h"

Kim::Kim(const std::string modelname) : modelname(modelname),kim_ptr(nullptr),computeArguments(nullptr)
{
	std::string message= "Connecting to model: " + modelname;
	MY_HEADING(message.c_str());

	int error, requestedUnitsAccepted;
	error = KIM::Model::Create(KIM::NUMBERING::zeroBased,
							   KIM::LENGTH_UNIT::A,
							   KIM::ENERGY_UNIT::eV,
							   KIM::CHARGE_UNIT::e,
							   KIM::TEMPERATURE_UNIT::K,
							   KIM::TIME_UNIT::ps,
							   modelname,
							   &requestedUnitsAccepted,
							   &kim_ptr);
	if (error) { MY_ERROR("KIM::Model::Create()"); }

	// Check for compatibility of units with the model's units
	if (!requestedUnitsAccepted) { MY_ERROR("Must Adapt to model units"); }

	error = kim_ptr->ComputeArgumentsCreate(&computeArguments);
	if (error) { MY_ERROR("Unable to create a ComputeArguments object."); }

	kim_ptr->GetInfluenceDistance(&influenceDistance);
	std::cout << "Influence distance = " << influenceDistance << std::endl;
}

Kim::~Kim()
{
    int error = kim_ptr->ComputeArgumentsDestroy(&computeArguments);
    if (error) { MY_ERROR("Unable to destroy compute arguments"); }
	/* call model destroy */
	KIM::Model::Destroy(&kim_ptr);
}

const double* Kim::getCutoffs() const
{
	const double* cutoffs;
	const int* modelWillNotRequestNeighborsOfNoncontributingParticles;
	int numberOfNeighborLists;
	kim_ptr->GetNeighborListPointers(&numberOfNeighborLists,
									 &cutoffs,
									 &modelWillNotRequestNeighborsOfNoncontributingParticles);
	return cutoffs;
}

int Kim::getNumberOfNeighborLists()  const
{
	const double* cutoffs;
	const int* modelWillNotRequestNeighborsOfNoncontributingParticles;
	int numberOfNeighborLists;
	kim_ptr->GetNeighborListPointers(&numberOfNeighborLists,
									 &cutoffs,
									 &modelWillNotRequestNeighborsOfNoncontributingParticles);
	return numberOfNeighborLists;

}

void Kim::queryModel()
{
	int error;

	// Check that we know about all required routines
	int numberOfModelRoutineNames;
	KIM::MODEL_ROUTINE_NAME::GetNumberOfModelRoutineNames(
	  &numberOfModelRoutineNames);
	for (int i = 0; i < numberOfModelRoutineNames; ++i)
	{
		KIM::ModelRoutineName modelRoutineName;
		int error
			= KIM::MODEL_ROUTINE_NAME::GetModelRoutineName(i, &modelRoutineName);
		if (error) { MY_ERROR("Unable to get ModelRoutineName."); }
		int present;
		int required;
		error = kim_ptr->IsRoutinePresent(
			modelRoutineName, &present, &required);
		if (error) { MY_ERROR("Unable to get routine present/required."); }

		std::cout << "Model routine name \"" << modelRoutineName.ToString()
				  << "\" has present = " << present
				  << " and required = " << required << "." << std::endl;

		if ((present == true) && (required == true))
		{
			using namespace KIM::MODEL_ROUTINE_NAME;
			if (!((modelRoutineName == Create)
				|| (modelRoutineName == ComputeArgumentsCreate)
				|| (modelRoutineName == Compute) || (modelRoutineName == Refresh)
				|| (modelRoutineName == ComputeArgumentsDestroy)
				|| (modelRoutineName == Destroy)))
			{
				MY_ERROR("Unknown Routine \"" + modelRoutineName.ToString()
						 + "\" is required by model.");
			}
		}
	}

	// print model units
	KIM::LengthUnit lengthUnit;
	KIM::EnergyUnit energyUnit;
	KIM::ChargeUnit chargeUnit;
	KIM::TemperatureUnit temperatureUnit;
	KIM::TimeUnit timeUnit;

	kim_ptr->GetUnits(
	    &lengthUnit, &energyUnit, &chargeUnit, &temperatureUnit, &timeUnit);

	std::cout << "LengthUnit is \"" << lengthUnit.ToString() << "\"" << std::endl
	          << "EnergyUnit is \"" << energyUnit.ToString() << "\"" << std::endl
	          << "ChargeUnit is \"" << chargeUnit.ToString() << "\"" << std::endl
	          << "TemperatureUnit is \"" << temperatureUnit.ToString() << "\""
	          << std::endl
	          << "TimeUnit is \"" << timeUnit.ToString() << "\"" << std::endl;

	int number_of_neighbor_lists;
	const double* cutoff_base;

	// print information about neighbor lists
	int const * modelWillNotRequestNeighborsOfNoncontributingParticles;
	kim_ptr->GetNeighborListPointers(&number_of_neighbor_lists,
									  &cutoff_base,
									  &modelWillNotRequestNeighborsOfNoncontributingParticles);
	std::cout << "Model has numberOfNeighborLists : " << number_of_neighbor_lists
	          << std::endl;
	for (int i = 0; i < number_of_neighbor_lists; ++i)
	{
	  std::cout << "\t"
	            << "Neighbor list " << i << " has cutoff "
	            << cutoff_base[i]
	            << " with "
	               "modelWillNotRequestNeighborsOfNoncontributingParticles "
	            << modelWillNotRequestNeighborsOfNoncontributingParticles[i]
	            << std::endl;
	}
//	ignoring hints from here on...
//	if (number_of_neighbor_lists != 1) MY_ERROR("too many neighbor lists");

	// check compute arguments
	int numberOfComputeArgumentNames;
	KIM::COMPUTE_ARGUMENT_NAME::GetNumberOfComputeArgumentNames(
	    &numberOfComputeArgumentNames);
	for (int i = 0; i < numberOfComputeArgumentNames; ++i)
	{
	  KIM::ComputeArgumentName computeArgumentName;
	  KIM::SupportStatus supportStatus;
	  KIM::COMPUTE_ARGUMENT_NAME::GetComputeArgumentName(i, &computeArgumentName);
	  KIM::DataType dataType;
	  KIM::COMPUTE_ARGUMENT_NAME::GetComputeArgumentDataType(computeArgumentName,
	                                                         &dataType);
	  error = computeArguments->GetArgumentSupportStatus(computeArgumentName,
	                                                     &supportStatus);
	  if (error) MY_ERROR("unable to get ComputeArgument SupportStatus");

	  std::cout << "ComputeArgument Name \"" << computeArgumentName.ToString()
	            << "\""
	            << " is of type \"" << dataType.ToString() << "\""
	            << " and has supportStatus \"" << supportStatus.ToString() << "\""
	            << std::endl;

	  // can only handle energy and force as a required arg
	  if (supportStatus == KIM::SUPPORT_STATUS::required)
	  {
	    if ((computeArgumentName != KIM::COMPUTE_ARGUMENT_NAME::partialEnergy)
	        && (computeArgumentName != KIM::COMPUTE_ARGUMENT_NAME::partialForces))
	    { MY_ERROR("unsupported required ComputeArgument"); }
	  }

	  // must have energy and forces
	  if ((computeArgumentName == KIM::COMPUTE_ARGUMENT_NAME::partialEnergy)
	      || (computeArgumentName == KIM::COMPUTE_ARGUMENT_NAME::partialForces))
	  {
	    if (!((supportStatus == KIM::SUPPORT_STATUS::required)
	          || (supportStatus == KIM::SUPPORT_STATUS::optional)))
	    { MY_ERROR("energy or forces not available"); }
	  }
	}

	// check compute callbacks
	int numberOfComputeCallbackNames;
	KIM::COMPUTE_CALLBACK_NAME::GetNumberOfComputeCallbackNames(
	    &numberOfComputeCallbackNames);
	for (int i = 0; i < numberOfComputeCallbackNames; ++i)
	{
	  KIM::ComputeCallbackName computeCallbackName;
	  KIM::COMPUTE_CALLBACK_NAME::GetComputeCallbackName(i, &computeCallbackName);
	  KIM::SupportStatus supportStatus;
	  computeArguments->GetCallbackSupportStatus(computeCallbackName,
	                                             &supportStatus);

	  std::cout << "ComputeCallback Name \"" << computeCallbackName.ToString()
	            << "\""
	            << " has supportStatus \"" << supportStatus.ToString() << "\""
	            << std::endl;

	  // cannot handle any "required" callbacks
	  if (supportStatus == KIM::SUPPORT_STATUS::required)
	  { MY_ERROR("unsupported required ComputeCallback"); }
	}

	int numberOfParameters;
	kim_ptr->GetNumberOfParameters(&numberOfParameters);
	for (int i = 0; i < numberOfParameters; ++i)
	{
	  KIM::DataType dataType;
	  std::string const * strName;
	  std::string const * strDesc;
	  int extent;
	  kim_ptr->GetParameterMetadata(
	      i, &dataType, &extent, &strName, &strDesc);
	  std::cout << "Parameter No. " << i << " has" << std::endl
	            << " data type   : \"" << dataType.ToString() << "\"" << std::endl
	            << " extent      : " << extent << std::endl
	            << " name        : " << *strName << std::endl
	            << " description : " << *strDesc << std::endl;
	}

	// Check supported extensions, if any
	int present;
	error = kim_ptr->IsRoutinePresent(
	    KIM::MODEL_ROUTINE_NAME::Extension, &present, NULL);
	if (error) { MY_ERROR("Unable to get Extension present/required."); }
	if (present)
	{
	  KIM::SupportedExtensions supportedExtensions;
	  error = kim_ptr->Extension(KIM_SUPPORTED_EXTENSIONS_ID,
	                                       &supportedExtensions);
	  if (error) { MY_ERROR("Error returned from KIM::Model::Extension()."); }
	  std::cout << "Model Supports "
	            << supportedExtensions.numberOfSupportedExtensions
	            << " Extensions:" << std::endl;
	  for (int i = 0; i < supportedExtensions.numberOfSupportedExtensions; ++i)
	  {
	    std::cout << " spportedExtensionID[" << std::setw(2) << i << "] = \""
	              << supportedExtensions.supportedExtensionID[i] << "\" "
	              << "which has required = "
	              << supportedExtensions.supportedExtensionRequired[i] << "."
	              << std::endl;
	  }
	}
}

void Kim::broadcastToModel(const Configuration* config_ptr,
						   const VectorXi& particleContributing,
                           const MatrixXd* forces_ptr,
					       NeighList* nl_ptr,
					       KIM::Function* get_neigh_ptr,
						   InteratomicForces* bonds,
						   KIM::Function* processDEDr_ptr)
{
	int error;
	// check species
	this->speciesCode.resize(config_ptr->numberOfParticles);
	int isSpeciesSupported;
	for (int i=0; i<config_ptr->numberOfParticles; i++)
	{
		KIM::SpeciesName speciesNameObject(config_ptr->species[i]);
		error = kim_ptr->GetSpeciesSupportAndCode(speciesNameObject, &isSpeciesSupported, speciesCode.data()+i);
		if ((error) || (!isSpeciesSupported))
		  { MY_ERROR("Species " + config_ptr->species[i] + " of particle " + std::to_string(i+1) + " not supported"); }
   	}


	error = computeArguments->SetArgumentPointer(KIM::COMPUTE_ARGUMENT_NAME::numberOfParticles, &(config_ptr->numberOfParticles))
	        ||
			computeArguments->SetArgumentPointer(KIM::COMPUTE_ARGUMENT_NAME::particleSpeciesCodes,this->speciesCode.data())
	        ||
			computeArguments->SetArgumentPointer(KIM::COMPUTE_ARGUMENT_NAME::particleContributing,particleContributing.data())
	        ||
			computeArguments->SetArgumentPointer(KIM::COMPUTE_ARGUMENT_NAME::coordinates, config_ptr->coordinates.at(Current).data())
            ||
            (forces_ptr == nullptr ? 0 :
             computeArguments->SetArgumentPointer(KIM::COMPUTE_ARGUMENT_NAME::partialForces, (*forces_ptr).data())
            );

	if (error) MY_ERROR("KIM_API_set_data");
	error = computeArguments->SetCallbackPointer(KIM::COMPUTE_CALLBACK_NAME::GetNeighborList,
												 KIM::LANGUAGE_NAME::cpp, get_neigh_ptr, nl_ptr)
			||

            (processDEDr_ptr== nullptr ? 0 :
                computeArguments->SetCallbackPointer(KIM::COMPUTE_CALLBACK_NAME::ProcessDEDrTerm,
                                                     KIM::LANGUAGE_NAME::cpp, processDEDr_ptr, bonds));
	if (error) MY_ERROR("set_call_back");

}

void Kim::compute()
{
	int error= kim_ptr->Compute(computeArguments);
    //if (error) MY_ERROR("compute");
    if (error) throw(std::runtime_error("compute"));
}

