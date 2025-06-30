/*
 * calculateStress.h
 *
 *  Created on: Nov 7, 2019
 *      Author: Nikhil
 */

#ifndef CALCULATESTRESS_H_
#define CALCULATESTRESS_H_

int process_DEDr(const void* dataObject, const double de, const double r, const double* const dx, const int i, const int j);
#include "calculateStress.cpp"

/*!
 * \example{lineno} crack/main.cpp
 * Compares Hardy stress computation using
 * projected forces vs. process_dEdr with multiple OpenKIM models.
 * The system is a single crystal Si with a crack under mode 1 loading
 *
 *  - **Projected forces**: Interatomic forces derived from forces.
 *  - **process_dEdr**: Interatomic forces derived from analytical derivatives of
 *                      energy w.r.t distances, provided by KIM models (if available).
 *
 * For each model:
 * - We attempt to compute stress using **projected forces** first.
 * - Subsequently, we compute stress using **process_dEdr** if the model supports it.
 * - Both approaches use the `MethodSphere` class to compute Hardy stress on a spatial grid.
 *
 *
 * -# **Query KIM models**: For each model in the list, extract its lattice constant and
 * elastic constants. In this demo, the list contains only one model.
 * @snippet{lineno} crack/main.cpp QueryModels
 *
 * -# **Load atomic configuration**: Read a MDStressLab-style configuration (`config.data`).
 * @snippet{lineno} crack/main.cpp ReadConfiguration
 *
 * -# **Load grid**: Read grid points for stress evaluation from `grid_cauchy.data`.
 * @snippet{lineno} crack/main.cpp LoadGrid
 *
 * -# **Compute and compare Hardy stress**:
 *     - Start with stress computation using projected forces (always possible).
 *     - Try `process_dEdr` method, if available.
 *     - Write output stress fields for both methods (if available).
 * @snippet{lineno} crack/main.cpp ComputeStress
 *
 * #### Output:
 * For each KIM model, this example will generate:
 * - `project_hardy_<modelname>`: Stress computed using projected forces.
 * - `hardy_<modelname>`: Stress computed using process_dEdr.
 *
 * #### Visual comparison
 *
 *  * See @ref visualization_intro "Visualization Utilities" for contour plotting
 * <table>
 *   <tr>
 *     <td align="center">
 *       <img src="project_hardy_Tersoff_LAMMPS_Tersoff_1988T3_Si__MO_186459956893_004_xx.pdf" width="300px"/><br/>
 *       Projected force–based stress
 *     </td>
 *     <td align="center">
 *       <img src="hardy_Tersoff_LAMMPS_Tersoff_1988T3_Si__MO_186459956893_004_xx.pdf" width="300px"/><br/>
 *       process_dEdr–based stress
 *     </td>
 *   </tr>
 *   <tr>
 *     <td align="center">
 *       <img src="project_hardy_Tersoff_LAMMPS_Tersoff_1988T3_Si__MO_186459956893_004_yy.pdf" width="300px"/><br/>
 *       Projected force–based stress
 *     </td>
 *     <td align="center">
 *       <img src="hardy_Tersoff_LAMMPS_Tersoff_1988T3_Si__MO_186459956893_004_yy.pdf" width="300px"/><br/>
 *       process_dEdr–based stress
 *     </td>
 *   </tr>
 *   <tr>
 *     <td align="center">
 *       <img src="project_hardy_Tersoff_LAMMPS_Tersoff_1988T3_Si__MO_186459956893_004_xy.pdf" width="300px"/><br/>
 *       Projected force–based stress
 *     </td>
 *     <td align="center">
 *       <img src="hardy_Tersoff_LAMMPS_Tersoff_1988T3_Si__MO_186459956893_004_xy.pdf" width="300px"/><br/>
 *       process_dEdr–based stress
 *     </td>
 *   </tr>
 * </table>
 *
 * Full code:
 */


#endif /* CALCULATESTRESS_H_ */
