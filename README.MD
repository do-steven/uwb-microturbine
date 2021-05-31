# Microturbine Simulation Capstone - University of Washington Bothell
Created by: Casey Barnett, Steven Do, and Jason Ginos

## Background
This repository was made for our senior capstone project for the 2020-21 UWB academic school year. It contains MATLAB, Simulink, and Excel files that are used in conjunction to simulate a home-installed microturbine generator. The main goal for this simulation was to provide a visual demonstration for future engineering students on how power generation system components interact with each other to produce usable output power. Specifically, the simulated microturbine operates on the Brayton cycle to produce input mechanical power, which is then fed into a synchronous generator to allow for the conversion to output electrical power. Concepts involved for this project include: simulation, thermal fluids, electric power, and mechanical systems design. See the following PDF file for prior documentation of this project, which includes more detail of its background and development. 

[Microturbine System Design - Senior Capstone Project (BME 494-495).pdf](https://github.com/do-steven/uwb-microturbine/files/6446747/Microturbine.System.Design.-.Senior.Capstone.Project.BME.494-495.pdf)

## How to Install
To successfully run this simulation, the following will be needed to be downloaded:

- Files contained in the repository
- MATLAB 2020 or newer (older versions have not been tested)
- Simulink
- All Simscape packages

**NOTE:** All of the files downloaded from this repository will need to be placed in the same folder, as the simulation will reference and pull necessary data from the files in the folder. Once this has been completed, run MATLAB and open the **Microturbine_Master_v3.slx** file. If all of the previous steps have been done properly, the simulation should run and allow for changes to be made at the users discretion. 

## Instructions for Use

## Description
Five different mechanical components layout the foundation for this simulation. These components include a compressor, regenerator, combustor, turbine, and synchronous generator. The following flowchart shows the process and general logic that the simulation follows. 

**INSERT FINALIZED VERSION OF PROCESS FLOWCHART HERE**

Descriptions for individuals components are listed below:

### 1. Compressor
*Relevant File(s):* Compressor.mlx | Properties of Air.xlsx

**Assumptions:**
- Changes in kinetic and potential energy are negligible
- Compressor is adiabatic (no heat transfer in/out)
- Input combustion by products are assumed to be air as an ideal gas
- Can be modeled either as isentropic or not

**Summary of Operation:**

The compressor subsystem uses rotational mechanical energy (provided by the Turbine) to compress air from the surrounding and pass it into the combustion chamber. Based on the required output power, the turbine subsystem calculates a total mass flow rate (air + fuel). This can be combined with the air/fuel mass ratio (which is known for combustion) in order to determine the required air mass flow rate.

### 2. Regenerator
*Relevant File(s):* Microturbine_Master_v3.slx

**Assumptions:**
- Regenerator is a counterflow heat exchanger
- Changes in kinetic and potential energy are negligible
- Maximum outlet temperature is equal to the turbine exhaust temperature
- Regenerator is well insulated

**Summary of Operation:**

A regenerator operates similarly to a heat exchanger and is connected between the compressor outlet and the exhaust of the system. The purpose for having a regenerator is to extract disposed heat from the exhaust and using it to further heat the fluid exiting the compressor. Having a regenerator is beneficial since it increases the thermal efficiency of the system by requiring less work to be done by the compressor, and it reduces fuel consumption since the working fluid can be heated from other sources prior to combustion. It is important to note that regeneration is only effective when the exhaust temperature is greater than the compressor outlet temperature; otherwise, heat transfer is reversed and energy is inadvertently lost to the surrounding. To prevent this from occurring, gas-turbine systems with low pressure ratios will benefit more from this component than those with high pressure ratios. The implementation of the regenerator within the code requires the regenerator outlet temperature to be defined, which serves as a "design temperature." Ideally, this temperature should be as close to the exhaust temperature for optimal performance. Without a heat exchanger designed, the analysis of the regenerator can only assess the expected heat transfer through this component, which could then provide insight to the size of the heat exhcanger needed for the application.

### 3. Combustor
*Relevant File(s):* combustion.mlx | Properties of Carbon Dioxide.xlsx | Properties of Nitrogen.xlsx | Properties of Oxygen.xlsx | Properties of Water Vapor.xlsx

**Assumptions:**
- Changes in kinetic and potential energy are negligible
- Combustion is complete
- Air is dry
- Combustion products are assumed to be air
- Only excess or stoichiometric air is considered

**Summary of Operation:**

The combustor is where fuel and air are reacted to extract chemical energy to do mechanical work downstream inside the turbine. A key metric when analyzing combustion is the adiabatic flame temperature. This temperature is useful since it is the maximum temperature that combustion occurs at given the composition of the reactants. Reactions occurring at the adiabatic flame temperature lose zero heat to the surrounding; therefore, this temperature is also the most efficient temperature for combustion to take place since all of the energy is used internally to do work downstream. The implementation of this code is focused toward solving for the adiabatic flame temperature. This is done by balancing the chemical equation for the reactants and quantities involved, calculating the enthalpy of the reactants, and then iteratively calculating the enthalpy of the products until the enthalpies are equal. The corresponding temperature is the adiabatic flame temperature.

### 4. Turbine
*Relevant File(s):* Turbine.mlx

**Assumptions:**
- Changes in kinetic and potential energy are negligible
- Turbine is adiabatic (no heat transfer in/out)
- Input combustion by products are assumed to be air as an ideal gas

**Summary of Operation:**

The turbine subsystem converts heat/pressure energy into rotational mechanical energy. This is done by funneling high pressure/temperature gases from the combustion chamber through a series of expanding rotor vanes. As the high temperature/pressure gases pass through the turbine stages, they expand and cool while transferring their heat/pressure energy into rotational mechanical energy (causing the rotors to spin). The spinning rotors are connected to a shaft which is in turn connected to the generator so that rotational energy can be transferred. 

### 5. Synchronous Generator
*Relevant File(s):* Microturbine_Master_v3.slx

**Assumptions:**
- 60Hz operating frequency (U.S. standard)
- Generator will run on three phase power at either 120/220V to accommodate small/large at-home appliances and devices
- Synchronous machine object accounts for stator and copper losses (mechanical to electrical conversion loss)
- Load will be wye connected

**Summary of Operation:**

The synchronous generator portion of the simulation operates in a reverse fashion (i.e. the simulation outputs the input power required to accommodate the loading conditions.) It is composed of three main components which are in the form of synchronous machine, loading, and excitation/voltage regulation objects. Essentially, the synchronous machine provides the necessary electrical power to the load based on its requirements, while the voltage regulation system maintains a steady voltage over time. The synchronous machine object is initialized through the "powergui", which solidifies the system's parameters and solving method that allows the simulation to be run. Once initialized, the simulation can be run using the green "start" button at the top of the Simulink window. A plot of the voltage vs. time of each phase can then be viewed through its respective scope. Additional plots of other parameters (output power, rotor speed, etc.) vs. time can be chosen through customization of the bus bar, which are then fed into the connecting scopes. The required input power to the generator gained from the simulation is sent to the *Microturbine.mlx* MATLAB script as a constant, which is then used for calculations involved with the turbine. 

## Further Discussion

Although we were able to create a fully functioning simulation of this power cycle, there is still much more work that can be done. First, this simulation could always be further refined to make the overall system more realistic by adding additional constraints and eliminating simplifying assumptions. However, this would require a much deeper understanding of the in depth thermodynamic and electric power concepts, likely with the help of experienced mentors/professionals. Once further refinement has been achieved, then individual component design would take place. The five different system components can be separately designed by different future capstone teams, using appropriate communicational skills especially when working with component integration. Concepts of machine design, finite element analysis, computational fluid dynamics, combustion analysis, control systems design, etc. would be necessary to build a full fledged microturbine generator that can supply the desired output power through the use of renewable energy. The overarching goal of this power generation system would be to have it be capable of providing back up electrical power to the user, which allows them to sustain comfort during a power outage. 
