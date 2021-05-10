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

**NOTE:** All of the files downloaded from this repository will need to be placed in the same folder, as the simulation will reference and pull necessary data from the files in the folder. Once this has been completed, run MATLAB and open the **Microturbine_Master_v2.slx** file. If all of the previous steps have been done properly, the simulation should run and allow for changes to be made at the users discretion. 

## Instructions for Use

## Description
Five different mechanical components layout the foundation for this simulation. These components include a compressor, regenerator, combustor, turbine, and synchronous generator. The following flowchart shows the process and general logic that the simulation follows. 

**INSERT FINALIZED VERSION OF PROCESS FLOWCHART HERE**

Descriptions for individuals components are listed below:

### 1. Compressor
*Relevent File(s):* Compressor.mlx | Properties of Air.xlsx

### 2. Regenerator
*Relevent File(s):* Microturbine_Master_v2.slx

### 3. Combustor
*Relevent File(s):* combustion.mlx | Properties of Carbon Dioxide.xlsx | Properties of Nitrogen.xlsx | Properties of Oxygen.xlsx | Properties of Water Vapor.xlsx

### 4. Turbine
*Relevent File(s):* Turbine.mlx 

### 5. Synchronous Generator
*Relevent File(s):* Microturbine_Master_v2.slx

  
